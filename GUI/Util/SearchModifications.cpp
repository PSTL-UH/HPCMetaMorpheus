#include "SearchModifications.h"
#include "ModTypeForTreeView.h"


namespace MetaMorpheusGUI
{

DispatcherTimer *SearchModifications::Timer;
bool SearchModifications::FixedSearch = false;
bool SearchModifications::VariableSearch = false;
bool SearchModifications::GptmdSearch = false;

	void SearchModifications::SetUpModSearchBoxes()
	{
		Timer = new DispatcherTimer();
		Timer->Interval = TimeSpan::FromMilliseconds(300);
	}

	void SearchModifications::SetTimer()
	{
		// Reset the timer
		Timer->Stop();
		Timer->Start();
	}

	void SearchModifications::FilterTree(TextBox *textbox, TreeView *tree, ObservableCollection<ModTypeForTreeView*> *collection)
	{
		std::wstring key = textbox->Text->ToLower();
		if (key.empty())
		{
			tree->DataContext = collection; // shows full tree if nothing is searched
			return;
		}

		auto modTypesWithMatchingMods = collection->Where([&] (std::any p)
		{
			p::Children->Any([&] (std::any c)
			{
				c::ModName->ToLower()->Contains(key);
			});
		}); // parent of child mods that match key

		auto modsThatMatchSearchString = new ObservableCollection<ModTypeForTreeView*>(); // new collection containing expanded mod types that match key

		for (auto modType : modTypesWithMatchingMods)
		{
			auto textFilteredModType = new ModTypeForTreeView(modType->getDisplayName(), false);
			modsThatMatchSearchString->Add(textFilteredModType);
			textFilteredModType->setExpanded(true);
			textFilteredModType->setUse(modType->getUse());

			auto matchingChildren = modType->getChildren()->Where([&] (std::any p)
			{
				p::ModName->ToLower()->Contains(key);
			});
			for (auto mod : matchingChildren)
			{
				textFilteredModType->getChildren()->Add(mod);
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete textFilteredModType' statement was not added since textFilteredModType was passed to a method or constructor. Handle memory management manually.
		}

		tree->DataContext = modsThatMatchSearchString;

//C# TO C++ CONVERTER TODO TASK: A 'delete modsThatMatchSearchString' statement was not added since modsThatMatchSearchString was assigned to another object. Handle memory management manually.
	}
}
