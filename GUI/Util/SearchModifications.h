#pragma once

#include <string>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace MetaMorpheusGUI { class ModTypeForTreeView; }


namespace MetaMorpheusGUI
{
	class SearchModifications
	{
	public:
		static DispatcherTimer *Timer;
		static bool FixedSearch;
		static bool VariableSearch;
		static bool GptmdSearch;

		static void SetUpModSearchBoxes();

		// starts timer to keep track of user keystrokes
		static void SetTimer();

		// filters and expands tree according to user mod search
		static void FilterTree(TextBox *textbox, TreeView *tree, ObservableCollection<ModTypeForTreeView*> *collection);
	};
}
