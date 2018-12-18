#include "MetaMorpheusGUI.SearchTaskWindow.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "SearchTaskWindow.xaml.h"
#include "ForDisplayingInDataGrids/SearchModeForDataGrid.h"
#include "Util/ModTypeForTreeView.h"
#include "Util/ModTypeForLoc.h"
#include "Util/ModTypeForGrid.h"
#include "Util/SearchModifications.h"
#include "Util/GlobalGuiSettings.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace MassSpectrometry;

namespace MetaMorpheusGUI
{

	SearchTaskWindow::SearchTaskWindow() : SearchTaskWindow(nullptr)
	{
	}

	SearchTaskWindow::SearchTaskWindow(SearchTask *task) : DataContextForSearchTaskWindow(new DataContextForSearchTaskWindow { ExpanderTitle = std::wstring::Join(L", ", SearchModesForThisTask->Where([&] (std::any b)
	{
			b::Use;
	})->Select([&] (std::any b)
	{
			b->Name;
		})), AnalysisExpanderTitle = L"Some analysis properties...", SearchModeExpanderTitle = L"Some search properties..." })
		{
		InitializeComponent();
		SearchTask tempVar();
		setTheTask((task != nullptr) ? task : &tempVar);
		PopulateChoices();
		UpdateFieldsFromTask(getTheTask());

		if (task == nullptr)
		{
			this->saveButton->Content = L"Add the Search Task";
		}

		this->DataContext = DataContextForSearchTaskWindow;
		SearchModifications::Timer->Tick += new EventHandler(TextChangeTimerHandler);
	}

	SearchTask *SearchTaskWindow::getTheTask() const
	{
		return privateTheTask;
	}

	void SearchTaskWindow::setTheTask(SearchTask *value)
	{
		privateTheTask = value;
	}

	void SearchTaskWindow::CheckIfNumber(std::any sender, TextCompositionEventArgs *e)
	{
		e->Handled = !GlobalGuiSettings::CheckIsNumber(e->Text);
	}

	void SearchTaskWindow::Row_DoubleClick(std::any sender, MouseButtonEventArgs *e)
	{
		auto ye = dynamic_cast<DataGridCell*>(sender);
//C# TO C++ CONVERTER TODO TASK: C++ has no equivalent to C# pattern variables in 'is' expressions:
//ORIGINAL LINE: if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
		if (dynamic_cast<TextBlock*>(ye->Content) != nullptr hm && !hm->Text->empty())
		{
			System::Diagnostics::Process::Start(hm->Text);
		}
	}

	void SearchTaskWindow::PopulateChoices()
	{
		for (Protease *protease : ProteaseDictionary::Dictionary->Values)
		{
			proteaseComboBox::Items->Add(protease);
		}
		Protease *trypsin = ProteaseDictionary::Dictionary[L"trypsin"];
		proteaseComboBox->SelectedItem = trypsin;

		for (auto initiatior_methionine_behavior : Enum::GetNames(InitiatorMethionineBehavior::typeid))
		{
			initiatorMethionineBehaviorComboBox::Items->Add(initiatior_methionine_behavior);
		}

		for (auto dissassociationType : GlobalVariables::getAllSupportedDissociationTypes())
		{
			dissociationTypeComboBox::Items->Add(dissassociationType.first);
		}

		productMassToleranceComboBox::Items->Add(L"Da");
		productMassToleranceComboBox::Items->Add(L"ppm");

		precursorMassToleranceComboBox::Items->Add(L"Da");
		precursorMassToleranceComboBox::Items->Add(L"ppm");

		for (auto hm : GlobalVariables::getAllModsKnown().Where([&] (std::any b)
		{
			return b->ValidModification == true;
		}).GroupBy([&] (std::any b)
		{
			b::ModificationType;
		}))
		{
			ModSelectionGridItems->Add(theModType);
		}
	}
}
