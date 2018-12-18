#include "MetaMorpheusGUI.CalibrateTaskWindow.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
#include "Util/ModTypeForTreeView.h"
#include "Util/ModTypeForLoc.h"
#include "Util/SearchModifications.h"
#include "Util/ModForTreeView.h"
#include "../EngineLayer/GlobalVariables.h"
#include "Util/GlobalGuiSettings.h"
#include "../EngineLayer/CommonParameters.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace TaskLayer;
using namespace Proteomics::ProteolyticDigestion;
using namespace MassSpectrometry;

namespace MetaMorpheusGUI
{

	CalibrateTaskWindow::CalibrateTaskWindow() : CalibrateTaskWindow(nullptr)
	{
	}

	CalibrateTaskWindow::CalibrateTaskWindow(CalibrationTask *myCalibrateTask)
	{
		InitializeComponent();
		PopulateChoices();
		CalibrationTask tempVar();
		setTheTask((myCalibrateTask != nullptr) ? myCalibrateTask : &tempVar);
		UpdateFieldsFromTask(getTheTask());

		if (myCalibrateTask == nullptr)
		{
			this->saveButton->Content = L"Add the Calibration Task";
		}
		SearchModifications::Timer->Tick += new EventHandler(TextChangeTimerHandler);
	}

	CalibrationTask *CalibrateTaskWindow::getTheTask() const
	{
		return privateTheTask;
	}

	void CalibrateTaskWindow::setTheTask(CalibrationTask *value)
	{
		privateTheTask = value;
	}

	void CalibrateTaskWindow::UpdateFieldsFromTask(CalibrationTask *task)
	{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		missedCleavagesTextBox->Text = task->getCommonParameters()->getDigestionParams()->MaxMissedCleavages == std::numeric_limits<int>::max() ? L"" : task->getCommonParameters()->getDigestionParams()->MaxMissedCleavages.ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MinPeptideLengthTextBox->Text = task->getCommonParameters()->getDigestionParams()->MinPeptideLength.ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MaxPeptideLengthTextBox->Text = task->getCommonParameters()->getDigestionParams()->MaxPeptideLength == std::numeric_limits<int>::max() ? L"" : task->getCommonParameters()->getDigestionParams()->MaxPeptideLength.ToString(CultureInfo::InvariantCulture);
		proteaseComboBox->SelectedItem = task->getCommonParameters()->getDigestionParams()->Protease;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		maxModificationIsoformsTextBox->Text = task->getCommonParameters()->getDigestionParams()->MaxModificationIsoforms.ToString(CultureInfo::InvariantCulture);
		initiatorMethionineBehaviorComboBox->SelectedIndex = static_cast<int>(task->getCommonParameters()->getDigestionParams()->InitiatorMethionineBehavior);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		DissociationTypeComboBox->SelectedItem = task->getCommonParameters()->getDissociationType()->ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		maxThreadsTextBox->Text = task->getCommonParameters()->getMaxThreadsToUsePerFile().ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MinVariantDepthTextBox->Text = task->getCommonParameters()->getMinVariantDepth().ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MaxHeterozygousVariantsTextBox->Text = task->getCommonParameters()->getMaxHeterozygousVariants().ToString(CultureInfo::InvariantCulture);

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		productMassToleranceTextBox->Text = task->getCommonParameters()->getProductMassTolerance()->Value->ToString(CultureInfo::InvariantCulture);
		productMassToleranceComboBox->SelectedIndex = dynamic_cast<AbsoluteTolerance*>(task->getCommonParameters()->getProductMassTolerance()) != nullptr ? 0 : 1;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		precursorMassToleranceTextBox->Text = task->getCommonParameters()->getPrecursorMassTolerance()->Value->ToString(CultureInfo::InvariantCulture);
		precursorMassToleranceComboBox->SelectedIndex = dynamic_cast<AbsoluteTolerance*>(task->getCommonParameters()->getPrecursorMassTolerance()) != nullptr ? 0 : 1;

		//writeIntermediateFilesCheckBox.IsChecked = task.CalibrationParameters.WriteIntermediateFiles;

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		minScoreAllowed->Text = task->getCommonParameters()->getScoreCutoff().ToString(CultureInfo::InvariantCulture);

		OutputFileNameTextBox->Text = task->getCommonParameters()->getTaskDescriptor();

		for (auto mod : task->getCommonParameters()->ListOfModsFixed)
		{
			auto theModType = FixedModTypeForTreeViewObservableCollection->FirstOrDefault([&] (std::any b)
			{
				b::DisplayName->Equals(mod->Item1);
			});
			if (theModType != nullptr)
			{
				auto theMod = theModType->Children->FirstOrDefault([&] (std::any b)
				{
					b::ModName->Equals(mod->Item2);
				});
				if (theMod != nullptr)
				{
					theMod->Use = true;
				}
				else
				{
					ModForTreeView tempVar(L"UNKNOWN MODIFICATION!", true, mod->Item2, true, theModType);
					theModType->Children->Add(&tempVar);
				}
			}
			else
			{
				theModType = new ModTypeForTreeView(mod->Item1, true);
				FixedModTypeForTreeViewObservableCollection->Add(theModType);
				ModForTreeView tempVar2(L"UNKNOWN MODIFICATION!", true, mod->Item2, true, theModType);
				theModType->Children->Add(&tempVar2);
			}
		}
		for (auto mod : task->getCommonParameters()->ListOfModsVariable)
		{
			auto theModType = VariableModTypeForTreeViewObservableCollection->FirstOrDefault([&] (std::any b)
			{
				b::DisplayName->Equals(mod->Item1);
			});
			if (theModType != nullptr)
			{
				auto theMod = theModType->Children->FirstOrDefault([&] (std::any b)
				{
					b::ModName->Equals(mod->Item2);
				});
				if (theMod != nullptr)
				{
					theMod->Use = true;
				}
				else
				{
					ModForTreeView tempVar3(L"UNKNOWN MODIFICATION!", true, mod->Item2, true, theModType);
					theModType->Children->Add(&tempVar3);
				}
			}
			else
			{
				theModType = new ModTypeForTreeView(mod->Item1, true);
				VariableModTypeForTreeViewObservableCollection->Add(theModType);
				ModForTreeView tempVar4(L"UNKNOWN MODIFICATION!", true, mod->Item2, true, theModType);
				theModType->Children->Add(&tempVar4);
			}
		}

		for (auto heh : LocalizeModTypeForTreeViewObservableCollection)
		{
			heh->setUse(std::make_optional(false));
		}
		for (auto ye : VariableModTypeForTreeViewObservableCollection)
		{
			ye->VerifyCheckState();
		}
		for (auto ye : FixedModTypeForTreeViewObservableCollection)
		{
			ye->VerifyCheckState();
		}
	}

	void CalibrateTaskWindow::PopulateChoices()
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
			DissociationTypeComboBox::Items->Add(dissassociationType.first);
		}

		productMassToleranceComboBox::Items->Add(L"Da");
		productMassToleranceComboBox::Items->Add(L"ppm");
		precursorMassToleranceComboBox::Items->Add(L"Da");
		precursorMassToleranceComboBox::Items->Add(L"ppm");

		for (auto hm : GlobalVariables::getAllModsKnown().GroupBy([&] (std::any b)
		{
			b::ModificationType;
		}))
		{
			auto theModType = new ModTypeForTreeView(hm::Key, false);
			FixedModTypeForTreeViewObservableCollection->Add(theModType);

			for (auto uah : hm)
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				ModForTreeView tempVar(uah->ToString(), false, uah->IdWithMotif, false, theModType);
				theModType->getChildren()->Add(&tempVar);
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete theModType' statement was not added since theModType was passed to a method or constructor. Handle memory management manually.
		}

		fixedModsTreeView->DataContext = FixedModTypeForTreeViewObservableCollection;

		for (auto hm : GlobalVariables::getAllModsKnown().GroupBy([&] (std::any b)
		{
			b::ModificationType;
		}))
		{
			auto theModType = new ModTypeForTreeView(hm::Key, false);
			VariableModTypeForTreeViewObservableCollection->Add(theModType);

			for (auto uah : hm)
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				ModForTreeView tempVar2(uah->ToString(), false, uah->IdWithMotif, false, theModType);
				theModType->getChildren()->Add(&tempVar2);
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete theModType' statement was not added since theModType was passed to a method or constructor. Handle memory management manually.
		}

		variableModsTreeView->DataContext = VariableModTypeForTreeViewObservableCollection;
	}

	void CalibrateTaskWindow::CancelButton_Click(std::any sender, RoutedEventArgs *e)
	{
		DialogResult = false;
	}

	void CalibrateTaskWindow::SaveButton_Click(std::any sender, RoutedEventArgs *e)
	{
		std::wstring fieldNotUsed = L"1";

		if (!GlobalGuiSettings::CheckTaskSettingsValidity(precursorMassToleranceTextBox->Text, productMassToleranceTextBox->Text, missedCleavagesTextBox->Text, maxModificationIsoformsTextBox->Text, MinPeptideLengthTextBox->Text, MaxPeptideLengthTextBox->Text, maxThreadsTextBox->Text, minScoreAllowed->Text, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed))
		{
			return;
		}

		Protease *protease = static_cast<Protease*>(proteaseComboBox::SelectedItem);
		int MaxMissedCleavages = missedCleavagesTextBox->Text->empty() ? std::numeric_limits<int>::max() : (std::stoi(missedCleavagesTextBox->Text, CultureInfo::InvariantCulture));
		int MinPeptideLength = std::stoi(MinPeptideLengthTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture);
		int MaxPeptideLength = MaxPeptideLengthTextBox->Text->empty() ? std::numeric_limits<int>::max() : (std::stoi(MaxPeptideLengthTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture));
		int MinVariantDepth = std::stoi(MinVariantDepthTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture);
		int MaxHeterozygousVariants = std::stoi(MaxHeterozygousVariantsTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture);
		int MaxModificationIsoforms = std::stoi(maxModificationIsoformsTextBox->Text, CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		DissociationType *dissociationType = GlobalVariables::getAllSupportedDissociationTypes()[DissociationTypeComboBox::SelectedItem.ToString()];

		DigestionParams *digestionParamsToSave = new DigestionParams(protease: protease->Name, maxMissedCleavages: MaxMissedCleavages, minPeptideLength: MinPeptideLength, maxPeptideLength: MaxPeptideLength, maxModificationIsoforms: MaxModificationIsoforms);

		auto listOfModsVariable = std::vector<(std::wstring, std::wstring)*>();
		for (auto heh : VariableModTypeForTreeViewObservableCollection)
		{
			listOfModsVariable.AddRange(heh->getChildren()->Where([&] (std::any b)
			{
				b::Use;
			})->Select([&] (std::any b)
			{
				(b::Parent->DisplayName, b::ModName);
			}));
		}

		if (!GlobalGuiSettings::VariableModCheck(listOfModsVariable))
		{
			delete digestionParamsToSave;
			return;
		}

		auto listOfModsFixed = std::vector<(std::wstring, std::wstring)*>();
		for (auto heh : FixedModTypeForTreeViewObservableCollection)
		{
			listOfModsFixed.AddRange(heh->getChildren()->Where([&] (std::any b)
			{
				b::Use;
			})->Select([&] (std::any b)
			{
				(b::Parent->DisplayName, b::ModName);
			}));
		}
		Tolerance *ProductMassTolerance;
		if (productMassToleranceComboBox->SelectedIndex == 0)
		{
			ProductMassTolerance = new AbsoluteTolerance(std::stod(productMassToleranceTextBox->Text, CultureInfo::InvariantCulture));
		}
		else
		{
			ProductMassTolerance = new PpmTolerance(std::stod(productMassToleranceTextBox->Text, CultureInfo::InvariantCulture));
		}

		Tolerance *PrecursorMassTolerance;
		if (precursorMassToleranceComboBox->SelectedIndex == 0)
		{
			PrecursorMassTolerance = new AbsoluteTolerance(std::stod(precursorMassToleranceTextBox->Text, CultureInfo::InvariantCulture));
		}
		else
		{
			PrecursorMassTolerance = new PpmTolerance(std::stod(precursorMassToleranceTextBox->Text, CultureInfo::InvariantCulture));
		}

		bool parseMaxThreadsPerFile = std::stoi(maxThreadsTextBox->Text, CultureInfo::InvariantCulture) <= Environment::ProcessorCount && std::stoi(maxThreadsTextBox->Text, CultureInfo::InvariantCulture) > 0;

		CommonParameters tempVar();
		CommonParameters *commonParamsToSave = new CommonParameters(OutputFileNameTextBox->Text != L"" ? OutputFileNameTextBox->Text: L"CalibrateTask", dissociationType, , , , , , , , std::stod(minScoreAllowed->Text, CultureInfo::InvariantCulture), , , , , , , ProductMassTolerance, PrecursorMassTolerance, , parseMaxThreadsPerFile ? std::stoi(maxThreadsTextBox->Text, CultureInfo::InvariantCulture): (&tempVar)->getMaxThreadsToUsePerFile(), digestionParamsToSave, listOfModsVariable, listOfModsFixed, , protease->Name != L"top-down", MaxHeterozygousVariants, MinVariantDepth);

		getTheTask()->setCommonParameters(commonParamsToSave);

		DialogResult = true;

//C# TO C++ CONVERTER TODO TASK: A 'delete commonParamsToSave' statement was not added since commonParamsToSave was assigned to another object. Handle memory management manually.
		delete PrecursorMassTolerance;
		delete ProductMassTolerance;
		delete digestionParamsToSave;
	}

	void CalibrateTaskWindow::CheckIfNumber(std::any sender, TextCompositionEventArgs *e)
	{
		e->Handled = !GlobalGuiSettings::CheckIsNumber(e->Text);
	}

	void CalibrateTaskWindow::KeyPressed(std::any sender, KeyEventArgs *e)
	{
		if (e->Key == Key->Return)
		{
			SaveButton_Click(sender, e);
		}
		else if (e->Key == Key->Escape)
		{
			CancelButton_Click(sender, e);
		}
	}

	void CalibrateTaskWindow::TextChanged_Fixed(std::any sender, TextChangedEventArgs *args)
	{
		SearchModifications::SetTimer();
		SearchModifications::FixedSearch = true;
	}

	void CalibrateTaskWindow::TextChanged_Var(std::any sender, TextChangedEventArgs *args)
	{
		SearchModifications::SetTimer();
		SearchModifications::VariableSearch = true;
	}

	void CalibrateTaskWindow::TextChangeTimerHandler(std::any sender, EventArgs *e)
	{
		if (SearchModifications::FixedSearch)
		{
			SearchModifications::FilterTree(SearchFixMod, fixedModsTreeView, FixedModTypeForTreeViewObservableCollection);
			SearchModifications::FixedSearch = false;
		}

		if (SearchModifications::VariableSearch)
		{
			SearchModifications::FilterTree(SearchVarMod, variableModsTreeView, VariableModTypeForTreeViewObservableCollection);
			SearchModifications::VariableSearch = false;
		}
	}
}
