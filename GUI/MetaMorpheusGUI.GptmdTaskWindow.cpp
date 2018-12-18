#include "MetaMorpheusGUI.GptmdTaskWindow.h"
#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
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

	GptmdTaskWindow::GptmdTaskWindow() : GptmdTaskWindow(nullptr)
	{
	}

	GptmdTaskWindow::GptmdTaskWindow(GptmdTask *myGPTMDtask)
	{
		InitializeComponent();
		PopulateChoices();

		GptmdTask tempVar();
		setTheTask((myGPTMDtask != nullptr) ? myGPTMDtask : &tempVar);
		UpdateFieldsFromTask(getTheTask());

		if (myGPTMDtask == nullptr)
		{
			this->saveButton->Content = L"Add the GPTMD Task";
		}

		SearchModifications::Timer->Tick += new EventHandler(TextChangeTimerHandler);
	}

	GptmdTask *GptmdTaskWindow::getTheTask() const
	{
		return privateTheTask;
	}

	void GptmdTaskWindow::setTheTask(GptmdTask *value)
	{
		privateTheTask = value;
	}

	void GptmdTaskWindow::Row_DoubleClick(std::any sender, MouseButtonEventArgs *e)
	{
		auto ye = dynamic_cast<DataGridCell*>(sender);
//C# TO C++ CONVERTER TODO TASK: C++ has no equivalent to C# pattern variables in 'is' expressions:
//ORIGINAL LINE: if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
		if (dynamic_cast<TextBlock*>(ye->Content) != nullptr hm && !hm->Text->empty())
		{
			System::Diagnostics::Process::Start(hm->Text);
		}
	}

	void GptmdTaskWindow::UpdateFieldsFromTask(GptmdTask *task)
	{
		useProvidedPrecursor->IsChecked = task->getCommonParameters()->getUseProvidedPrecursorInfo();
		deconvolutePrecursors->IsChecked = task->getCommonParameters()->getDoPrecursorDeconvolution();
		DeconvolutionMaxAssumedChargeStateTextBox->Text = std::to_wstring(task->getCommonParameters()->getDeconvolutionMaxAssumedChargeState());
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
		productMassToleranceTextBox->Text = task->getCommonParameters()->getProductMassTolerance()->Value->ToString(CultureInfo::InvariantCulture);
		productMassToleranceComboBox->SelectedIndex = dynamic_cast<AbsoluteTolerance*>(task->getCommonParameters()->getProductMassTolerance()) != nullptr ? 0 : 1;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		precursorMassToleranceTextBox->Text = task->getCommonParameters()->getPrecursorMassTolerance()->Value->ToString(CultureInfo::InvariantCulture);
		precursorMassToleranceComboBox->SelectedIndex = dynamic_cast<AbsoluteTolerance*>(task->getCommonParameters()->getPrecursorMassTolerance()) != nullptr ? 0 : 1;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		minScoreAllowed->Text = task->getCommonParameters()->getScoreCutoff().ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		maxThreadsTextBox->Text = task->getCommonParameters()->getMaxThreadsToUsePerFile().ToString(CultureInfo::InvariantCulture);
		addCompIonCheckBox->IsChecked = task->getCommonParameters()->getAddCompIons();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MinVariantDepthTextBox->Text = task->getCommonParameters()->getMinVariantDepth().ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MaxHeterozygousVariantsTextBox->Text = task->getCommonParameters()->getMaxHeterozygousVariants().ToString(CultureInfo::InvariantCulture);

		OutputFileNameTextBox->Text = task->getCommonParameters()->getTaskDescriptor();

		for (auto mod : task->getCommonParameters()->ListOfModsFixed)
		{
			auto theModType = fixedModTypeForTreeViewObservableCollection->FirstOrDefault([&] (std::any b)
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
				fixedModTypeForTreeViewObservableCollection->Add(theModType);
				ModForTreeView tempVar2(L"UNKNOWN MODIFICATION!", true, mod->Item2, true, theModType);
				theModType->Children->Add(&tempVar2);
			}
		}
		for (auto mod : task->getCommonParameters()->ListOfModsVariable)
		{
			auto theModType = variableModTypeForTreeViewObservableCollection->FirstOrDefault([&] (std::any b)
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
				variableModTypeForTreeViewObservableCollection->Add(theModType);
				ModForTreeView tempVar4(L"UNKNOWN MODIFICATION!", true, mod->Item2, true, theModType);
				theModType->Children->Add(&tempVar4);
			}
		}

		for (auto heh : localizeModTypeForTreeViewObservableCollection)
		{
			heh->setUse(std::make_optional(false));
		}

		for (auto mod : task->getGptmdParameters()->ListOfModsGptmd)
		{
			auto theModType = gptmdModTypeForTreeViewObservableCollection->FirstOrDefault([&] (std::any b)
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
					ModForTreeView tempVar5(L"UNKNOWN MODIFICATION!", true, mod->Item2, true, theModType);
					theModType->Children->Add(&tempVar5);
				}
			}
			else
			{
				theModType = new ModTypeForTreeView(mod->Item1, true);
				gptmdModTypeForTreeViewObservableCollection->Add(theModType);
				ModForTreeView tempVar6(L"UNKNOWN MODIFICATION!", true, mod->Item2, true, theModType);
				theModType->Children->Add(&tempVar6);
			}
		}

		for (auto ye : variableModTypeForTreeViewObservableCollection)
		{
			ye->VerifyCheckState();
		}
		for (auto ye : fixedModTypeForTreeViewObservableCollection)
		{
			ye->VerifyCheckState();
		}

		for (auto ye : gptmdModTypeForTreeViewObservableCollection)
		{
			ye->VerifyCheckState();
		}
	}

	void GptmdTaskWindow::PopulateChoices()
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
			fixedModTypeForTreeViewObservableCollection->Add(theModType);
			for (auto uah : hm)
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				ModForTreeView tempVar(uah->ToString(), false, uah->IdWithMotif, false, theModType);
				theModType->getChildren()->Add(&tempVar);
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete theModType' statement was not added since theModType was passed to a method or constructor. Handle memory management manually.
		}
		fixedModsTreeView->DataContext = fixedModTypeForTreeViewObservableCollection;
		for (auto hm : GlobalVariables::getAllModsKnown().GroupBy([&] (std::any b)
		{
			b::ModificationType;
		}))
		{
			auto theModType = new ModTypeForTreeView(hm::Key, false);
			variableModTypeForTreeViewObservableCollection->Add(theModType);
			for (auto uah : hm)
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				ModForTreeView tempVar2(uah->ToString(), false, uah->IdWithMotif, false, theModType);
				theModType->getChildren()->Add(&tempVar2);
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete theModType' statement was not added since theModType was passed to a method or constructor. Handle memory management manually.
		}
		variableModsTreeView->DataContext = variableModTypeForTreeViewObservableCollection;

		for (auto hm : GlobalVariables::getAllModsKnown().GroupBy([&] (std::any b)
		{
			b::ModificationType;
		}))
		{
			auto theModType = new ModTypeForTreeView(hm::Key, false);
			gptmdModTypeForTreeViewObservableCollection->Add(theModType);
			for (auto uah : hm)
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				ModForTreeView tempVar3(uah->ToString(), false, uah->IdWithMotif, false, theModType);
				theModType->getChildren()->Add(&tempVar3);
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete theModType' statement was not added since theModType was passed to a method or constructor. Handle memory management manually.
		}
		gptmdModsTreeView->DataContext = gptmdModTypeForTreeViewObservableCollection;
	}

	void GptmdTaskWindow::CancelButton_Click(std::any sender, RoutedEventArgs *e)
	{
		DialogResult = false;
	}

	void GptmdTaskWindow::SaveButton_Click(std::any sender, RoutedEventArgs *e)
	{
		std::wstring fieldNotUsed = L"1";

		if (!GlobalGuiSettings::CheckTaskSettingsValidity(precursorMassToleranceTextBox->Text, productMassToleranceTextBox->Text, missedCleavagesTextBox->Text, maxModificationIsoformsTextBox->Text, MinPeptideLengthTextBox->Text, MaxPeptideLengthTextBox->Text, maxThreadsTextBox->Text, minScoreAllowed->Text, fieldNotUsed, fieldNotUsed, DeconvolutionMaxAssumedChargeStateTextBox->Text, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed))

		{
			return;
		}

		Protease *protease = static_cast<Protease*>(proteaseComboBox::SelectedItem);
		int MaxMissedCleavages = missedCleavagesTextBox->Text->empty() ? std::numeric_limits<int>::max() : (std::stoi(missedCleavagesTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture));
		int MinPeptideLength = std::stoi(MinPeptideLengthTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture);
		int MaxPeptideLength = MaxPeptideLengthTextBox->Text->empty() ? std::numeric_limits<int>::max() : (std::stoi(MaxPeptideLengthTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture));
		int MinVariantDepth = std::stoi(MinVariantDepthTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture);
		int MaxHeterozygousVariants = std::stoi(MaxHeterozygousVariantsTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture);
		int MaxModificationIsoforms = std::stoi(maxModificationIsoformsTextBox->Text, CultureInfo::InvariantCulture);
		InitiatorMethionineBehavior *InitiatorMethionineBehavior = static_cast<InitiatorMethionineBehavior*>(initiatorMethionineBehaviorComboBox::SelectedIndex);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		DissociationType *dissociationType = GlobalVariables::getAllSupportedDissociationTypes()[DissociationTypeComboBox::SelectedItem.ToString()];

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

		auto listOfModsVariable = std::vector<(std::wstring, std::wstring)*>();
		for (auto heh : variableModTypeForTreeViewObservableCollection)
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
			delete PrecursorMassTolerance;
			delete ProductMassTolerance;
			return;
		}

		auto listOfModsFixed = std::vector<(std::wstring, std::wstring)*>();
		for (auto heh : fixedModTypeForTreeViewObservableCollection)
		{
			listOfModsFixed.AddRange(heh->getChildren()->Where([&] (std::any b)
			{
				b::Use;
			})->Select([&] (std::any b)
			{
				(b::Parent->DisplayName, b::ModName);
			}));
		}
		bool parseMaxThreadsPerFile = std::stoi(maxThreadsTextBox->Text, CultureInfo::InvariantCulture) <= Environment::ProcessorCount && std::stoi(maxThreadsTextBox->Text, CultureInfo::InvariantCulture) > 0;

		CommonParameters tempVar();
		DigestionParams tempVar2(protease: protease->Name, maxMissedCleavages: MaxMissedCleavages, minPeptideLength: MinPeptideLength, maxPeptideLength: MaxPeptideLength, maxModificationIsoforms: MaxModificationIsoforms, initiatorMethionineBehavior: InitiatorMethionineBehavior);
		CommonParameters *commonParamsToSave = new CommonParameters(OutputFileNameTextBox->Text != L"" ? OutputFileNameTextBox->Text: L"GPTMDTask", dissociationType, deconvolutePrecursors::IsChecked->Value, useProvidedPrecursor::IsChecked->Value, , std::stoi(DeconvolutionMaxAssumedChargeStateTextBox->Text, CultureInfo::InvariantCulture), , addCompIonCheckBox::IsChecked->Value, , std::stod(minScoreAllowed->Text, CultureInfo::InvariantCulture), , , , , , , ProductMassTolerance, PrecursorMassTolerance, , parseMaxThreadsPerFile ? std::stoi(maxThreadsTextBox->Text, CultureInfo::InvariantCulture): (&tempVar)->getMaxThreadsToUsePerFile(), &tempVar2, listOfModsVariable, listOfModsFixed, , protease->Name != L"top-down", MaxHeterozygousVariants, MinVariantDepth);

		getTheTask()->getGptmdParameters()->ListOfModsGptmd = std::vector<(std::wstring, std::wstring)*>();
		for (auto heh : gptmdModTypeForTreeViewObservableCollection)
		{
			getTheTask()->getGptmdParameters()->ListOfModsGptmd->AddRange(heh->getChildren()->Where([&] (std::any b)
			{
				b::Use;
			})->Select([&] (std::any b)
			{
				(b::Parent->DisplayName, b::ModName);
			}));
		}

		getTheTask()->setCommonParameters(commonParamsToSave);

		DialogResult = true;

//C# TO C++ CONVERTER TODO TASK: A 'delete commonParamsToSave' statement was not added since commonParamsToSave was assigned to another object. Handle memory management manually.
		delete PrecursorMassTolerance;
		delete ProductMassTolerance;
	}

	void GptmdTaskWindow::CheckIfNumber(std::any sender, TextCompositionEventArgs *e)
	{
		e->Handled = !GlobalGuiSettings::CheckIsNumber(e->Text);
	}

	void GptmdTaskWindow::KeyPressed(std::any sender, KeyEventArgs *e)
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

	void GptmdTaskWindow::TextChanged_Fixed(std::any sender, TextChangedEventArgs *args)
	{
		SearchModifications::SetTimer();
		SearchModifications::FixedSearch = true;
	}

	void GptmdTaskWindow::TextChanged_Var(std::any sender, TextChangedEventArgs *args)
	{
		SearchModifications::SetTimer();
		SearchModifications::VariableSearch = true;
	}

	void GptmdTaskWindow::TextChanged_GPTMD(std::any sender, TextChangedEventArgs *args)
	{
		SearchModifications::SetTimer();
		SearchModifications::GptmdSearch = true;
	}

	void GptmdTaskWindow::TextChangeTimerHandler(std::any sender, EventArgs *e)
	{
		if (SearchModifications::FixedSearch)
		{
			SearchModifications::FilterTree(SearchFixMod, fixedModsTreeView, fixedModTypeForTreeViewObservableCollection);
			SearchModifications::FixedSearch = false;
		}

		if (SearchModifications::VariableSearch)
		{
			SearchModifications::FilterTree(SearchVarMod, variableModsTreeView, variableModTypeForTreeViewObservableCollection);
			SearchModifications::VariableSearch = false;
		}

		if (SearchModifications::GptmdSearch)
		{
			SearchModifications::FilterTree(SearchGPTMD, gptmdModsTreeView, gptmdModTypeForTreeViewObservableCollection);
			SearchModifications::GptmdSearch = false;
		}
	}
}
