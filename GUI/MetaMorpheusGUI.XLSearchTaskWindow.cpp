#include "MetaMorpheusGUI.XLSearchTaskWindow.h"
#include "../TaskLayer/XLSearchTask/TaskLayer.XLSearchTask.h"
#include "SearchTaskWindow.xaml.h"
#include "ForDisplayingInDataGrids/SearchModeForDataGrid.h"
#include "Util/ModTypeForTreeView.h"
#include "Util/SearchModifications.h"
#include "Util/GlobalGuiSettings.h"
#include "../EngineLayer/GlobalVariables.h"
#include "Util/ModForTreeView.h"
#include "../EngineLayer/CommonParameters.h"

using namespace EngineLayer;
using namespace EngineLayer::CrosslinkSearch;
using namespace MzLibUtil;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;
using namespace Proteomics::ProteolyticDigestion;
using namespace MassSpectrometry;

namespace MetaMorpheusGUI
{

	XLSearchTaskWindow::XLSearchTaskWindow() : XLSearchTaskWindow(nullptr)
	{
	}

	XLSearchTaskWindow::XLSearchTaskWindow(XLSearchTask *task) : DataContextForSearchTaskWindow(new DataContextForSearchTaskWindow() { ExpanderTitle = std::wstring::Join(L", ", SearchModesForThisTask->Where([&] (std::any b)
	{
			b::Use;
	})->Select([&] (std::any b)
	{
			b->Name;
		})), AnalysisExpanderTitle = L"Some analysis properties...", SearchModeExpanderTitle = L"Some search properties..." })
		{
		InitializeComponent();
		PopulateChoices();
		XLSearchTask tempVar();
		setTheTask((task != nullptr) ? task : &tempVar);
		UpdateFieldsFromTask(getTheTask());

		if (task == nullptr)
		{
			this->saveButton->Content = L"Add the XLSearch Task";
		}
		this->DataContext = DataContextForSearchTaskWindow;
		SearchModifications::Timer->Tick += new EventHandler(TextChangeTimerHandler);
	}

	XLSearchTask *XLSearchTaskWindow::getTheTask() const
	{
		return privateTheTask;
	}

	void XLSearchTaskWindow::setTheTask(XLSearchTask *value)
	{
		privateTheTask = value;
	}

	void XLSearchTaskWindow::CheckIfNumber(std::any sender, TextCompositionEventArgs *e)
	{
		e->Handled = !GlobalGuiSettings::CheckIsNumber(e->Text);
	}

	void XLSearchTaskWindow::PopulateChoices()
	{
		for (auto crosslinkerName : Enum::GetNames(CrosslinkerType::typeid))
		{
			cbCrosslinker::Items->Add(crosslinkerName);
		}

		for (auto dissassociationType : GlobalVariables::getAllSupportedDissociationTypes())
		{
			DissociationTypeComboBox::Items->Add(dissassociationType.first);
		}

		cbbXLprecusorMsTl::Items->Add(L"Da");
		cbbXLprecusorMsTl::Items->Add(L"ppm");

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

		productMassToleranceComboBox::Items->Add(L"Da");
		productMassToleranceComboBox::Items->Add(L"ppm");

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

	void XLSearchTaskWindow::UpdateFieldsFromTask(XLSearchTask *task)
	{
		//Crosslink search para
		//RbSearchCrosslink.IsChecked = !task.XlSearchParameters.SearchGlyco;
		//RbSearchGlyco.IsChecked = task.XlSearchParameters.SearchGlyco;
		//CkbSearchGlycoWithBgYgIndex.IsChecked = task.XlSearchParameters.SearchGlycoWithBgYgIndex;
		cbCrosslinker->SelectedIndex = static_cast<int>(task->getXlSearchParameters()->getCrosslinkerType());
		ckbXLTopNum->IsChecked = task->getXlSearchParameters()->getRestrictToTopNHits();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		txtXLTopNum->Text = task->getXlSearchParameters()->getCrosslinkSearchTopNum().ToString(CultureInfo::InvariantCulture);
		ckbQuenchH2O->IsChecked = task->getXlSearchParameters()->getXlQuench_H2O();
		ckbQuenchNH2->IsChecked = task->getXlSearchParameters()->getXlQuench_NH2();
		ckbQuenchTris->IsChecked = task->getXlSearchParameters()->getXlQuench_Tris();
		txtUdXLKerName->Text = task->getXlSearchParameters()->getCrosslinkerName();
		ckbUdXLkerCleavable->IsChecked = task->getXlSearchParameters()->getIsCleavable();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		txtUdXLkerTotalMs->Text = task->getXlSearchParameters()->getCrosslinkerTotalMass().HasValue ? task->getXlSearchParameters()->getCrosslinkerTotalMass().Value.ToString(CultureInfo::InvariantCulture) : L"";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		txtUdXLkerShortMass->Text = task->getXlSearchParameters()->getCrosslinkerShortMass().HasValue ? task->getXlSearchParameters()->getCrosslinkerShortMass().Value.ToString(CultureInfo::InvariantCulture) : L"";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		txtUdXLkerLongMass->Text = task->getXlSearchParameters()->getCrosslinkerLongMass().HasValue ? task->getXlSearchParameters()->getCrosslinkerLongMass().Value.ToString(CultureInfo::InvariantCulture) : L"";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		txtH2OQuenchMass->Text = task->getXlSearchParameters()->getCrosslinkerDeadEndMassH2O().HasValue ? task->getXlSearchParameters()->getCrosslinkerDeadEndMassH2O().Value.ToString(CultureInfo::InvariantCulture) : L"";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		txtNH2QuenchMass->Text = task->getXlSearchParameters()->getCrosslinkerDeadEndMassNH2().HasValue ? task->getXlSearchParameters()->getCrosslinkerDeadEndMassNH2().Value.ToString(CultureInfo::InvariantCulture) : L"";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		txtTrisQuenchMass->Text = task->getXlSearchParameters()->getCrosslinkerDeadEndMassTris().HasValue ? task->getXlSearchParameters()->getCrosslinkerDeadEndMassTris().Value.ToString(CultureInfo::InvariantCulture) : L"";

		txtUdXLkerAminoAcids->Text = task->getXlSearchParameters()->getCrosslinkerResidues();
		txtUdXLkerAminoAcids2->Text = task->getXlSearchParameters()->getCrosslinkerResidues2();
		cbbXLprecusorMsTl->SelectedIndex = dynamic_cast<AbsoluteTolerance*>(task->getCommonParameters()->getPrecursorMassTolerance()) != nullptr ? 0 : 1;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		XLPrecusorMsTlTextBox->Text = task->getCommonParameters()->getPrecursorMassTolerance()->Value->ToString(CultureInfo::InvariantCulture);
		trimMs1->IsChecked = task->getCommonParameters()->getTrimMs1Peaks();
		trimMsMs->IsChecked = task->getCommonParameters()->getTrimMsMsPeaks();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		TopNPeaksTextBox->Text = task->getCommonParameters()->getTopNpeaks().ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MinRatioTextBox->Text = task->getCommonParameters()->getMinRatio().ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		DissociationTypeComboBox->SelectedItem = task->getCommonParameters()->getDissociationType()->ToString();

		//ckbCharge_2_3.IsChecked = task.XlSearchParameters.XlCharge_2_3;
		checkBoxDecoy->IsChecked = task->getXlSearchParameters()->getDecoyType() != DecoyType::None;
		deconvolutePrecursors->IsChecked = task->getCommonParameters()->getDoPrecursorDeconvolution();
		useProvidedPrecursor->IsChecked = task->getCommonParameters()->getUseProvidedPrecursorInfo();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		missedCleavagesTextBox->Text = task->getCommonParameters()->getDigestionParams()->MaxMissedCleavages.ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MinPeptideLengthTextBox->Text = task->getCommonParameters()->getDigestionParams()->MinPeptideLength.ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MaxPeptideLengthTextBox->Text = task->getCommonParameters()->getDigestionParams()->MaxPeptideLength == std::numeric_limits<int>::max() ? L"" : task->getCommonParameters()->getDigestionParams()->MaxPeptideLength.ToString(CultureInfo::InvariantCulture);
		proteaseComboBox->SelectedItem = task->getCommonParameters()->getDigestionParams()->Protease;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		maxModificationIsoformsTextBox->Text = task->getCommonParameters()->getDigestionParams()->MaxModificationIsoforms.ToString(CultureInfo::InvariantCulture);
		initiatorMethionineBehaviorComboBox->SelectedIndex = static_cast<int>(task->getCommonParameters()->getDigestionParams()->InitiatorMethionineBehavior);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		productMassToleranceTextBox->Text = task->getCommonParameters()->getProductMassTolerance()->Value->ToString(CultureInfo::InvariantCulture);
		productMassToleranceComboBox->SelectedIndex = dynamic_cast<AbsoluteTolerance*>(task->getCommonParameters()->getProductMassTolerance()) != nullptr ? 0 : 1;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		DissociationTypeComboBox->SelectedItem = task->getCommonParameters()->getDissociationType()->ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		minScoreAllowed->Text = task->getCommonParameters()->getScoreCutoff().ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		numberOfDatabaseSearchesTextBox->Text = task->getCommonParameters()->getTotalPartitions().ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		maxThreadsTextBox->Text = task->getCommonParameters()->getMaxThreadsToUsePerFile().ToString(CultureInfo::InvariantCulture);

		ckbPercolator->IsChecked = task->getXlSearchParameters()->getWriteOutputForPercolator();
		ckbPepXML->IsChecked = task->getXlSearchParameters()->getWritePepXml();

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

		for (auto ye : VariableModTypeForTreeViewObservableCollection)
		{
			ye->VerifyCheckState();
		}
		for (auto ye : FixedModTypeForTreeViewObservableCollection)
		{
			ye->VerifyCheckState();
		}
	}

	void XLSearchTaskWindow::CancelButton_Click(std::any sender, RoutedEventArgs *e)
	{
		DialogResult = false;
	}

	void XLSearchTaskWindow::SaveButton_Click(std::any sender, RoutedEventArgs *e)
	{
		std::wstring fieldNotUsed = L"1";

		if (!GlobalGuiSettings::CheckTaskSettingsValidity(XLPrecusorMsTlTextBox->Text, productMassToleranceTextBox->Text, missedCleavagesTextBox->Text, maxModificationIsoformsTextBox->Text, MinPeptideLengthTextBox->Text, MaxPeptideLengthTextBox->Text, maxThreadsTextBox->Text, minScoreAllowed->Text, fieldNotUsed, fieldNotUsed, fieldNotUsed, TopNPeaksTextBox->Text, MinRatioTextBox->Text, numberOfDatabaseSearchesTextBox->Text, fieldNotUsed, fieldNotUsed, fieldNotUsed))
		{
			return;
		}

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		DissociationType *dissociationType = GlobalVariables::getAllSupportedDissociationTypes()[DissociationTypeComboBox::SelectedItem.ToString()];
		//TheTask.XlSearchParameters.SearchGlyco = RbSearchGlyco.IsChecked.Value;
		//TheTask.XlSearchParameters.SearchGlycoWithBgYgIndex = CkbSearchGlycoWithBgYgIndex.IsChecked.Value;
		getTheTask()->getXlSearchParameters()->setRestrictToTopNHits(ckbXLTopNum::IsChecked->Value);
		getTheTask()->getXlSearchParameters()->setCrosslinkSearchTopNum(std::stoi(txtXLTopNum->Text, CultureInfo::InvariantCulture));
		getTheTask()->getXlSearchParameters()->setCrosslinkerType(static_cast<CrosslinkerType>(cbCrosslinker::SelectedIndex));

		//TheTask.XlSearchParameters.XlCharge_2_3 = ckbCharge_2_3.IsChecked.Value;
		getTheTask()->getXlSearchParameters()->setXlQuench_H2O(ckbQuenchH2O::IsChecked->Value);
		getTheTask()->getXlSearchParameters()->setXlQuench_NH2(ckbQuenchNH2::IsChecked->Value);
		getTheTask()->getXlSearchParameters()->setXlQuench_Tris(ckbQuenchTris::IsChecked->Value);

		if (getTheTask()->getXlSearchParameters()->getCrosslinkerType() == CrosslinkerType::UserDefined)
		{
			getTheTask()->getXlSearchParameters()->setCrosslinkerName(txtUdXLKerName->Text);
			getTheTask()->getXlSearchParameters()->setIsCleavable(ckbUdXLkerCleavable::IsChecked->Value);
			getTheTask()->getXlSearchParameters()->setCrosslinkerResidues(txtUdXLkerAminoAcids->Text);
			getTheTask()->getXlSearchParameters()->setCrosslinkerResidues2(txtUdXLkerAminoAcids2->Text);
			getTheTask()->getXlSearchParameters()->setCrosslinkerLongMass(txtUdXLkerLongMass->Text->empty() ? static_cast<std::optional<double>>(nullptr): std::stod(txtUdXLkerLongMass->Text, CultureInfo::InvariantCulture));

			getTheTask()->getXlSearchParameters()->setCrosslinkerShortMass(txtUdXLkerShortMass->Text->empty() ? static_cast<std::optional<double>>(nullptr): std::stod(txtUdXLkerShortMass->Text, CultureInfo::InvariantCulture));

			getTheTask()->getXlSearchParameters()->setCrosslinkerTotalMass(txtUdXLkerTotalMs->Text->empty() ? static_cast<std::optional<double>>(nullptr): std::stod(txtUdXLkerTotalMs->Text, CultureInfo::InvariantCulture));

			getTheTask()->getXlSearchParameters()->setCrosslinkerDeadEndMassH2O(txtH2OQuenchMass->Text->empty() ? static_cast<std::optional<double>>(nullptr): std::stod(txtH2OQuenchMass->Text, CultureInfo::InvariantCulture));

			getTheTask()->getXlSearchParameters()->setCrosslinkerDeadEndMassNH2(txtNH2QuenchMass->Text->empty() ? static_cast<std::optional<double>>(nullptr): std::stod(txtNH2QuenchMass->Text, CultureInfo::InvariantCulture));

			getTheTask()->getXlSearchParameters()->setCrosslinkerDeadEndMassTris(txtTrisQuenchMass->Text->empty() ? static_cast<std::optional<double>>(nullptr): std::stod(txtTrisQuenchMass->Text, CultureInfo::InvariantCulture));
		}

		getTheTask()->getXlSearchParameters()->setDecoyType(checkBoxDecoy::IsChecked->Value ? DecoyType::Reverse : DecoyType::None);

		Protease *protease = static_cast<Protease*>(proteaseComboBox::SelectedItem);
		int MaxMissedCleavages = missedCleavagesTextBox->Text->empty() ? std::numeric_limits<int>::max() : (std::stoi(missedCleavagesTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture));
		int MinPeptideLength = (std::stoi(MinPeptideLengthTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture));
		int MaxPeptideLength = MaxPeptideLengthTextBox->Text->empty() ? std::numeric_limits<int>::max() : (std::stoi(MaxPeptideLengthTextBox->Text, NumberStyles::Any, CultureInfo::InvariantCulture));
		int MaxModificationIsoforms = (std::stoi(maxModificationIsoformsTextBox->Text, CultureInfo::InvariantCulture));
		InitiatorMethionineBehavior *InitiatorMethionineBehavior = (static_cast<InitiatorMethionineBehavior*>(initiatorMethionineBehaviorComboBox::SelectedIndex));
		DigestionParams *digestionParamsToSave = new DigestionParams(protease: protease->Name, maxMissedCleavages: MaxMissedCleavages, minPeptideLength: MinPeptideLength, maxPeptideLength: MaxPeptideLength, maxModificationIsoforms: MaxModificationIsoforms, initiatorMethionineBehavior: InitiatorMethionineBehavior);

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
		if (cbbXLprecusorMsTl->SelectedIndex == 0)
		{
			PrecursorMassTolerance = new AbsoluteTolerance(std::stod(XLPrecusorMsTlTextBox->Text, CultureInfo::InvariantCulture));
		}
		else
		{
			PrecursorMassTolerance = new PpmTolerance(std::stod(XLPrecusorMsTlTextBox->Text, CultureInfo::InvariantCulture));
		}

		getTheTask()->getXlSearchParameters()->setWriteOutputForPercolator(ckbPercolator::IsChecked->Value);
		getTheTask()->getXlSearchParameters()->setWritePepXml(ckbPepXML::IsChecked->Value);
		//TheTask.UseProvidedPrecursorInfo = useProvidedPrecursor.IsChecked.Value;

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

		CommonParameters *commonParamsToSave = new CommonParameters(OutputFileNameTextBox->Text != L"" ? OutputFileNameTextBox->Text: L"XLSearchTask", dissociationType, deconvolutePrecursors::IsChecked->Value, useProvidedPrecursor::IsChecked->Value, , , , , std::stoi(numberOfDatabaseSearchesTextBox->Text, CultureInfo::InvariantCulture), std::stod(minScoreAllowed->Text, CultureInfo::InvariantCulture), std::stoi(TopNPeaksTextBox->Text), std::stod(MinRatioTextBox->Text), trimMs1::IsChecked->Value, trimMsMs::IsChecked->Value, , , ProductMassTolerance, PrecursorMassTolerance, , , digestionParamsToSave, listOfModsVariable, listOfModsFixed, , protease->Name != L"top-down");

		getTheTask()->setCommonParameters(commonParamsToSave);

		DialogResult = true;

//C# TO C++ CONVERTER TODO TASK: A 'delete commonParamsToSave' statement was not added since commonParamsToSave was assigned to another object. Handle memory management manually.
		delete PrecursorMassTolerance;
		delete ProductMassTolerance;
		delete digestionParamsToSave;
	}

	void XLSearchTaskWindow::ApmdExpander_Collapsed(std::any sender, RoutedEventArgs *e)
	{
		DataContextForSearchTaskWindow->setExpanderTitle(std::wstring::Join(L", ", SearchModesForThisTask->Where([&] (std::any b)
		{
			b::Use;
		})->Select([&] (std::any b)
		{
			b->Name;
		})));
		DataContextForSearchTaskWindow->setAnalysisExpanderTitle(L"Some analysis properties...");
		DataContextForSearchTaskWindow->setSearchModeExpanderTitle(L"Some search properties...");
	}

	void XLSearchTaskWindow::ModExpander_Expanded(std::any sender, RoutedEventArgs *e)
	{
		//Envent Handler not used #1
	}

	void XLSearchTaskWindow::ModificationsDataGrid_Loaded(std::any sender, RoutedEventArgs *e)
	{
		//Envent Handler not used #2
	}

	void XLSearchTaskWindow::ModificationsDataGrid_DataContextChanged(std::any sender, DependencyPropertyChangedEventArgs *e)
	{
		//Envent Handler not used #3
	}

	void XLSearchTaskWindow::ModificationsDataGrid_AutoGeneratedColumns(std::any sender, EventArgs *e)
	{
		//if (!TheTask.WritePrunedDatabase)
		//    modificationsDataGrid.Columns[3].Visibility = Visibility.Collapsed;
	}

	void XLSearchTaskWindow::KeyPressed(std::any sender, KeyEventArgs *e)
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

	void XLSearchTaskWindow::TextChanged_Fixed(std::any sender, TextChangedEventArgs *args)
	{
		SearchModifications::SetTimer();
		SearchModifications::FixedSearch = true;
	}

	void XLSearchTaskWindow::TextChanged_Var(std::any sender, TextChangedEventArgs *args)
	{
		SearchModifications::SetTimer();
		SearchModifications::VariableSearch = true;
	}

	void XLSearchTaskWindow::TextChangeTimerHandler(std::any sender, EventArgs *e)
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
