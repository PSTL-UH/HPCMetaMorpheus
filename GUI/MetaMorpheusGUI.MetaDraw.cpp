#include "MetaMorpheusGUI.MetaDraw.h"
#include "MetaDraw/PsmAnnotationViewModel.h"
#include "../TaskLayer/MyFileManager.h"
#include "../EngineLayer/MetaDraw/MetaDrawPsm.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/MetaDraw/TsvResultReader.h"
#include "MetaDraw/BaseSequenceAnnotation.h"

using namespace ViewModels;
using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace TaskLayer;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace MetaMorpheusGUI
{

	MetaDraw::MetaDraw() : peptideSpectralMatches(new ObservableCollection<MetaDrawPsm*>()), propertyView(new DataTable())
	{
		InitializeComponent();

		mainViewModel = new PsmAnnotationViewModel();
		plotView->DataContext = mainViewModel;
		propertyView->Columns->Add(L"Name", std::wstring::typeid);
		propertyView->Columns->Add(L"Value", std::wstring::typeid);
		peptideSpectralMatchesView = CollectionViewSource::GetDefaultView(peptideSpectralMatches);
		dataGridScanNums->DataContext = peptideSpectralMatchesView;
		dataGridProperties->DataContext = propertyView->DefaultView;
		Title = L"MetaDraw: version " + GlobalVariables::getMetaMorpheusVersion();
		spectraFileManager = new MyFileManager(true);
		SetUpDictionaries();
		modificationAnnotationColor = Brushes::Yellow;
	}

	void MetaDraw::SetUpDictionaries()
	{
		// colors of each fragment to annotate on base sequence
		productTypeToColor = (static_cast<std::vector<ProductType*>>(Enum::GetValues(ProductType::typeid)))->ToDictionary([&] (std::any p)
		{
			return p;
		}, [&] (std::any p)
		{
			Colors::Aqua;
		});
		productTypeToColor[ProductType::b] = Colors::Blue;
		productTypeToColor[ProductType::y] = Colors::Purple;
		productTypeToColor[ProductType::zDot] = Colors::Orange;
		productTypeToColor[ProductType::c] = Colors::Gold;

		// offset for annotation on base sequence
		productTypeToYOffset = (static_cast<std::vector<ProductType*>>(Enum::GetValues(ProductType::typeid)))->ToDictionary([&] (std::any p)
		{
			return p;
		}, [&] (std::any p)
		{
			0.0;
		});
		productTypeToYOffset[ProductType::b] = 50;
		productTypeToYOffset[ProductType::y] = 0;
		productTypeToYOffset[ProductType::c] = 50;
		productTypeToYOffset[ProductType::zDot] = 0;
	}

	void MetaDraw::Window_Drop(std::any sender, DragEventArgs *e)
	{
		std::vector<std::wstring> files = (static_cast<std::vector<std::wstring>>(e->Data->GetData(DataFormats::FileDrop)))->OrderBy([&] (std::any p)
		{
			return p;
		})->ToArray();

		if (files.size() > 0)
		{
			for (auto draggedFilePath : files)
			{
				if (FileSystem::fileExists(draggedFilePath))
				{
					LoadFile(draggedFilePath);
				}
			}
		}
	}

	void MetaDraw::LoadFile(const std::wstring &filePath)
	{
		auto theExtension = StringHelper::toLower(Path::GetExtension(filePath));

//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//		switch (theExtension)
//ORIGINAL LINE: case ".raw":
		if (theExtension == L".raw" || theExtension == L".mzml" || theExtension == L".mgf")
		{
				spectraFilePath = filePath;
				spectraFileNameLabel->Text = filePath;
		}
//ORIGINAL LINE: case ".psmtsv":
		else if (theExtension == L".psmtsv" || theExtension == L".tsv")
		{
				tsvResultsFilePath = filePath;
				psmFileNameLabel->Text = filePath;
		}
		else
		{
				MessageBox::Show(L"Cannot read file type: " + theExtension);
		}
	}

	void MetaDraw::LoadPsms(const std::wstring &filename)
	{
		std::wstring fileNameWithExtension = FileSystem::getFileName(spectraFilePath);
		std::wstring fileNameWithoutExtension = Path::GetFileNameWithoutExtension(spectraFilePath);

		try
		{
			std::vector<std::wstring> warnings; // TODO: print warnings
			for (auto psm : TsvResultReader::ReadTsv(filename, warnings))
			{
				if (psm->getFilename() == fileNameWithExtension || psm->getFilename() == fileNameWithoutExtension || psm->getFilename().find(fileNameWithoutExtension) != std::wstring::npos)
				{
					std::function<void()> tempVar([&] ()
					{
						peptideSpectralMatches->Add(psm);
					});
					Dispatcher::BeginInvoke(&tempVar);
				}
			}
		}
		catch (const std::runtime_error &e)
		{
			MessageBox::Show(L"Could not open PSM file:\n" + e.what());
		}
	}

	void MetaDraw::DrawPsm(int oneBasedScanNumber, const std::wstring &fullSequence)
	{
		MsDataScan *msDataScanToDraw = MsDataFile->GetOneBasedScan(oneBasedScanNumber);
		std::vector<MetaDrawPsm*> scanPsms = peptideSpectralMatches->Where([&] (std::any p)
		{
			return p->Ms2ScanNumber == oneBasedScanNumber;
		});

		if (fullSequence != L"")
		{
			scanPsms = scanPsms.Where([&] (std::any p)
			{
				return p->FullSequence == fullSequence;
			});
		}

		MetaDrawPsm *psmToDraw = scanPsms.FirstOrDefault();

		// draw annotated spectrum
		mainViewModel->DrawPeptideSpectralMatch(msDataScanToDraw, psmToDraw);

		// draw annotated base sequence
		//TO DO: Annotate crosslinked peptide sequence           
		if (psmToDraw->getCrossType() == L"") // if the psm is single peptide (not crosslinked).
		{
			DrawAnnotatedBaseSequence(psmToDraw);
		}
	}

	void MetaDraw::dataGridScanNums_SelectedCellsChanged(std::any sender, SelectedCellsChangedEventArgs *e)
	{
		if (dataGridScanNums->SelectedItem == nullptr)
		{
			return;
		}

		// draw the selected PSM
		propertyView->Clear();
		MetaDrawPsm *row = static_cast<MetaDrawPsm*>(dataGridScanNums::SelectedItem);
		std::vector<System::Reflection::PropertyInfo*> temp = row->GetType()->GetProperties();

		for (int i = 0; i < temp.size(); i++)
		{
			if (temp[i]->Name == L"MatchedIons")
			{
				propertyView->Rows->Add(temp[i]->Name, std::wstring::Join(L", ", row->getMatchedIons().Select([&] (std::any p)
				{
					p::Annotation;
				})));
			}
			else
			{
				propertyView->Rows->Add(temp[i]->Name, temp[i]->GetValue(row, nullptr));
			}
		}
		dataGridProperties::Items->Refresh();
		DrawPsm(row->getMs2ScanNumber(), row->getFullSequence());
	}

	void MetaDraw::selectSpectraFileButton_Click(std::any sender, RoutedEventArgs *e)
	{
		Microsoft::Win32::OpenFileDialog *openFileDialog1 = new Microsoft::Win32::OpenFileDialog();
		openFileDialog1->Filter = L"Spectra Files(*.raw;*.mzML)|*.raw;*.mzML";
		openFileDialog1->FilterIndex = 1;
		openFileDialog1->RestoreDirectory = true;
		openFileDialog1->Multiselect = false;
		if (openFileDialog1->ShowDialog() == true)
		{
			for (auto filePath : openFileDialog1->FileNames.OrderBy([&] (std::any p)
			{
			delete openFileDialog1;
				return p;
			}))
			{
				LoadFile(filePath);
			}
		}

		delete openFileDialog1;
	}

	void MetaDraw::selectPsmFileButton_Click(std::any sender, RoutedEventArgs *e)
	{
		Microsoft::Win32::OpenFileDialog *openFileDialog1 = new Microsoft::Win32::OpenFileDialog();
		openFileDialog1->Filter = L"Result Files(*.psmtsv)|*.psmtsv";
		openFileDialog1->FilterIndex = 1;
		openFileDialog1->RestoreDirectory = true;
		openFileDialog1->Multiselect = false;
		if (openFileDialog1->ShowDialog() == true)
		{
			for (auto filePath : openFileDialog1->FileNames.OrderBy([&] (std::any p)
			{
			delete openFileDialog1;
				return p;
			}))
			{
				LoadFile(filePath);
			}
		}

		delete openFileDialog1;
	}

//C# TO C++ CONVERTER TODO TASK: There is no equivalent in C++ to the 'async' keyword:
//ORIGINAL LINE: private async void loadFilesButton_Click(object sender, RoutedEventArgs e)
	void MetaDraw::loadFilesButton_Click(std::any sender, RoutedEventArgs *e)
	{
		// check for validity
		propertyView->Clear();
		if (spectraFilePath == L"")
		{
			MessageBox::Show(L"Please add a spectra file.");
			return;
		}

		if (tsvResultsFilePath == L"")
		{
			MessageBox::Show(L"Please add a search result file.");
			return;
		}

		// load the spectra file
		(dynamic_cast<Button*>(sender))->IsEnabled = false;
		selectSpectraFileButton->IsEnabled = false;
		selectPsmFileButton->IsEnabled = false;
		prgsFeed->IsOpen = true;
		prgsText->Content = L"Loading spectra file...";

		auto slowProcess = Task<MsDataFile*>::Factory->StartNew([&] ()
		{
			CommonParameters tempVar();
			spectraFileManager->LoadFile(spectraFilePath, std::nullopt, std::nullopt, false, false, &tempVar);
		});
//C# TO C++ CONVERTER TODO TASK: There is no equivalent to 'await' in C++:
		await slowProcess;
		MsDataFile = slowProcess->Result;

		// load the PSMs
		this->prgsText->Content = L"Loading PSMs...";
//C# TO C++ CONVERTER TODO TASK: There is no equivalent to 'await' in C++:
		await Task::Run([&] ()
		{
			LoadPsms(tsvResultsFilePath);
		});

		// done loading - restore controls
		this->prgsFeed->IsOpen = false;
		(dynamic_cast<Button*>(sender))->IsEnabled = true;
		selectSpectraFileButton->IsEnabled = true;
		selectPsmFileButton->IsEnabled = true;
	}

	void MetaDraw::TextBox_TextChanged(std::any sender, TextChangedEventArgs *e)
	{
		std::wstring txt = (dynamic_cast<TextBox*>(sender))->Text;
		if (txt == L"")
		{
			peptideSpectralMatchesView->Filter = nullptr;
		}
		else
		{
			peptideSpectralMatchesView->Filter = [&] (std::any obj)
			{
				MetaDrawPsm *psm = dynamic_cast<MetaDrawPsm*>(obj);
				return ((std::to_wstring(psm->getMs2ScanNumber()))->StartsWith(txt) || StringHelper::toUpper(psm->getFullSequence()).find(txt.ToUpper()) != std::wstring::npos);
			};
		}
	}

	void MetaDraw::DrawAnnotatedBaseSequence(MetaDrawPsm *psm)
	{
		double spacing = 22;
		BaseDraw::clearCanvas(canvas);

		// don't draw ambiguous sequences
		if (psm->getFullSequence().find(L"|") != std::wstring::npos)
		{
			return;
		}

		// draw base sequence
		for (int r = 0; r < psm->getBaseSeq().length(); r++)
		{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			Point tempVar(r * spacing + 10, 10);
			BaseDraw::txtDrawing(canvas, &tempVar, psm->getBaseSeq()[r].ToString(), Brushes::Black);
		}

		// draw the fragment ion annotations on the base sequence
		for (auto ion : psm->getMatchedIons())
		{
			int residue = ion->NeutralTheoreticalProduct.TerminusFragment.AminoAcidPosition;
			std::wstring annotation = ion->NeutralTheoreticalProduct.ProductType + L"" + ion->NeutralTheoreticalProduct.TerminusFragment.FragmentNumber;

			if (ion->NeutralTheoreticalProduct.NeutralLoss != 0)
			{
				annotation += L"-" + ion->NeutralTheoreticalProduct.NeutralLoss;
			}

			if (ion->NeutralTheoreticalProduct.TerminusFragment->Terminus == FragmentationTerminus::C)
			{
				Point tempVar2(residue * spacing + 8, productTypeToYOffset[ion->NeutralTheoreticalProduct.ProductType]);
				BaseDraw::topSplittingDrawing(canvas, &tempVar2, productTypeToColor[ion->NeutralTheoreticalProduct.ProductType], annotation);
			}
			else if (ion->NeutralTheoreticalProduct.TerminusFragment->Terminus == FragmentationTerminus::N)
			{
				Point tempVar3(residue * spacing + 8, productTypeToYOffset[ion->NeutralTheoreticalProduct.ProductType]);
				BaseDraw::botSplittingDrawing(canvas, &tempVar3, productTypeToColor[ion->NeutralTheoreticalProduct.ProductType], annotation);
			}
			// don't draw diagnostic ions, precursor ions, etc
		}

		// draw modifications
		auto peptide = new PeptideWithSetModifications(psm->getFullSequence(), GlobalVariables::getAllModsKnownDictionary());
		for (auto mod : peptide->AllModsOneIsNterminus)
		{
			Point tempVar4((mod->Key - 1) * spacing - 17, 12);
			BaseDraw::circledTxtDraw(canvas, &tempVar4, modificationAnnotationColor);
		}

		delete peptide;
	}

	void MetaDraw::dataGridProperties_SelectedCellsChanged(std::any sender, SelectedCellsChangedEventArgs *e)
	{
		(dynamic_cast<DataGrid*>(sender))->UnselectAll();
	}
}
