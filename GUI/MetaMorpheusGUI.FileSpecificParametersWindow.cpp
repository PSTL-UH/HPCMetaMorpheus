#include "MetaMorpheusGUI.FileSpecificParametersWindow.h"
#include "ForDisplayingInDataGrids/RawDataForDataGrid.h"
#include "../TaskLayer/FileSpecificParameters.h"
#include "Util/GlobalGuiSettings.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../EngineLayer/CommonParameters.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Nett;
using namespace TaskLayer;
using namespace Proteomics::ProteolyticDigestion;

namespace MetaMorpheusGUI
{

	FileSpecificParametersWindow::FileSpecificParametersWindow(ObservableCollection<RawDataForDataGrid*> *selectedSpectraFiles)
	{
		setSelectedSpectra(selectedSpectraFiles);
		InitializeComponent();
		PopulateChoices();
	}

	ObservableCollection<RawDataForDataGrid*> *FileSpecificParametersWindow::getSelectedSpectra() const
	{
		return privateSelectedSpectra;
	}

	void FileSpecificParametersWindow::setSelectedSpectra(ObservableCollection<RawDataForDataGrid*> *value)
	{
		privateSelectedSpectra = value;
	}

	void FileSpecificParametersWindow::Save_Click(std::any sender, RoutedEventArgs *e)
	{
		auto parametersToWrite = new FileSpecificParameters();

		// parse the file-specific parameters to text
		int paramsToSaveCount = 0;
		if (fileSpecificPrecursorMassTolEnabled::IsChecked->Value)
		{
			paramsToSaveCount++;
			if (GlobalGuiSettings::CheckPrecursorMassTolerance(precursorMassToleranceTextBox->Text))
			{
				double value = std::stod(precursorMassToleranceTextBox->Text);
				if (precursorMassToleranceComboBox->SelectedIndex == 0)
				{
					AbsoluteTolerance tempVar(value);
					parametersToWrite->setPrecursorMassTolerance(&tempVar);
				}
				else
				{
					PpmTolerance tempVar2(value);
					parametersToWrite->setPrecursorMassTolerance(&tempVar2);
				}
			}
			else
			{
				delete parametersToWrite;
				return;
			}
		}
		if (fileSpecificProductMassTolEnabled::IsChecked->Value)
		{
			paramsToSaveCount++;
			if (GlobalGuiSettings::CheckProductMassTolerance(productMassToleranceTextBox->Text))
			{
				double value = std::stod(productMassToleranceTextBox->Text);
				if (productMassToleranceComboBox->SelectedIndex == 0)
				{
					AbsoluteTolerance tempVar3(value);
					parametersToWrite->setProductMassTolerance(&tempVar3);
				}
				else
				{
					PpmTolerance tempVar4(value);
					parametersToWrite->setProductMassTolerance(&tempVar4);
				}
			}
			else
			{
				delete parametersToWrite;
				return;
			}
		}
		if (fileSpecificProteaseEnabled::IsChecked->Value)
		{
			paramsToSaveCount++;
			parametersToWrite->setProtease(static_cast<Protease*>(fileSpecificProtease::SelectedItem));
		}
		if (fileSpecificMinPeptideLengthEnabled::IsChecked->Value)
		{
			paramsToSaveCount++;
			int i;
			if (int::TryParse(MinPeptideLengthTextBox->Text, i) && i > 0)
			{
				parametersToWrite->setMinPeptideLength(i);
			}
			else
			{
				MessageBox::Show(L"The minimum peptide length must be a positive integer");

				delete parametersToWrite;
				return;
			}
		}
		if (fileSpecificMaxPeptideLengthEnabled::IsChecked->Value)
		{
			paramsToSaveCount++;
			std::wstring lengthMaxPeptide = GlobalGuiSettings::MaxValueConversion(MaxPeptideLengthTextBox->Text);
			if (GlobalGuiSettings::CheckPeptideLength(MinPeptideLengthTextBox->Text, lengthMaxPeptide))
			{
				parametersToWrite->setMaxPeptideLength(std::make_optional(std::stoi(lengthMaxPeptide)));
			}
			else
			{
				delete parametersToWrite;
				return;
			}
		}
		if (fileSpecificMissedCleavagesEnabled::IsChecked->Value)
		{
			paramsToSaveCount++;
			std::wstring lengthCleavage = GlobalGuiSettings::MaxValueConversion(missedCleavagesTextBox->Text);
			if (GlobalGuiSettings::CheckMaxMissedCleavages(lengthCleavage))
			{
				parametersToWrite->setMaxMissedCleavages(std::make_optional(std::stoi(lengthCleavage)));
			}
			else
			{
				delete parametersToWrite;
				return;
			}
		}
		if (fileSpecificMaxModNumEnabled::IsChecked->Value)
		{
			paramsToSaveCount++;
			if (GlobalGuiSettings::CheckMaxModsPerPeptide(MaxModNumTextBox->Text))
			{
				parametersToWrite->setMaxModsForPeptide(std::make_optional(std::stoi(MaxModNumTextBox->Text)));
			}
			else
			{
				delete parametersToWrite;
				return;
			}
		}
		//if (fileSpecificIonTypesEnabled.IsChecked.Value)
		//{
		//    paramsToSaveCount++;

		//    // don't think there's any way to mess up checkboxes... no error message needed
		//    parametersToWrite.BIons = bCheckBox.IsChecked;
		//    parametersToWrite.YIons = yCheckBox.IsChecked;
		//    parametersToWrite.CIons = cCheckBox.IsChecked;
		//    parametersToWrite.ZdotIons = zdotCheckBox.IsChecked;
		//}

		// write parameters to toml files for the selected spectra files


		auto tomlPathsForSelectedFiles = getSelectedSpectra()->Select([&] (std::any p)
		{
		delete parametersToWrite;
			return FileSystem::combine(Directory::GetParent(p::FilePath)->ToString(), Path::GetFileNameWithoutExtension(p::FileName)) + L".toml";
		});
		for (auto tomlToWrite : tomlPathsForSelectedFiles)
		{
			if (paramsToSaveCount > 0)
			{
				Toml::WriteFile(parametersToWrite, tomlToWrite, MetaMorpheusTask::tomlConfig);

				// make sure the settings are able to be parsed...
				auto tempTomlTable = Toml::ReadFile(tomlToWrite, MetaMorpheusTask::tomlConfig);
				FileSpecificParameters *tempParams = new FileSpecificParameters(tempTomlTable);

				delete tempParams;
			}
			else
			{
				// user has specified that no file-specific settings should be used; delete the file-specific toml if it exists
				File::Delete(tomlToWrite);
			}
		}

		// done
		DialogResult = true;

//C# TO C++ CONVERTER TODO TASK: A 'delete parametersToWrite' statement was not added since parametersToWrite was passed to a method or constructor. Handle memory management manually.
	}

	void FileSpecificParametersWindow::Cancel_Click(std::any sender, RoutedEventArgs *e)
	{
		DialogResult = false;
	}

	void FileSpecificParametersWindow::PopulateChoices()
	{
		// use default settings to populate
		auto defaultParams = new CommonParameters();
		Protease *tempProtease = defaultParams->getDigestionParams()->Protease;
		int tempMinPeptideLength = defaultParams->getDigestionParams()->MinPeptideLength;
		int tempMaxPeptideLength = defaultParams->getDigestionParams()->MaxPeptideLength;
		int tempMaxMissedCleavages = defaultParams->getDigestionParams()->MaxMissedCleavages;
		int tempMaxModsForPeptide = defaultParams->getDigestionParams()->MaxModsForPeptide;
		auto tempPrecursorMassTolerance = defaultParams->getPrecursorMassTolerance();
		auto tempProductMassTolerance = defaultParams->getProductMassTolerance();

		// do any of the selected files already have file-specific parameters specified?
		auto spectraFiles = getSelectedSpectra()->Select([&] (std::any p)
		{
			p::FilePath;
		});
		for (auto file : spectraFiles)
		{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			std::wstring tomlPath = FileSystem::combine(Directory::GetParent(file)->ToString(), Path::GetFileNameWithoutExtension(file)) + L".toml";

			if (FileSystem::fileExists(tomlPath))
			{
				TomlTable *tomlTable = Toml::ReadFile(tomlPath, MetaMorpheusTask::tomlConfig);
				FileSpecificParameters *fileSpecificParams = new FileSpecificParameters(tomlTable);

				if (fileSpecificParams->getPrecursorMassTolerance() != nullptr)
				{
					tempPrecursorMassTolerance = fileSpecificParams->getPrecursorMassTolerance();
					fileSpecificPrecursorMassTolEnabled->IsChecked = true;
				}
				if (fileSpecificParams->getProductMassTolerance() != nullptr)
				{
					tempProductMassTolerance = fileSpecificParams->getProductMassTolerance();
					fileSpecificProductMassTolEnabled->IsChecked = true;
				}
				if (fileSpecificParams->getProtease() != nullptr)
				{
					tempProtease = (fileSpecificParams->getProtease());
					fileSpecificProteaseEnabled->IsChecked = true;
				}
				if (fileSpecificParams->getMinPeptideLength() != nullptr)
				{
					tempMinPeptideLength = (fileSpecificParams->getMinPeptideLength().Value);
					fileSpecificMinPeptideLengthEnabled->IsChecked = true;
				}
				if (fileSpecificParams->getMaxPeptideLength() != nullptr)
				{
					tempMaxPeptideLength = (fileSpecificParams->getMaxPeptideLength().Value);
					fileSpecificMaxPeptideLengthEnabled->IsChecked = true;
				}
				if (fileSpecificParams->getMaxMissedCleavages() != nullptr)
				{
					tempMaxMissedCleavages = (fileSpecificParams->getMaxMissedCleavages().Value);
					fileSpecificMissedCleavagesEnabled->IsChecked = true;
				}
				if (fileSpecificParams->getMaxModsForPeptide() != nullptr)
				{
					tempMaxModsForPeptide = (fileSpecificParams->getMaxMissedCleavages().Value);
					fileSpecificMaxModNumEnabled->IsChecked = true;
				}

//C# TO C++ CONVERTER TODO TASK: A 'delete fileSpecificParams' statement was not added since fileSpecificParams was assigned to an outer scope variable. Handle memory management manually.
			}
		}

		DigestionParams *digestParams = new DigestionParams(protease: tempProtease->Name, maxMissedCleavages: tempMaxMissedCleavages, minPeptideLength: tempMinPeptideLength, maxPeptideLength: tempMaxPeptideLength, maxModsForPeptides: tempMaxModsForPeptide);

		// populate the GUI
		for (Protease *protease : ProteaseDictionary::Dictionary->Values)
		{
			fileSpecificProtease::Items->Add(protease);
		}

		fileSpecificProtease->SelectedItem = digestParams->Protease;

		productMassToleranceComboBox::Items->Add(L"Da");
		productMassToleranceComboBox::Items->Add(L"ppm");
		productMassToleranceComboBox->SelectedIndex = 1;

		precursorMassToleranceComboBox::Items->Add(L"Da");
		precursorMassToleranceComboBox::Items->Add(L"ppm");
		precursorMassToleranceComboBox->SelectedIndex = 1;

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		precursorMassToleranceTextBox->Text = tempPrecursorMassTolerance->Value->ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		productMassToleranceTextBox->Text = tempProductMassTolerance->Value->ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MinPeptideLengthTextBox->Text = digestParams->MinPeptideLength.ToString();

		if (std::numeric_limits<int>::max() != digestParams->MaxPeptideLength)
		{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			MaxPeptideLengthTextBox->Text = digestParams->MaxPeptideLength.ToString();
		}

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MaxModNumTextBox->Text = digestParams->MaxModsForPeptide.ToString();
		if (std::numeric_limits<int>::max() != digestParams->MaxMissedCleavages)
		{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			missedCleavagesTextBox->Text = digestParams->MaxMissedCleavages.ToString();
		}
		//yCheckBox.IsChecked = tempCommonParams.YIons;
		//bCheckBox.IsChecked = tempCommonParams.BIons;
		//cCheckBox.IsChecked = tempCommonParams.CIons;
		//zdotCheckBox.IsChecked = tempCommonParams.ZdotIons;

		delete digestParams;
		delete defaultParams;
	}

	void FileSpecificParametersWindow::CheckIfNumber(std::any sender, TextCompositionEventArgs *e)
	{
		e->Handled = !GlobalGuiSettings::CheckIsNumber(e->Text);
	}

	void FileSpecificParametersWindow::KeyPressed(std::any sender, KeyEventArgs *e)
	{
		if (e->Key == Key->Return)
		{
			Save_Click(sender, e);
		}
		else if (e->Key == Key->Escape)
		{
			Cancel_Click(sender, e);
		}
	}
}
