#include "MetaMorpheusGUI.CustomModButtonWindow.h"
#include "../EngineLayer/GlobalVariables.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace Proteomics;

namespace MetaMorpheusGUI
{

std::unordered_map<std::wstring, std::wstring> CustomModButtonWindow::locationRestrictions;

	CustomModButtonWindow::CustomModButtonWindow()
	{
		InitializeComponent();

		if (locationRestrictions.empty())
		{
			locationRestrictions = std::unordered_map<std::wstring, std::wstring>();
			locationRestrictions.emplace(L"Anywhere", L"Anywhere.");
			locationRestrictions.emplace(L"Peptide N-Terminus", L"Peptide N-terminal.");
			locationRestrictions.emplace(L"Peptide C-Terminus", L"Peptide C-terminal.");
			locationRestrictions.emplace(L"Protein N-Terminus", L"N-terminal.");
			locationRestrictions.emplace(L"Protein C-Terminus", L"C-terminal.");
		}

		for (auto locationRestriction : locationRestrictions)
		{
			locationRestrictionComboBox::Items->Add(locationRestriction.first);
		}

		for (auto type : GlobalVariables::getAllSupportedDissociationTypes())
		{
			dissociationTypeComboBox::Items->Add(type->second);
		}

		locationRestrictionComboBox->SelectedItem = L"Anywhere";
		dissociationTypeComboBox->SelectedItem = DissociationType::HCD;
	}

	void CustomModButtonWindow::SaveCustomMod_Click(std::any sender, RoutedEventArgs *e)
	{
		std::wstring modsDirectory = FileSystem::combine(GlobalVariables::getDataDir(), LR"(Mods)");
		std::wstring customModsPath = FileSystem::combine(modsDirectory, LR"(CustomModifications.txt)");
		std::vector<std::wstring> customModsText;

		if (!FileSystem::fileExists(customModsPath))
		{
			customModsText.push_back(L"Custom Modifications");
		}
		else
		{
			customModsText = File::ReadAllLines(customModsPath).ToList();
		}

		std::wstring idText = originalIdTextBox->Text;
		std::wstring motifText = motifTextBox->Text;
		std::wstring chemicalFormulaText = chemicalFormulaTextBox->Text;
		std::wstring modMassText = modMassTextBox->Text;
		std::wstring neutralLossText = neutralLossTextBox->Text;
		std::wstring diagnosticIonText = diagnosticIonTextBox->Text;
		std::wstring modificationTypeText = modificationTypeTextBox->Text;
		std::wstring locationRestriction = locationRestrictions[locationRestrictionComboBox->Text];
		DissociationType *disType = GlobalVariables::getAllSupportedDissociationTypes()[dissociationTypeComboBox->Text];

		if (ErrorsDetected(idText, motifText, modMassText, chemicalFormulaText, neutralLossText, modificationTypeText, diagnosticIonText))
		{
			return;
		}

		// create custom mod
		std::unordered_map<DissociationType*, std::vector<double>> neutralLosses;
		if (!neutralLossText.empty())
		{
			neutralLosses = std::unordered_map<DissociationType*, std::vector<double>>
			{
				{disType, StringHelper::split(neutralLossText, L',')->Select(double::Parse).ToList()}
			};
		}

		std::unordered_map<DissociationType*, std::vector<double>> diagnosticIons;
		if (!diagnosticIonText.empty())
		{
			diagnosticIons = std::unordered_map<DissociationType*, std::vector<double>>() {{disType, StringHelper::split(diagnosticIonText, L',')->Select([&] (std::any p)
			{
				std::stod(p).ToMass(1);
			}).ToList()}};
		}

		ModificationMotif finalMotif;
		ModificationMotif::TryGetMotif(motifText, finalMotif);

		ChemicalFormula *chemicalFormula = nullptr;
		if (!chemicalFormulaText.empty())
		{
			chemicalFormula = ChemicalFormula::ParseFormula(chemicalFormulaText);
		}

		std::optional<double> modMass = std::nullopt;
		if (!modMassText.empty())
		{
			modMass = std::make_optional(std::stod(modMassText));
		}

		Modification *modification = new Modification(_originalId: idText, _modificationType: modificationTypeText, _target: finalMotif, _locationRestriction: locationRestriction, _chemicalFormula: chemicalFormula, _monoisotopicMass: modMass, _neutralLosses: neutralLosses, _diagnosticIons: diagnosticIons);

		if (GlobalVariables::getAllModsKnownDictionary().find(modification->IdWithMotif) != GlobalVariables::getAllModsKnownDictionary().end())
		{
			MessageBox::Show(L"A modification already exists with the name: " + modification->IdWithMotif, L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);

			delete modification;
			return;
		}

		// write custom mod to mods file

		// write/read temp file to make sure the mod is readable, then delete it
		std::wstring tempPath = FileSystem::combine(modsDirectory, LR"(temp.txt)");
		try
		{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			std::vector<std::wstring> temp = {modification->ToString(), LR"(//)"};
			File::WriteAllLines(tempPath, temp);
			std::any errors;
			auto parsedMods = UsefulProteomicsDatabases::PtmListLoader::ReadModsFromFile(tempPath, errors);

			if (parsedMods->Count() != 1)
			{
				MessageBox::Show(L"Problem parsing custom mod: One mod was expected, a different number was generated", L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);

				delete modification;
				return;
			}

			if (errors::Any())
			{
				std::wstring concatErrors = std::wstring::Join(L"\r\n", errors->Select([&] (std::any p)
				{
					p::Item2;
				}));
				MessageBox::Show(L"Problem(s) parsing custom mod: " + L"\r\n" + concatErrors, L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);

				delete modification;
				return;
			}

			File::Delete(tempPath);
		}
		catch (const std::runtime_error &ex)
		{
			MessageBox::Show(L"Problem parsing custom mod: " + ex.what(), L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);
			File::Delete(tempPath);

			delete modification;
			return;
		}

		// delete old custom mods file, write new one
		try
		{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			customModsText.push_back(modification->ToString());
			customModsText.push_back(LR"(//)");
			File::Delete(customModsPath);
			File::WriteAllLines(customModsPath, customModsText);
		}
		catch (const std::runtime_error &ex)
		{
			MessageBox::Show(L"Problem saving custom mod to file: " + ex.what(), L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);

			delete modification;
			return;
		}

		GlobalVariables::AddMods({modification}, false);

		DialogResult = true;

//C# TO C++ CONVERTER TODO TASK: A 'delete modification' statement was not added since modification was passed to a method or constructor. Handle memory management manually.
	}

	void CustomModButtonWindow::CancelCustomMod_Click(std::any sender, RoutedEventArgs *e)
	{
		DialogResult = false;
	}

	bool CustomModButtonWindow::ErrorsDetected(const std::wstring &myModName, const std::wstring &motif, const std::wstring &mass, const std::wstring &chemFormula, const std::wstring &neutralLoss, const std::wstring &modType, const std::wstring &diagnosticIon)
	{
		// parse input
		if (myModName.empty() || modType.empty())
		{
			MessageBox::Show(L"The mod name and type need to be specified", L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);
			return true;
		}

		if (modType.find(L':') != std::wstring::npos)
		{
			MessageBox::Show(L"Modification Type cannot contain ':'", L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);
			return true;
		}

		if (!motif.empty())
		{
			if (motif.Count(wchar_t::IsUpper) != 1)
			{
				MessageBox::Show(L"Motif must contain exactly one uppercase letter", L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);
				return true;
			}
		}
		else
		{
			MessageBox::Show(L"Motif must be defined", L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);
			return true;
		}

		if (mass.empty() && chemFormula.empty())
		{
			MessageBox::Show(L"Either the mass or chemical formula needs to be specified", L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);
			return true;
		}

		if (!chemFormula.empty())
		{
			try
			{
				ChemicalFormula *cf = ChemicalFormula::ParseFormula(chemFormula);
			}
			catch (...)
			{
				MessageBox::Show(L"Could not parse chemical formula", L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);
				return true;
			}
		}

		double dmass;
		if (!mass.empty() && !double::TryParse(mass, dmass))
		{
			MessageBox::Show(L"Could not parse modification mass", L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);
			return true;
		}

		try
		{
			if (!neutralLoss.empty())
			{
				StringHelper::split(neutralLoss, L',')->Select(double::Parse).ToList();
			}
			if (!diagnosticIon.empty())
			{
				StringHelper::split(diagnosticIon, L',')->Select(double::Parse).ToList();
			}
		}
		catch (...)
		{
			MessageBox::Show(L"Neutral losses and diagnostic ions must be entered as numbers separated by ','", L"Error", MessageBoxButton::OK, MessageBoxImage::Hand);
			return true;
		}

		return false;
	}
}
