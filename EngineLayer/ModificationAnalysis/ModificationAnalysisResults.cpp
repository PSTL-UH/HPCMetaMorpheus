#include "ModificationAnalysisResults.h"
#include "ModificationAnalysisEngine.h"

using namespace Chemistry;
namespace EngineLayer
{
	namespace ModificationAnalysis
	{

		ModificationAnalysisResults::ModificationAnalysisResults(ModificationAnalysisEngine *modificationAnalysisEngine) : MetaMorpheusEngineResults(modificationAnalysisEngine)
		{
		}

		std::unordered_map<std::wstring, int> ModificationAnalysisResults::getCountOfAmbiguousButLocalizedModsSeen() const
		{
			return privateCountOfAmbiguousButLocalizedModsSeen;
		}

		void ModificationAnalysisResults::setCountOfAmbiguousButLocalizedModsSeen(const std::unordered_map<std::wstring, int> &value)
		{
			privateCountOfAmbiguousButLocalizedModsSeen = value;
		}

		std::unordered_map<std::wstring, int> ModificationAnalysisResults::getCountOfModsSeenAndLocalized() const
		{
			return privateCountOfModsSeenAndLocalized;
		}

		void ModificationAnalysisResults::setCountOfModsSeenAndLocalized(const std::unordered_map<std::wstring, int> &value)
		{
			privateCountOfModsSeenAndLocalized = value;
		}

		std::unordered_map<std::wstring, int> ModificationAnalysisResults::getCountOfEachModSeenOnProteins() const
		{
			return privateCountOfEachModSeenOnProteins;
		}

		void ModificationAnalysisResults::setCountOfEachModSeenOnProteins(const std::unordered_map<std::wstring, int> &value)
		{
			privateCountOfEachModSeenOnProteins = value;
		}

		std::unordered_map<std::wstring, int> ModificationAnalysisResults::getCountOfUnlocalizedMods() const
		{
			return privateCountOfUnlocalizedMods;
		}

		void ModificationAnalysisResults::setCountOfUnlocalizedMods(const std::unordered_map<std::wstring, int> &value)
		{
			privateCountOfUnlocalizedMods = value;
		}

		std::unordered_map<ChemicalFormula*, int> ModificationAnalysisResults::getCountOfUnlocalizedFormulas() const
		{
			return privateCountOfUnlocalizedFormulas;
		}

		void ModificationAnalysisResults::setCountOfUnlocalizedFormulas(const std::unordered_map<ChemicalFormula*, int> &value)
		{
			privateCountOfUnlocalizedFormulas = value;
		}

		std::wstring ModificationAnalysisResults::ToString()
		{
			auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->appendLine(MetaMorpheusEngineResults::ToString());
			sb->appendLine(L"Localized mods seen below q-value 0.01:");
			sb->appendLine(std::wstring::Join(L"\r\n", getCountOfModsSeenAndLocalized().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return L"\t" + b::Key + L"\t" + b->Value;
			})));
			sb->appendLine(L"(Approx) Additional localized but protein ambiguous mods seen below q-value 0.01:");
			sb->appendLine(std::wstring::Join(L"\r\n", getCountOfAmbiguousButLocalizedModsSeen().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return L"\t" + b::Key + L"\t" + b->Value;
			})));
			sb->appendLine(L"(Approx) Additional unlocalized mods seen below q-value 0.01:");
			sb->appendLine(std::wstring::Join(L"\r\n", getCountOfUnlocalizedMods().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return L"\t" + b::Key + L"\t" + b->Value;
			})));
			sb->appendLine(L"(Approx) Additional unlocalized modification formulas seen below q-value 0.01:");
			sb->appendLine(std::wstring::Join(L"\r\n", getCountOfUnlocalizedFormulas().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return L"\t" + b::Key->Formula + L"\t" + b->Value;
			})));
			sb->appendLine();
			sb->appendLine(L"All mods in database limited to peptides observed in the results:");
			sb->appendLine(std::wstring::Join(L"\r\n", getCountOfEachModSeenOnProteins().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return L"\t" + b::Key + L"\t" + b->Value;
			})));

			delete sb;
			return sb->toString();
		}
	}
}
