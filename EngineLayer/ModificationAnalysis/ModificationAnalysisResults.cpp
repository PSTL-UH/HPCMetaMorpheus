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

		std::unordered_map<std::string, int> ModificationAnalysisResults::getCountOfAmbiguousButLocalizedModsSeen() const
		{
			return privateCountOfAmbiguousButLocalizedModsSeen;
		}

		void ModificationAnalysisResults::setCountOfAmbiguousButLocalizedModsSeen(const std::unordered_map<std::string, int> &value)
		{
			privateCountOfAmbiguousButLocalizedModsSeen = value;
		}

		std::unordered_map<std::string, int> ModificationAnalysisResults::getCountOfModsSeenAndLocalized() const
		{
			return privateCountOfModsSeenAndLocalized;
		}

		void ModificationAnalysisResults::setCountOfModsSeenAndLocalized(const std::unordered_map<std::string, int> &value)
		{
			privateCountOfModsSeenAndLocalized = value;
		}

		std::unordered_map<std::string, int> ModificationAnalysisResults::getCountOfEachModSeenOnProteins() const
		{
			return privateCountOfEachModSeenOnProteins;
		}

		void ModificationAnalysisResults::setCountOfEachModSeenOnProteins(const std::unordered_map<std::string, int> &value)
		{
			privateCountOfEachModSeenOnProteins = value;
		}

		std::unordered_map<std::string, int> ModificationAnalysisResults::getCountOfUnlocalizedMods() const
		{
			return privateCountOfUnlocalizedMods;
		}

		void ModificationAnalysisResults::setCountOfUnlocalizedMods(const std::unordered_map<std::string, int> &value)
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

		std::string ModificationAnalysisResults::ToString()
		{
			auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->appendLine(MetaMorpheusEngineResults::ToString());
			sb->appendLine("Localized mods seen below q-value 0.01:");
			sb->appendLine(std::string::Join("\r\n", getCountOfModsSeenAndLocalized().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return "\t" + b::Key + "\t" + b->Value;
			})));
			sb->appendLine("(Approx) Additional localized but protein ambiguous mods seen below q-value 0.01:");
			sb->appendLine(std::string::Join("\r\n", getCountOfAmbiguousButLocalizedModsSeen().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return "\t" + b::Key + "\t" + b->Value;
			})));
			sb->appendLine("(Approx) Additional unlocalized mods seen below q-value 0.01:");
			sb->appendLine(std::string::Join("\r\n", getCountOfUnlocalizedMods().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return "\t" + b::Key + "\t" + b->Value;
			})));
			sb->appendLine("(Approx) Additional unlocalized modification formulas seen below q-value 0.01:");
			sb->appendLine(std::string::Join("\r\n", getCountOfUnlocalizedFormulas().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return "\t" + b::Key->Formula + "\t" + b->Value;
			})));
			sb->appendLine();
			sb->appendLine("All mods in database limited to peptides observed in the results:");
			sb->appendLine(std::string::Join("\r\n", getCountOfEachModSeenOnProteins().OrderBy([&] (std::any b)
			{
				-b->Value;
			})->Select([&] (std::any b)
			{
			delete sb;
				return "\t" + b::Key + "\t" + b->Value;
			})));

			delete sb;
			return sb->toString();
		}
	}
}
