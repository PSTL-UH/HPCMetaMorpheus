#pragma once

#include "../MetaMorpheusEngineResults.h"
#include <string>
#include <unordered_map>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { namespace ModificationAnalysis { class ModificationAnalysisEngine; } }

using namespace Chemistry;

namespace EngineLayer
{
	namespace ModificationAnalysis
	{
		class ModificationAnalysisResults : public MetaMorpheusEngineResults
		{
		private:
			std::unordered_map<std::wstring, int> privateCountOfAmbiguousButLocalizedModsSeen;
			std::unordered_map<std::wstring, int> privateCountOfModsSeenAndLocalized;
			std::unordered_map<std::wstring, int> privateCountOfEachModSeenOnProteins;
			std::unordered_map<std::wstring, int> privateCountOfUnlocalizedMods;
			std::unordered_map<ChemicalFormula*, int> privateCountOfUnlocalizedFormulas;

		public:
			ModificationAnalysisResults(ModificationAnalysisEngine *modificationAnalysisEngine);

			/// <summary>
			/// String is the mod ID, integer is the count of that mod observed
			/// </summary>
				std::unordered_map<std::wstring, int> getCountOfAmbiguousButLocalizedModsSeen() const;
				void setCountOfAmbiguousButLocalizedModsSeen(const std::unordered_map<std::wstring, int> &value);
				std::unordered_map<std::wstring, int> getCountOfModsSeenAndLocalized() const;
				void setCountOfModsSeenAndLocalized(const std::unordered_map<std::wstring, int> &value);
				std::unordered_map<std::wstring, int> getCountOfEachModSeenOnProteins() const;
				void setCountOfEachModSeenOnProteins(const std::unordered_map<std::wstring, int> &value);
				std::unordered_map<std::wstring, int> getCountOfUnlocalizedMods() const;
				void setCountOfUnlocalizedMods(const std::unordered_map<std::wstring, int> &value);
				std::unordered_map<ChemicalFormula*, int> getCountOfUnlocalizedFormulas() const;
				void setCountOfUnlocalizedFormulas(const std::unordered_map<ChemicalFormula*, int> &value);

			std::wstring ToString() override;
		};
	}
}
