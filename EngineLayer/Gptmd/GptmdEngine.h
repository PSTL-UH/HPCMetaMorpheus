#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cctype>
#include <tuple>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }

using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace Gptmd
	{
		class GptmdEngine : public MetaMorpheusEngine
		{
		private:
			const std::vector<PeptideSpectralMatch*> AllIdentifications;
			const std::vector<std::tuple<double, double>> Combos;
			const std::vector<Modification*> GptmdModifications;
			const std::unordered_map<std::string, Tolerance*> FilePathToPrecursorMassTolerance; // this exists because of file-specific tolerances

		public:
			GptmdEngine(std::vector<PeptideSpectralMatch*> &allIdentifications, std::vector<Modification*> &gptmdModifications, std::vector<std::tuple<double, double>> &combos, std::unordered_map<std::string, Tolerance*> &filePathToPrecursorMassTolerance, CommonParameters *commonParameters, std::vector<std::string> &nestedIds);

			static bool ModFits(Modification *attemptToLocalize, Protein *protein, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;

		private:
			static void AddIndexedMod(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>> &modDict, const std::string &proteinAccession, std::tuple<int, Modification*> &indexedMod);

			static std::vector<Modification*> GetPossibleMods(double totalMassToGetTo, std::vector<Modification*> &allMods, std::vector<std::tuple<double, double>> &combos, Tolerance *precursorTolerance, PeptideWithSetModifications *peptideWithSetModifications);
		};
	}
}
