#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <tuple>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { namespace CrosslinkSearch { class Crosslinker; } }

using namespace MassSpectrometry;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace CrosslinkSearch
	{
		class CrosslinkedPeptide
		{
		public:
			static std::vector<std::tuple<int, std::vector<Product*>>> XlGetTheoreticalFragments(DissociationType *dissociationType, Crosslinker *crosslinker, std::vector<int> &possibleCrosslinkerPositions, double otherPeptideMass, PeptideWithSetModifications *peptide);

			static std::unordered_map<std::tuple<int, int>, std::vector<Product*>> XlLoopGetTheoreticalFragments(DissociationType *dissociationType, Modification *loopMass, std::vector<int> &modPos, PeptideWithSetModifications *peptide);
		};
	}
}
