#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <cmath>
#include "exceptionhelper.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class Ms2ScanWithSpecificMass; }
namespace EngineLayer { class MassDiffAcceptor; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace ModernSearch
	{
		class ModernSearchEngine : public MetaMorpheusEngine
		{
		protected:
			static constexpr int FragmentBinsPerDalton = 1000;
			std::vector<std::vector<int>> const FragmentIndex;
			std::vector<PeptideSpectralMatch*> const PeptideSpectralMatches;
			std::vector<Ms2ScanWithSpecificMass*> const ListOfSortedMs2Scans;
			const std::vector<PeptideWithSetModifications*> PeptideIndex;
			const int CurrentPartition;
			MassDiffAcceptor *const MassDiffAcceptor;
			DissociationType *const DissociationType;
			const double MaxMassThatFragmentIonScoreIsDoubled;

		public:
			virtual ~ModernSearchEngine()
			{
				delete MassDiffAcceptor;
				delete DissociationType;
			}

			ModernSearchEngine(std::vector<PeptideSpectralMatch*> &globalPsms, std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans, std::vector<PeptideWithSetModifications*> &peptideIndex, std::vector<std::vector<int>&> &fragmentIndex, int currentPartition, CommonParameters *commonParameters, MassDiffAcceptor *massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled, std::vector<std::wstring> &nestedIds);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;

			std::vector<int> GetBinsToSearch(Ms2ScanWithSpecificMass *scan);

			static int BinarySearchBinForPrecursorIndex(std::vector<int> &peptideIdsInThisBin, double peptideMassToLookFor, std::vector<PeptideWithSetModifications*> &peptideIndex);

			void IndexedScoring(std::vector<int> &binsToSearch, std::vector<unsigned char> &scoringTable, unsigned char byteScoreCutoff, std::vector<int> &idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor, std::vector<PeptideWithSetModifications*> &peptideIndex, MassDiffAcceptor *massDiffAcceptor, double maxMassThatFragmentIonScoreIsDoubled);
		};
	}
}
