#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <limits>
#include <tuple>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { namespace CrosslinkSearch { class CrosslinkSpectralMatch; } }
namespace EngineLayer { namespace CrosslinkSearch { class Crosslinker; } }
namespace EngineLayer { class MassDiffAcceptor; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class Ms2ScanWithSpecificMass; }
namespace EngineLayer { class MetaMorpheusEngineResults; }
namespace EngineLayer { namespace CrosslinkSearch { class BestPeptideScoreNotch; } }

using namespace EngineLayer::ModernSearch;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace CrosslinkSearch
	{
		class CrosslinkSearchEngine : public ModernSearchEngine
		{
		protected:
			std::vector<CrosslinkSpectralMatch*> const GlobalCsms;

			// crosslinker molecule
		private:
			Crosslinker *const Crosslinker;

			const bool CrosslinkSearchTopN;
			const int TopN;
			const bool QuenchH2O;
			const bool QuenchNH2;
			const bool QuenchTris;
			MassDiffAcceptor *XLPrecusorSearchMode;
			Modification *TrisDeadEnd;
			Modification *H2ODeadEnd;
			Modification *NH2DeadEnd;
			Modification *Loop;
			std::vector<wchar_t> AllCrosslinkerSites;

		public:
			virtual ~CrosslinkSearchEngine()
			{
				delete Crosslinker;
				delete XLPrecusorSearchMode;
				delete TrisDeadEnd;
				delete H2ODeadEnd;
				delete NH2DeadEnd;
				delete Loop;
			}

			CrosslinkSearchEngine(std::vector<CrosslinkSpectralMatch*> &globalCsms, std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans, std::vector<PeptideWithSetModifications*> &peptideIndex, std::vector<std::vector<int>&> &fragmentIndex, int currentPartition, CommonParameters *commonParameters, Crosslinker *crosslinker, bool CrosslinkSearchTop, int CrosslinkSearchTopNum, bool quench_H2O, bool quench_NH2, bool quench_Tris, std::vector<std::string> &nestedIds);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;

			/// <summary>
			/// 
			/// </summary>
		private:
			void GenerateCrosslinkModifications(Crosslinker *crosslinker);

			/// <summary>
			/// 
			/// </summary>
			CrosslinkSpectralMatch *FindCrosslinkedPeptide(Ms2ScanWithSpecificMass *theScan, std::vector<BestPeptideScoreNotch*> &theScanBestPeptide, int scanIndex);

			/// <summary>
			/// Localizes the crosslink position on the alpha and beta peptides
			/// </summary>
			CrosslinkSpectralMatch *LocalizeCrosslinkSites(Ms2ScanWithSpecificMass *theScan, BestPeptideScoreNotch *alphaPeptide, BestPeptideScoreNotch *betaPeptide, Crosslinker *crosslinker, int ind, int inx);

			/// <summary>
			/// Localizes the deadend mod to a residue
			/// </summary>
			CrosslinkSpectralMatch *LocalizeDeadEndSite(PeptideWithSetModifications *originalPeptide, Ms2ScanWithSpecificMass *theScan, CommonParameters *commonParameters, std::vector<int> &possiblePositions, Modification *deadEndMod, int notch, int scanIndex, int peptideIndex);

			/// <summary>
			/// Localizes the loop to a begin and end residue
			/// </summary>
			CrosslinkSpectralMatch *LocalizeLoopSites(PeptideWithSetModifications *originalPeptide, Ms2ScanWithSpecificMass *theScan, CommonParameters *commonParameters, std::vector<int> &possiblePositions, Modification *loopMod, int notch, int scanIndex, int peptideIndex);
		};
	}
}
