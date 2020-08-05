#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <limits>
#include <tuple>

#include "CrosslinkSpectralMatch.h"
#include "Crosslinker.h"
#include "../PrecursorSearchModes/MassDiffAcceptor.h"

#include "../CommonParameters.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../MetaMorpheusEngineResults.h"
#include "BestPeptideScoreNotch.h"

#include "../ModernSearch/ModernSearchEngine.h"
using namespace EngineLayer::ModernSearch;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    namespace CrosslinkSearch
    {
        class CrosslinkSearchEngine : public ModernSearchEngine
        {
            //protected:
        public:
            std::vector<CrosslinkSpectralMatch*> &GlobalCsms;
            
            // crosslinker molecule
        private:
            Crosslinker *const privateCrosslinker;
            
            const bool CrosslinkSearchTopN;
            const int TopN;
            const bool QuenchH2O;
            const bool QuenchNH2;
            const bool QuenchTris;
            MassDiffAcceptor *XLPrecusorSearchMode=nullptr;
            Modification *TrisDeadEnd=nullptr;
            Modification *H2ODeadEnd=nullptr;
            Modification *NH2DeadEnd=nullptr;
            Modification *Loop=nullptr;
            std::vector<char> AllCrosslinkerSites;
            
        public:
            virtual ~CrosslinkSearchEngine()
                {
                    //delete Crosslinker;
                    //delete XLPrecusorSearchMode;
                    //delete TrisDeadEnd;
                    //delete H2ODeadEnd;
                    //delete NH2DeadEnd;
                    if ( Loop != nullptr ) {
                        delete Loop;
                    }
                }
            
            CrosslinkSearchEngine(std::vector<CrosslinkSpectralMatch*> &globalCsms,
                                  std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans,
                                  std::vector<PeptideWithSetModifications*> &peptideIndex,
                                  std::vector<std::vector<int>> &fragmentIndex,
                                  int currentPartition,
                                  CommonParameters *commonParameters,
                                  Crosslinker *crosslinker,
                                  bool CrosslinkSearchTop,
                                  int CrosslinkSearchTopNum,
                                  bool quench_H2O,
                                  bool quench_NH2,
                                  bool quench_Tris,
                                  std::vector<std::string> &nestedIds);
            
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
