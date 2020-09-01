#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <cmath>
#include "exceptionhelper.h"

#include "../PeptideSpectralMatch.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../PrecursorSearchModes/MassDiffAcceptor.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "Proteomics/Proteomics.h"
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
            std::vector<PeptideSpectralMatch*> &PeptideSpectralMatches;
            std::vector<Ms2ScanWithSpecificMass*> const ListOfSortedMs2Scans;
            std::vector<PeptideWithSetModifications*> PeptideIndex;
            const int CurrentPartition;
            MassDiffAcceptor *const massDiffAcceptor;
            const DissociationType dissociationType;
            const double MaxMassThatFragmentIonScoreIsDoubled;
            
        public:

            virtual ~ModernSearchEngine()
            {
                delete massDiffAcceptor;
                //delete dissociationType;
            }
            
            ModernSearchEngine(std::vector<PeptideSpectralMatch*> &globalPsms,
                               std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans,
                               std::vector<PeptideWithSetModifications*> &peptideIndex,
                               std::vector<std::vector<int>> &fragmentIndex,
                               int currentPartition,
                               CommonParameters *commonParameters,
                               MassDiffAcceptor *massDiffAcceptor,
                               double maximumMassThatFragmentIonScoreIsDoubled,
                               std::vector<std::string> &nestedIds,
                               int verbosityLevel=0);

            
        protected:
            MetaMorpheusEngineResults *RunSpecific() override;
            
            std::vector<int> GetBinsToSearch(Ms2ScanWithSpecificMass *scan);
            
            static int BinarySearchBinForPrecursorIndex(std::vector<int> &peptideIdsInThisBin,
                                                        double peptideMassToLookFor,
                                                        std::vector<PeptideWithSetModifications*> &peptideIndex);
            
            void IndexedScoring(std::vector<int> &binsToSearch, std::vector<unsigned char> &scoringTable,
                                unsigned char byteScoreCutoff, std::vector<int> &idsOfPeptidesPossiblyObserved,
                                double scanPrecursorMass, double lowestMassPeptideToLookFor,
                                double highestMassPeptideToLookFor, std::vector<PeptideWithSetModifications*> &peptideIndex,
                                MassDiffAcceptor *massDiffAcceptor, double maxMassThatFragmentIonScoreIsDoubled);
        };
    }
}
