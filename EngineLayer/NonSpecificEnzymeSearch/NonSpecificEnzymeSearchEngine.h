#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../PrecursorSearchModes/MassDiffAcceptor.h"
#include "../MetaMorpheusEngineResults.h"

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

#include "../FdrAnalysis/FdrAnalysisEngine.h"
using namespace EngineLayer::FdrAnalysis;

#include "../ModernSearch/ModernSearchEngine.h"
using namespace EngineLayer::ModernSearch;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    namespace NonSpecificEnzymeSearch
    {
        class NonSpecificEnzymeSearchEngine : public ModernSearchEngine
        {
        private:
            static const double WaterMonoisotopicMass;
            
            std::vector<std::vector<int>> const PrecursorIndex;
            const int MinimumPeptideLength;
            std::vector<std::vector<PeptideSpectralMatch*>> GlobalCategorySpecificPsms;
            CommonParameters *ModifiedParametersNoComp;
            std::vector<ProductType> ProductTypesToSearch;
            std::vector<PeptideSpectralMatch*> unusedPsms;  // Just to keep the ModernSearchEngine constructor happy

        public:
            virtual ~NonSpecificEnzymeSearchEngine()
            {
                delete ModifiedParametersNoComp;
            }
            
            NonSpecificEnzymeSearchEngine(std::vector<std::vector<PeptideSpectralMatch*>> &globalPsms,
                                          std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans,
                                          std::vector<PeptideWithSetModifications*> &peptideIndex,
                                          std::vector<std::vector<int>> &fragmentIndex,
                                          std::vector<std::vector<int>> &precursorIndex, int currentPartition,
                                          CommonParameters *CommonParameters,
                                          MassDiffAcceptor *massDiffAcceptor,
                                          double maximumMassThatFragmentIonScoreIsDoubled,
                                          std::vector<std::string> &nestedIds);
            
        protected:
            MetaMorpheusEngineResults *RunSpecific() override;
            
        private:
            std::tuple<int, PeptideWithSetModifications*> Accepts(std::vector<Product*> &fragments,
                                                                  double scanPrecursorMass,
                                                                  PeptideWithSetModifications *peptide,
                                                                  FragmentationTerminus fragmentationTerminus,
                                                                  MassDiffAcceptor *searchMode);
            
        public:
            static std::vector<PeptideSpectralMatch*> ResolveFdrCategorySpecificPsms(
                std::vector<std::vector<PeptideSpectralMatch*>> &AllPsms, int numNotches,
                const std::string &taskId, CommonParameters *commonParameters);
        };
    }
}
