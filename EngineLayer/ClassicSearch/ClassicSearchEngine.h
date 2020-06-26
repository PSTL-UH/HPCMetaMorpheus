#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <any>
#include <mutex>

#include "../PrecursorSearchModes/MassDiffAcceptor.h"
#include "../PeptideSpectralMatch.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "../ScanWithIndexAndNotchInfo.h"

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    namespace ClassicSearch
    {
        class ClassicSearchEngine : public MetaMorpheusEngine
        {
        private:
            MassDiffAcceptor *const SearchMode;
            const std::vector<Protein*> Proteins;
            const std::vector<Modification*> FixedModifications;
            const std::vector<Modification*> VariableModifications;
            std::vector<PeptideSpectralMatch*> &PeptideSpectralMatches;
            std::vector<Ms2ScanWithSpecificMass*> const ArrayOfSortedMS2Scans;
            std::vector<double> MyScanPrecursorMasses;
            
        public:
            virtual ~ClassicSearchEngine()
            {
                //delete SearchMode;
            }
            
            ClassicSearchEngine(std::vector<PeptideSpectralMatch*> &globalPsms,
                                std::vector<Ms2ScanWithSpecificMass*> &arrayOfSortedMS2Scans,
                                std::vector<Modification*> &variableModifications,
                                std::vector<Modification*> &fixedModifications,
                                std::vector<Protein*> &proteinList,
                                MassDiffAcceptor *searchMode,
                                CommonParameters *commonParameters,
                                std::vector<std::string> &nestedIds);
            
        protected:
            MetaMorpheusEngineResults *RunSpecific() override;
            
        private:
            std::vector<ScanWithIndexAndNotchInfo*> GetAcceptableScans(double peptideMonoisotopicMass, MassDiffAcceptor *searchMode);
            
            int GetFirstScanWithMassOverOrEqual(double minimum);
        };
    }
}
