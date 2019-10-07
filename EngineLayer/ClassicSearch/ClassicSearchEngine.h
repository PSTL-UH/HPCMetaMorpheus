#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <any>
#include <mutex>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MassDiffAcceptor; }
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class Ms2ScanWithSpecificMass; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }
namespace EngineLayer { class ScanWithIndexAndNotchInfo; }

using namespace MassSpectrometry;
using namespace MzLibUtil;
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
            std::vector<PeptideSpectralMatch*> const PeptideSpectralMatches;
            std::vector<Ms2ScanWithSpecificMass*> const ArrayOfSortedMS2Scans;
            std::vector<double> const MyScanPrecursorMasses;
            
        public:
            virtual ~ClassicSearchEngine()
            {
                //delete SearchMode;
            }
            
            ClassicSearchEngine(std::vector<PeptideSpectralMatch*> &globalPsms, std::vector<Ms2ScanWithSpecificMass*> &arrayOfSortedMS2Scans, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<Protein*> &proteinList, MassDiffAcceptor *searchMode, CommonParameters *commonParameters, std::vector<std::string> &nestedIds);
            
        protected:
            MetaMorpheusEngineResults *RunSpecific() override;
            
        private:
            std::vector<ScanWithIndexAndNotchInfo*> GetAcceptableScans(double peptideMonoisotopicMass, MassDiffAcceptor *searchMode);
            
            int GetFirstScanWithMassOverOrEqual(double minimum);
        };
    }
}
