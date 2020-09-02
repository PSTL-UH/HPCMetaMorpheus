#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <cmath>
#include <mutex>
#include "exceptionhelper.h"
#include "stringhelper.h"
#include "stringbuilder.h"

#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
using namespace UsefulProteomicsDatabases;

namespace EngineLayer
{
    namespace Indexing
    {
        class IndexingEngine : public MetaMorpheusEngine
        {
        private:
            static constexpr int FragmentBinsPerDalton = 1000;
            const std::vector<Protein*> ProteinList;
            const std::vector<Modification*> FixedModifications;
            const std::vector<Modification*> VariableModifications;
            const int CurrentPartition;
            DecoyType const decoyType;
            const double MaxFragmentSize;

        public:
            const bool GeneratePrecursorIndex;
            std::vector<std::string> ProteinDatabases;
            
            virtual ~IndexingEngine()
            {
            }
            
            IndexingEngine(std::vector<Protein*> &proteinList, std::vector<Modification*> &variableModifications,
                           std::vector<Modification*> &fixedModifications, int currentPartition,
                           DecoyType decoyType, CommonParameters *commonParams, double maxFragmentSize,
                           bool generatePrecursorIndex, std::vector<std::string> &proteinDatabases,
                           std::vector<std::string> nestedIds, int verbosityLevel=0);
            
            std::string ToString();
            
        protected:
            MetaMorpheusEngineResults *RunSpecific() override;
        };
    }
}
