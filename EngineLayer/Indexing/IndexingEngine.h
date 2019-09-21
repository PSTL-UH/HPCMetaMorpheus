#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <cmath>
#include <mutex>
#include "exceptionhelper.h"
#include "stringhelper.h"
#include "stringbuilder.h"
#include "fileinfo.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
//namespace EngineLayer { class CommonParameters; }
#include "../CommonParameters.h"

//namespace EngineLayer { class MetaMorpheusEngineResults; }
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
            DecoyType *const decoyType;
            const double MaxFragmentSize;

        public:
            const bool GeneratePrecursorIndex;
            const std::vector<FileInfo*> ProteinDatabases;
            
            virtual ~IndexingEngine()
            {
                delete decoyType;
            }
            
            IndexingEngine(std::vector<Protein*> &proteinList, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, int currentPartition, DecoyType *decoyType, CommonParameters *commonParams, double maxFragmentSize, bool generatePrecursorIndex, std::vector<FileInfo*> &proteinDatabases, std::vector<std::string> &nestedIds);
            
            std::string ToString();
            
        protected:
            MetaMorpheusEngineResults *RunSpecific() override;
        };
    }
}
