#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "tangible_filesystem.h"

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

#include "../EngineLayer/EngineLayer.h"
using namespace EngineLayer;

#include "../EngineLayer/CrosslinkSearch/CrosslinkSearchEngine.h"
using namespace EngineLayer::CrosslinkSearch;

#include "../EngineLayer/Indexing/IndexingResults.h"
using namespace EngineLayer::Indexing;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::AminoAcidPolymer;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

#include "../TaskLayer/TaskLayer.h"
#include "../TaskLayer/XLSearchTask/XLSearchParameters.h"
using namespace TaskLayer;

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
using namespace UsefulProteomicsDatabases;

namespace Test
{
    class CSMSerializationTest final
    {
    public:
        static void CSMSerializationTest_BSA_DSSO();
        static void TestDeadendTrisSerialized();
    };

    class XLTestDataFile : public MsDataFile
    {
        //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
    public:
        XLTestDataFile();
        
        std::string getFilePath() const;
        
        std::string getName() const;
        
        void ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities);
    };

}
