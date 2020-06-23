#pragma once

#include "../MetaMorpheusTask.h"
#include "MassDiffAcceptorType.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <limits>
#include <stdexcept>
#include <any>
#include <mutex>
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_filesystem.h"

#include "SearchParameters.h"
#include "../../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"
#include "../MyTaskResults.h"

using namespace EngineLayer;

#include "../../EngineLayer/ClassicSearch/ClassicSearchEngine.h"
using namespace EngineLayer::ClassicSearch;

#include "../../EngineLayer/NonSpecificEnzymeSearch/NonSpecificEnzymeSearchEngine.h"
using namespace EngineLayer::NonSpecificEnzymeSearch;

#include "../../EngineLayer/Indexing/IndexingEngine.h"
using namespace EngineLayer::Indexing;

#include "../../EngineLayer/ModernSearch/ModernSearchEngine.h"
using namespace EngineLayer::ModernSearch;


#include "FlashLFQ/FlashLFQ.h"
using namespace FlashLFQ;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include  "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{
    class SearchTask : public MetaMorpheusTask
    {
    private:
        TaskLayer::SearchParameters *privateSearchParameters;
        
    public:
        SearchTask();
        
        TaskLayer::SearchParameters *getSearchParameters() const;
        void setSearchParameters(TaskLayer::SearchParameters *value);

        void writeTomlConfig(std::string &filename, std::ofstream &tomlFd );
        
        static MassDiffAcceptor *GetMassDiffAcceptor(Tolerance *precursorMassTolerance,
                                                     MassDiffAcceptorType massDiffAcceptorType,
                                                     const std::string &customMdac);
        
    protected:
        MyTaskResults *RunSpecific(const std::string &OutputFolder, std::vector<DbForTask*> &dbFilenameList,
                                   std::vector<std::string> &currentRawFileList, const std::string &taskId,
                                   std::vector<FileSpecificParameters*> &fileSettingsList) override;
        
    private:
        int GetNumNotches(MassDiffAcceptorType massDiffAcceptorType, const std::string &customMdac);
        
        static MassDiffAcceptor *ParseSearchMode(const std::string &text);
    };
}
