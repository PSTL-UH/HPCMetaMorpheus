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

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class SearchParameters; }
namespace EngineLayer { class MassDiffAcceptor; }
namespace TaskLayer { class DbForTask; }
namespace TaskLayer { class FileSpecificParameters; }
namespace TaskLayer { class MyTaskResults; }

using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::Indexing;
using namespace EngineLayer::ModernSearch;
using namespace EngineLayer::NonSpecificEnzymeSearch;
using namespace FlashLFQ;
using namespace MassSpectrometry;
using namespace MzLibUtil;
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

		static MassDiffAcceptor *GetMassDiffAcceptor(Tolerance *precursorMassTolerance, MassDiffAcceptorType massDiffAcceptorType, const std::wstring &customMdac);

	protected:
		MyTaskResults *RunSpecific(const std::wstring &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::wstring> &currentRawFileList, const std::wstring &taskId, std::vector<FileSpecificParameters*> &fileSettingsList) override;

	private:
		int GetNumNotches(MassDiffAcceptorType massDiffAcceptorType, const std::wstring &customMdac);

		static MassDiffAcceptor *ParseSearchMode(const std::wstring &text);
	};
}
