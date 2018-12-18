#pragma once

#include "../MetaMorpheusTask.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <limits>
#include <any>
#include <tuple>
#include <optional>
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class GptmdParameters; }
namespace TaskLayer { class DbForTask; }
namespace TaskLayer { class FileSpecificParameters; }
namespace TaskLayer { class MyTaskResults; }

using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::Gptmd;

#if defined(NETFRAMEWORK)

using namespace IO::Thermo;

#else
#endif

using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace UsefulProteomicsDatabases;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{
	class GptmdTask : public MetaMorpheusTask
	{
	private:
		EngineLayer::GptmdParameters *privateGptmdParameters;

		static constexpr double tolForComboLoading = 1e-3;

	public:
		GptmdTask();

		EngineLayer::GptmdParameters *getGptmdParameters() const;
		void setGptmdParameters(EngineLayer::GptmdParameters *value);

	protected:
		MyTaskResults *RunSpecific(const std::wstring &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::wstring> &currentRawFileList, const std::wstring &taskId, std::vector<FileSpecificParameters*> &fileSettingsList) override;

	private:
		static std::vector<std::tuple<double, double>> LoadCombos(std::vector<Modification*> &modificationsThatCanBeCombined);

		static std::vector<double> GetAcceptableMassShifts(std::vector<Modification*> &fixedMods, std::vector<Modification*> &variableMods, std::vector<Modification*> &gptmdMods, std::vector<std::tuple<double, double>> &combos);

		static std::vector<double> GetObservedMasses(std::vector<Modification*> &enumerable, std::vector<Modification*> &gptmdModifications);
	};
}
