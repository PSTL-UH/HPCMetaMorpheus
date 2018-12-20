#pragma once

#include "../MetaMorpheusTask.h"
#include <string>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <limits>
#include <any>
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class XlSearchParameters; }
namespace TaskLayer { class DbForTask; }
namespace TaskLayer { class FileSpecificParameters; }
namespace TaskLayer { class MyTaskResults; }

using namespace EngineLayer;
using namespace EngineLayer::CrosslinkSearch;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace MzLibUtil;
using namespace EngineLayer::FdrAnalysis;

using namespace Proteomics::Fragmentation;

namespace TaskLayer
{
	class XLSearchTask : public MetaMorpheusTask
	{
	private:
		TaskLayer::XlSearchParameters *privateXlSearchParameters;

	public:
		XLSearchTask();

		TaskLayer::XlSearchParameters *getXlSearchParameters() const;
		void setXlSearchParameters(TaskLayer::XlSearchParameters *value);

	protected:
		MyTaskResults *RunSpecific(const std::wstring &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::wstring> &currentRawFileList, const std::wstring &taskId, std::vector<FileSpecificParameters*> &fileSettingsList) override;

		//Calculate the FDR of single peptide FP/TP
	private:
		void SingleFDRAnalysis(std::vector<CrosslinkSpectralMatch*> &items, std::vector<std::wstring> &taskIds);

		//Calculate the FDR of crosslinked peptide FP/TP
		void DoCrosslinkFdrAnalysis(std::vector<CrosslinkSpectralMatch*> &csms);

		//Generate user defined crosslinker
	public:
		static Crosslinker *GenerateUserDefinedCrosslinker(TaskLayer::XlSearchParameters *xlSearchParameters);


		void WritePsmCrossToTsv(std::vector<CrosslinkSpectralMatch*> &items, const std::wstring &filePath, int writeType);

		void WriteCrosslinkToTxtForPercolator(std::vector<CrosslinkSpectralMatch*> &items, const std::wstring &outputFolder, const std::wstring &fileName, Crosslinker *crosslinker, std::vector<std::wstring> &nestedIds);

		void WritePepXML_xl(std::vector<CrosslinkSpectralMatch*> &items, std::vector<Protein*> &proteinList, const std::wstring &databasePath, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<std::wstring> &localizeableModificationTypes, const std::wstring &outputFolder, const std::wstring &fileName, std::vector<std::wstring> &nestedIds);
	};
}
