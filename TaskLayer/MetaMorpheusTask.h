#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <any>
#include <typeinfo>
#include "exceptionhelper.h"
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_event.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class CommonParameters; }
namespace TaskLayer { class MyTaskResults; }
namespace EngineLayer { class Ms2ScanWithSpecificMass; }
namespace TaskLayer { class FileSpecificParameters; }
namespace TaskLayer { class DbForTask; }
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class ProgressEventArgs; }
namespace EngineLayer { class SingleEngineFinishedEventArgs; }

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Nett;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{
	enum class MyTask
	{
		Search,
		Gptmd,
		Calibrate,
		XLSearch
	};

	class MetaMorpheusTask
	{
	private:
		MyTask privateTaskType = static_cast<MyTask>(0);
		EngineLayer::CommonParameters *privateCommonParameters;

	public:
		static TomlSettings *const tomlConfig;

	protected:
		StringBuilder *const ProseCreatedWhileRunning = new StringBuilder();

		MyTaskResults *MyTaskResults;

	public:
		virtual ~MetaMorpheusTask()
		{
			delete ProseCreatedWhileRunning;
			delete MyTaskResults;
		}

	protected:
		MetaMorpheusTask(MyTask taskType);

	public:
		static TangibleEvent<EventHandler<SingleTaskEventArgs>> *FinishedSingleTaskHandler = new TangibleEvent<EventHandler<SingleTaskEventArgs>>();

		static TangibleEvent<EventHandler<SingleFileEventArgs>> *FinishedWritingFileHandler = new TangibleEvent<EventHandler<SingleFileEventArgs>>();

		static TangibleEvent<EventHandler<SingleTaskEventArgs>> *StartingSingleTaskHander = new TangibleEvent<EventHandler<SingleTaskEventArgs>>();

		static TangibleEvent<EventHandler<StringEventArgs>> *StartingDataFileHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler<StringEventArgs>> *FinishedDataFileHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler<StringEventArgs>> *OutLabelStatusHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler<StringEventArgs>> *WarnHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler<StringEventArgs>> *LogHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler<StringEventArgs>> *NewCollectionHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler<ProgressEventArgs>> *OutProgressHandler = new TangibleEvent<EventHandler<ProgressEventArgs>>();

		MyTask getTaskType() const;
		void setTaskType(MyTask value);

		EngineLayer::CommonParameters *getCommonParameters() const;
		void setCommonParameters(EngineLayer::CommonParameters *value);

		static const std::wstring IndexFolderName;

		static std::vector<Ms2ScanWithSpecificMass*> GetMs2Scans(MsDataFile *myMSDataFile, const std::wstring &fullFilePath, EngineLayer::CommonParameters *commonParameters);

		static EngineLayer::CommonParameters *SetAllFileSpecificCommonParams(EngineLayer::CommonParameters *commonParams, FileSpecificParameters *fileSpecificParams);

		MyTaskResults *RunTask(const std::wstring &output_folder, std::vector<DbForTask*> &currentProteinDbFilenameList, std::vector<std::wstring> &currentRawDataFilepathList, const std::wstring &displayName);

	protected:
		std::vector<Protein*> LoadProteins(const std::wstring &taskId, std::vector<DbForTask*> &dbFilenameList, bool searchTarget, DecoyType *decoyType, std::vector<std::wstring> &localizeableModificationTypes, EngineLayer::CommonParameters *commonParameters);

		static std::vector<Protein*> LoadProteinDb(const std::wstring &fileName, bool generateTargets, DecoyType *decoyType, std::vector<std::wstring> &localizeableModificationTypes, bool isContaminant, std::unordered_map<std::wstring, Modification*> &um, int &emptyEntriesCount, EngineLayer::CommonParameters *commonParameters);

		void LoadModifications(const std::wstring &taskId, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<std::wstring> &localizableModificationTypes);

		static void WritePsmsToTsv(std::vector<PeptideSpectralMatch*> &psms, const std::wstring &filePath, IReadOnlyDictionary<std::wstring, int> *modstoWritePruned);

		void ReportProgress(ProgressEventArgs *v);

		virtual MyTaskResults *RunSpecific(const std::wstring &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::wstring> &currentRawFileList, const std::wstring &taskId, std::vector<FileSpecificParameters*> &fileSettingsList) = 0;

		void FinishedWritingFile(const std::wstring &path, std::vector<std::wstring> &nestedIDs);

		void StartingDataFile(const std::wstring &v, std::vector<std::wstring> &nestedIDs);

		void FinishedDataFile(const std::wstring &v, std::vector<std::wstring> &nestedIDs);

		void Status(const std::wstring &v, const std::wstring &id);

		void Status(const std::wstring &v, std::vector<std::wstring> &nestedIds);

		static void Warn(const std::wstring &v);

		void Log(const std::wstring &v, std::vector<std::wstring> &nestedIds);

		void NewCollection(const std::wstring &displayName, std::vector<std::wstring> &nestedIds);

	private:
		static std::vector<std::wstring> GetModsTypesFromString(const std::wstring &value);

//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//		private static List<(string, string)> GetModsFromString(string value)
//		{
//			return value.Split(new string[] { "\t\t" }, StringSplitOptions.RemoveEmptyEntries).Select(b => (b.Split('\t').First(), b.Split('\t').Last())).ToList();
//		}

		void SingleEngineHandlerInTask(std::any sender, SingleEngineFinishedEventArgs *e);

		void FinishedSingleTask(const std::wstring &displayName);

		void StartingSingleTask(const std::wstring &displayName);

		static std::vector<std::type_info> GetSubclassesAndItself(std::type_info type);

		static bool SameSettings(const std::wstring &pathToOldParamsFile, IndexingEngine *indexEngine);

		static void WritePeptideIndex(std::vector<PeptideWithSetModifications*> &peptideIndex, const std::wstring &peptideIndexFile);

		static void WriteFragmentIndexNetSerializer(std::vector<std::vector<int>&> &fragmentIndex, const std::wstring &fragmentIndexFile);

		static std::wstring GetExistingFolderWithIndices(IndexingEngine *indexEngine, std::vector<DbForTask*> &dbFilenameList);

		static std::wstring CheckFiles(IndexingEngine *indexEngine, DirectoryInfo *folder);

		static void WriteIndexEngineParams(IndexingEngine *indexEngine, const std::wstring &fileName);

		static std::wstring GenerateOutputFolderForIndices(std::vector<DbForTask*> &dbFilenameList);

	public:
		void GenerateIndexes(IndexingEngine *indexEngine, std::vector<DbForTask*> &dbFilenameList, std::vector<PeptideWithSetModifications*> &peptideIndex, std::vector<std::vector<int>> &fragmentIndex, std::vector<std::vector<int>> &precursorIndex, std::vector<Protein*> &allKnownProteins, std::vector<Modification*> &allKnownModifications, const std::wstring &taskId);
	};
}
