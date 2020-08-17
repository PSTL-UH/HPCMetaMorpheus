#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <any>
#include <typeinfo>
#include <fstream>

#include "exceptionhelper.h"
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_filesystem.h"

#include "EventArgs/SingleTaskEventArgs.h"
#include "DbForTask.h"

#include "MyTaskResults.h"
#include "FileSpecificParameters.h"

#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineEventArgs.h"
#include "../EngineLayer/EventArgs/SingleFileEventArgs.h"

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

using namespace EngineLayer;

#include "../EngineLayer/Indexing/IndexingEngine.h"
using namespace EngineLayer::Indexing;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{
	enum class MyTask
	{
            Search = 0,
	    Gptmd,
	    Calibrate,
	    XLSearch
        };

	class MetaMorpheusTask
	{
	private:
            MyTask privateTaskType = static_cast<MyTask>(0);
            EngineLayer::CommonParameters *privateCommonParameters = nullptr;
            bool privateVerbose = false;
            
	public:
            //static TomlSettings *const tomlConfig;
            // static toml::Table tomlConfig;
            std::ofstream tomlFile;
            virtual void writeTomlConfig ( std::string &filename, std::ofstream &tomlFd);
            
	protected:
            StringBuilder *const ProseCreatedWhileRunning = new StringBuilder();
            MyTaskResults *myTaskResults=nullptr;
            
	public:
            virtual ~MetaMorpheusTask()
            {
                delete ProseCreatedWhileRunning;
                if (myTaskResults != nullptr ) {
                    delete myTaskResults;
                }
                if ( privateCommonParameters != nullptr ) {
                    delete privateCommonParameters;
                }
            }
            
	protected:
            MetaMorpheusTask(MyTask taskType);
            
	public:
            // Keeping one original line for revisiting the item later.
            //static TangibleEvent<EventHandler<SingleTaskEventArgs>> *FinishedSingleTaskHandler = new TangibleEvent<EventHandler<SingleTaskEventArgs>>();
            
#ifdef ORIG
            static EventHandler<SingleTaskEventArgs> *FinishedSingleTaskHandler = new EventHandler<SingleTaskEventArgs>();
            static EventHandler<SingleFileEventArgs> *FinishedWritingFileHandler = new EventHandler<SingleFileEventArgs>();
            static EventHandler<SingleTaskEventArgs> *StartingSingleTaskHandler = new EventHandler<SingleTaskEventArgs>();
            static EventHandler<StringEventArgs> *StartingDataFileHandler = new EventHandler<StringEventArgs>();
            static EventHandler<StringEventArgs> *FinishedDataFileHandler = new EventHandler<StringEventArgs>();
            static EventHandler<StringEventArgs> *OutLabelStatusHandler = new TangibleEvent<EventHandler<StringEventArgs>();
            static EventHandler<StringEventArgs> *WarnHandler = new EventHandler<StringEventArgs>>();
            static EventHandler<StringEventArgs> *LogHandler = new EventHandler<StringEventArgs>>();
            static EventHandler<StringEventArgs> *NewCollectionHandler = new EventHandler<StringEventArgs>();
            static EventHandler<ProgressEventArgs> *OutProgressHandler = new EventHandler<ProgressEventArgs>();
#endif

            EventHandler<SingleTaskEventArgs> *FinishedSingleTaskHandler=nullptr;
            EventHandler<SingleFileEventArgs> *FinishedWritingFileHandler=nullptr;
            EventHandler<SingleTaskEventArgs> *StartingSingleTaskHandler=nullptr;
            EventHandler<StringEventArgs> *StartingDataFileHandler=nullptr;
            EventHandler<StringEventArgs> *FinishedDataFileHandler=nullptr;
            EventHandler<StringEventArgs> *OutLabelStatusHandler=nullptr;
            EventHandler<StringEventArgs> *WarnHandler=nullptr;
            EventHandler<StringEventArgs> *LogHandler=nullptr;
            EventHandler<StringEventArgs> *NewCollectionHandler=nullptr;
            EventHandler<ProgressEventArgs> *OutProgressHandler=nullptr;
            
            MyTask getTaskType() const;
            void setTaskType(MyTask value);
            
            EngineLayer::CommonParameters *getCommonParameters() const;
            void setCommonParameters(EngineLayer::CommonParameters *value);

            bool getVerbose() const;
            void setVerbose ( bool verbose );
            
            static const std::string IndexFolderName;
            
            static std::vector<Ms2ScanWithSpecificMass*> GetMs2Scans(MsDataFile *myMSDataFile,
                                                                     const std::string &fullFilePath,
                                                                     EngineLayer::CommonParameters *commonParameters);
            
            static EngineLayer::CommonParameters *SetAllFileSpecificCommonParams(EngineLayer::CommonParameters *commonParams,
                                                                                 FileSpecificParameters *fileSpecificParams);
            
            MyTaskResults *RunTask(const std::string &output_folder, std::vector<DbForTask*> &currentProteinDbFilenameList,
                                   std::vector<std::string> &currentRawDataFilepathList, const std::string &displayName);

	protected:
            std::vector<Protein*> LoadProteins(const std::string &taskId, std::vector<DbForTask*> &dbFilenameList,
                                               bool searchTarget, DecoyType decoyType,
                                               std::vector<std::string> &localizeableModificationTypes,
                                               EngineLayer::CommonParameters *commonParameters);
            
            static std::vector<Protein*> LoadProteinDb(const std::string &fileName, bool generateTargets,
                                                       DecoyType decoyType,
                                                       std::vector<std::string> &localizeableModificationTypes,
                                                       bool isContaminant, std::unordered_map<std::string, Modification*> &um,
                                                       int &emptyEntriesCount,
                                                       EngineLayer::CommonParameters *commonParameters);

            void LoadModifications(const std::string &taskId, std::vector<Modification*> &variableModifications,
                                   std::vector<Modification*> &fixedModifications,
                                   std::vector<std::string> &localizableModificationTypes);

            static void WritePsmsToTsv(std::vector<PeptideSpectralMatch*> &psms, const std::string &filePath,
                                       std::unordered_map<std::string, int> *modstoWritePruned);

            void ReportProgress(ProgressEventArgs *v);
            
            virtual MyTaskResults *RunSpecific(const std::string &OutputFolder, std::vector<DbForTask*> &dbFilenameList,
                                               std::vector<std::string> &currentRawFileList, const std::string &taskId,
                                               std::vector<FileSpecificParameters*> &fileSettingsList) = 0;

            void FinishedWritingFile(const std::string &path, std::vector<std::string> &nestedIDs);
            
            void StartingDataFile(const std::string &v, std::vector<std::string> &nestedIDs);
            
            void FinishedDataFile(const std::string &v, std::vector<std::string> &nestedIDs);
            
            void Status(const std::string &v, const std::string &id);
            
            void Status(const std::string &v, std::vector<std::string> &nestedIds);
            
            void Warn(const std::string &v);
            
            void Log(const std::string &v, std::vector<std::string> &nestedIds);
            
            void NewCollection(const std::string &displayName, std::vector<std::string> &nestedIds);
            
	private:
            static std::vector<std::string> GetModsTypesFromString(const std::string &value);
            
            //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
            // private static List<(string, string)> GetModsFromString(string value)
            // {
            //      return value.Split(new string[] { "\t\t" }, StringSplitOptions.RemoveEmptyEntries).Select(
            //                                    b => (b.Split('\t').First(), b.Split('\t').Last())).ToList();
            // }
            static std::pair<std::string, std::string> GetModsFromString (std::string value);
            
            //void SingleEngineHandlerInTask(std::any sender, SingleEngineFinishedEventArgs *e);
            void SingleEngineHandlerInTask(SingleEngineFinishedEventArgs e);
            
            void FinishedSingleTask(const std::string &displayName);
            
            void StartingSingleTask(const std::string &displayName);
            
#ifdef NOT_NOW
            static std::vector<std::type_info> GetSubclassesAndItself(std::type_info type);
#endif
            static bool SameSettings(const std::string &pathToOldParamsFile, IndexingEngine *indexEngine);
            
            static void WritePeptideIndex(std::vector<PeptideWithSetModifications*> &peptideIndex,
                                          const std::string &peptideIndexFile);
            
            static void WriteFragmentIndexNetSerializer(std::vector<std::vector<int>> &fragmentIndex,
                                                        const std::string &fragmentIndexFile);
            
            static std::string GetExistingFolderWithIndices(IndexingEngine *indexEngine, std::vector<DbForTask*> &dbFilenameList);
            
            static std::string CheckFiles(IndexingEngine *indexEngine, std::string &folder);
            
            static void WriteIndexEngineParams(IndexingEngine *indexEngine, const std::string &fileName);
            
            static std::string GenerateOutputFolderForIndices(std::vector<DbForTask*> &dbFilenameList);
            
	public:
            void GenerateIndexes(IndexingEngine *indexEngine, std::vector<DbForTask*> &dbFilenameList,
                                 std::vector<PeptideWithSetModifications*> &peptideIndex,
                                 std::vector<std::vector<int>> &fragmentIndex,
                                 std::vector<std::vector<int>> &precursorIndex,
                                 std::vector<Protein*> &allKnownProteins,
                                 std::vector<Modification*> &allKnownModifications,
                                 const std::string &taskId);
	};
}
