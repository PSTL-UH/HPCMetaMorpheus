#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <typeinfo>
#include <fstream>
#include <tuple>

#include "exceptionhelper.h"
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_filesystem.h"

#include "EventArgs/SingleTaskEventArgs.h"
#include "DbForTask.h"

#include "MyTaskResults.h"
#include "MyFileManager.h"
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

#include "mpi.h"

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
            int privateVerbosityLevel = 0;
            
	public:
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
            
            MyTask getTaskType() const;
            void setTaskType(MyTask value);
            
            EngineLayer::CommonParameters *getCommonParameters() const;
            void setCommonParameters(EngineLayer::CommonParameters *value);

            int getVerbose() const;
            void setVerbose ( int verbosityLevel );
            
            static const std::string IndexFolderName;
            
            static std::vector<Ms2ScanWithSpecificMass*> GetMs2Scans(MsDataFile *myMSDataFile,
                                                                     const std::string &fullFilePath,
                                                                     EngineLayer::CommonParameters *commonParameters,
                                                                     int firstIndex=-1, int lastIndex=-1 );
            
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

            
            virtual MyTaskResults *RunSpecific(const std::string &OutputFolder, std::vector<DbForTask*> &dbFilenameList,
                                               std::vector<std::string> &currentRawFileList, const std::string &taskId,
                                               std::vector<FileSpecificParameters*> &fileSettingsList) = 0;

        public:
            static void ReportProgress(ProgressEventArgs *v,  int verbositLevel);
            static void ReportEngineProgress(std::string id, int value, int verbosityLevel);
            static void FinishedWritingFile(const std::string &path, std::vector<std::string> &nestedIDs,  int verbositLevel);            
            static void StartingDataFile(const std::string &v, std::vector<std::string> &nestedIDs,  int verbositLevel);
            static void FinishedDataFile(const std::string &v, std::vector<std::string> &nestedIDs,  int verbositLevel);
            static void StartingSingleTask(const std::string &taskName, int verbosityLevel);
            static void FinishedSingleTask(const std::string &displayName, int verbosityLevel);
            static void StartingSingleEngine(std::vector<std::string> &nestedIDs, int verbosityLevel);
            static void FinishedSingleEngine(std::vector<std::string> &nestedIDs, MetaMorpheusEngineResults *myResults, int verbosityLevel);
            static void Status(const std::string &v, const std::string &id,  int verbositLevel);
            static void Status(const std::string &v, std::vector<std::string> &nestedIds,  int verbositLevel);
            static void Warn(const std::string &v);
            static void Log(const std::string &v, std::vector<std::string> &nestedIds,  int verbositLevel);
            
            //static void NewCollection(const std::string &displayName, std::vector<std::string> &nestedIds);
            
	private:
            static std::vector<std::string> GetModsTypesFromString(const std::string &value);
            
            //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
            // private static List<(string, string)> GetModsFromString(string value)
            // {
            //      return value.Split(new string[] { "\t\t" }, StringSplitOptions.RemoveEmptyEntries).Select(
            //                                    b => (b.Split('\t').First(), b.Split('\t').Last())).ToList();
            // }
            static std::pair<std::string, std::string> GetModsFromString (std::string value);
                        
            static bool SameSettings(const std::string &pathToOldParamsFile, IndexingEngine *indexEngine);
            
            static void WritePeptideIndex(std::vector<PeptideWithSetModifications*> &peptideIndex,
                                          std::string &peptideIndexFile);
            
            static void WriteFragmentIndexSerializer(std::vector<std::vector<int>> &fragmentIndex,
                                                     const std::string &fragmentIndexFile);
            static void ReadFragmentIndexDeserializer(std::vector<std::vector<int>> &fragmentIndex,
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
                                 const std::string &taskId,
                                 MPI_Comm comm=MPI_COMM_WORLD);

            void DataFilePartitioning ( std::vector<std::string> &allFiles,
                                        MPI_Comm comm,
                                        MyFileManager *fileManager,                             
                                        std::vector<FileSpecificParameters*> &fileSettingsList, 
                                        std::vector<std::string> &myFile,                     
                                        std::vector<std::tuple<int, int>> &myFirstLastIndex );
        };
}
