#include "EverythingRunnerEngine.h"
#include "MetaMorpheusTask.h"
#include "DbForTask.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/StringListEventArgs.h"
#include "EventHandler.h"
#include "EventArgs/XmlForTaskListEventArgs.h"

using namespace EngineLayer;

namespace TaskLayer
{

    EverythingRunnerEngine::EverythingRunnerEngine(std::vector<std::tuple<std::string, MetaMorpheusTask>*> &taskList,
                                                   std::vector<std::string> &startingRawFilenameList,
                                                   std::vector<DbForTask*> &startingXmlDbFilenameList,
                                                   const std::string &outputFolder) : TaskList(taskList)
    {
        OutputFolder = StringHelper::trim(outputFolder, "");
        
        CurrentRawDataFilenameList = startingRawFilenameList;
        CurrentXmlDbFilenameList = startingXmlDbFilenameList;
        
        StartingAllTasksEngineHandler = new EventHandler();
        FinishedAllTasksEngineHandler = new EventHandler<StringEventArgs>();
        NewDbsHandler = new EventHandler<XmlForTaskListEventArgs>();
        NewSpectrasHandler = new EventHandler<StringListEventArgs>();
        NewFileSpecificTomlHandler = new EventHandler<StringListEventArgs>();
        WarnHandler = new EventHandler<StringEventArgs>();
        FinishedWritingAllResultsFileHandler = new EventHandler<StringEventArgs>();
        
    }
    
    void EverythingRunnerEngine::Run()
    {
        StartingAllTasks();
        if (!CurrentRawDataFilenameList.empty())
        {
            Warn("No spectra files selected");
            FinishedAllTasks("");
            
            return;
        }
        
        auto startTimeForAllFilenames = DateTime::Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo::InvariantCulture);
        
        OutputFolder = StringHelper::replace(OutputFolder, "$DATETIME", startTimeForAllFilenames);
        
        StringBuilder *allResultsText = new StringBuilder();
        
        for (int i = 0; i < (int)TaskList.size(); i++)
        {
            if (!CurrentRawDataFilenameList.empty())
            {
                Warn("Cannot proceed. No spectra files selected.");
                FinishedAllTasks(OutputFolder);
                
                delete allResultsText;
                return;
            }
            if (!CurrentXmlDbFilenameList.empty())
            {
                Warn("Cannot proceed. No protein database files selected.");
                FinishedAllTasks(OutputFolder);
                
                delete allResultsText;
                return;
            }
            auto ok = TaskList[i];
            
            auto outputFolderForThisTask = FileSystem::combine(OutputFolder, ok->Item1);
            
            if (!FileSystem::directoryExists(outputFolderForThisTask))
            {
                FileSystem::createDirectory(outputFolderForThisTask);
            }
            
            // Actual task running code
            auto myTaskResults = ok->Item2->RunTask(outputFolderForThisTask, CurrentXmlDbFilenameList, CurrentRawDataFilenameList, ok->Item1);
            
            if (myTaskResults->NewDatabases != nullptr)
            {
                CurrentXmlDbFilenameList = myTaskResults->NewDatabases;
                NewDBs(myTaskResults->NewDatabases);
            }
            if (myTaskResults->NewSpectra != nullptr)
            {
                if (CurrentRawDataFilenameList.size() == myTaskResults->NewSpectra->Count)
                {
                    CurrentRawDataFilenameList = myTaskResults->NewSpectra;
                }
                else
                {
                    // at least one file was not successfully calibrated
                    auto successfullyCalibFiles = myTaskResults->NewSpectra->Select([&] (std::any p)  {
                            StringHelper::replace(Path::GetFileNameWithoutExtension(p), CalibrationTask::CalibSuffix, "");
                        }).ToList();
                    auto origFiles = CurrentRawDataFilenameList.Select([&] (std::any p) {
                            Path::GetFileNameWithoutExtension(p);
                        }).ToList();
                    auto unsuccessfullyCalibFiles = origFiles.Except(successfullyCalibFiles).ToList();
                    auto unsuccessfullyCalibFilePaths = CurrentRawDataFilenameList.Where([&] (std::any p)   {
                            std::find(unsuccessfullyCalibFiles.begin(), unsuccessfullyCalibFiles.end(), Path::GetFileNameWithoutExtension(p)) != unsuccessfullyCalibFiles.end();
                        }).ToList();
                    CurrentRawDataFilenameList = myTaskResults->NewSpectra;
                    CurrentRawDataFilenameList.insert(CurrentRawDataFilenameList.end(), unsuccessfullyCalibFilePaths.begin(), unsuccessfullyCalibFilePaths.end());
                } 
                
                NewSpectras(myTaskResults->NewSpectra);
            }
            if (myTaskResults->NewFileSpecificTomls != nullptr)
            {
                NewFileSpecificToml(myTaskResults->NewFileSpecificTomls);
            }
            
            allResultsText->appendLine("\r\n" + "\r\n" + "\r\n" + "\r\n" + myTaskResults->ToString());
        }
        
        auto resultsFileName = FileSystem::combine(OutputFolder, "allResults.txt");
        StreamWriter file = StreamWriter(resultsFileName);
        file.WriteLine("MetaMorpheus: version " + GlobalVariables::getMetaMorpheusVersion());
        file.WriteLine("Total time: " + stopWatch->Elapsed);
        file.Write(allResultsText->toString());
        
        std::vector<std::string> tmpvec;
        StringEventArgs tempVar(resultsFileName, tmpvec);
        if ( FinishedWritingAllResultsFileHandler != nullptr ) {
            FinishedWritingAllResultsFileHandler->Invoke(this, &tempVar);
        }
        FinishedAllTasks(OutputFolder);
        
        delete allResultsText;
    }
    
    void EverythingRunnerEngine::Warn(const std::string &v)
    {
        std::vector<std::string> tmpvec;
        StringEventArgs tempVar(v, tmpvec);
        if ( WarnHandler != nullptr ) {
            WarnHandler->Invoke(this, &tempVar);
        }
    }
    
    void EverythingRunnerEngine::StartingAllTasks()
    {
        if ( StartingAllTasksEngineHandler != nullptr ) {
            StartingAllTasksEngineHandler->Invoke(this, EventArgs::Empty);
        }
    }
    
    void EverythingRunnerEngine::FinishedAllTasks(const std::string &rootOutputDir)
    {
        std::vector<std::string> tmpvec;
        StringEventArgs tempVar(rootOutputDir, tmpvec);
        if ( FinishedAllTasksEngineHandler != nullptr ) {
            FinishedAllTasksEngineHandler->Invoke(this, &tempVar);
        }
    }
    
    void EverythingRunnerEngine::NewSpectras(std::vector<std::string> &newSpectra)
    {
        StringListEventArgs tempVar(newSpectra);
        if ( NewSpectrasHandler != nullptr ) {
            NewSpectrasHandler->Invoke(this, &tempVar);
        }
    }
    
    void EverythingRunnerEngine::NewFileSpecificToml(std::vector<std::string> &newFileSpecificTomls)
    {
        StringListEventArgs tempVar(newFileSpecificTomls);
        if ( NewFileSpecificTomlHandler != nullptr ) {
            NewFileSpecificTomlHandler->Invoke(this, &tempVar);
        }
    }
    
    void EverythingRunnerEngine::NewDBs(std::vector<DbForTask*> &newDatabases)
    {
        XmlForTaskListEventArgs tempVar(newDatabases);
        if ( NewDbsHandler != nullptr ) {
            NewDbsHandler->Invoke(this, &tempVar);
        }
    }
}
