#include "EverythingRunnerEngine.h"
#include "MetaMorpheusTask.h"
#include "DbForTask.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/StringListEventArgs.h"
#include "EventHandler.h"
#include "EventArgs/XmlForTaskListEventArgs.h"
#include "CalibrationTask/CalibrationTask.h"

#include <experimental/filesystem>
#include <ctime>

using namespace EngineLayer;

namespace TaskLayer
{

    EverythingRunnerEngine::EverythingRunnerEngine(std::vector<std::tuple<std::string, MetaMorpheusTask*>> &taskList,
                                                   std::vector<std::string> &startingRawFilenameList,
                                                   std::vector<DbForTask*> &startingXmlDbFilenameList,
                                                   const std::string &outputFolder) : TaskList(taskList)
    {
        OutputFolder = StringHelper::trim(outputFolder, "");
        
        CurrentRawDataFilenameList = startingRawFilenameList;
        CurrentXmlDbFilenameList = startingXmlDbFilenameList;
        
        StartingAllTasksEngineHandler = new EventHandler<StringEventArgs>();
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
        time_t timer;
        time(&timer);
        struct tm *tmi_start = localtime(&timer);
        clock_t begin = clock();
        
        if (!CurrentRawDataFilenameList.empty())
        {
            Warn("No spectra files selected");
            FinishedAllTasks("");
            
            return;
        }
        
        //auto startTimeForAllFilenames = DateTime::Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo::InvariantCulture);
        std::string startTimeForAllFilenames = std::to_string(tmi_start->tm_year) + "-" +
            std::to_string(tmi_start->tm_mon) + "-" +  std::to_string(tmi_start->tm_mday) + "-" +
            std::to_string(tmi_start->tm_hour) + "-" + std::to_string(tmi_start->tm_min) + "-" +
            std::to_string(tmi_start->tm_sec);
        
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
            
            auto outputFolderForThisTask = OutputFolder + std::get<0>(ok);
            
            if (!std::experimental::filesystem::exists(outputFolderForThisTask))
            {
                std::experimental::filesystem::create_directory(outputFolderForThisTask);
            }
            
            // Actual task running code
            auto myTaskResults = std::get<1>(ok)->RunTask(outputFolderForThisTask,
                                                          CurrentXmlDbFilenameList,
                                                          CurrentRawDataFilenameList,
                                                          std::get<0>(ok));
            
            if (!myTaskResults->NewDatabases.empty() )
            {
                CurrentXmlDbFilenameList = myTaskResults->NewDatabases;
                NewDBs(myTaskResults->NewDatabases);
            }
            if (!myTaskResults->NewSpectra.empty())
            {
                if (CurrentRawDataFilenameList.size() == myTaskResults->NewSpectra.size() )
                {
                    CurrentRawDataFilenameList = myTaskResults->NewSpectra;
                }
                else
                {
                    // at least one file was not successfully calibrated
#ifdef ORIG
                    auto successfullyCalibFiles = myTaskResults->NewSpectra->Select([&] (std::any p)  {
                            StringHelper::replace(Path::GetFileNameWithoutExtension(p), CalibrationTask::CalibSuffix, "");
                        }).ToList();
#endif
                    std::vector<std::string> successfullyCalibFiles;
                    for (auto p: myTaskResults->NewSpectra ) {
                        std::string fname = p.substr(0, p.find_last_of("."));
                        successfullyCalibFiles.push_back( StringHelper::replace(fname, CalibrationTask::CalibSuffix, ""));
                    }

#ifdef ORIG
                    auto origFiles = CurrentRawDataFilenameList.Select([&] (std::any p) {
                            Path::GetFileNameWithoutExtension(p);
                        }).ToList();
#endif
                    std::vector<std::string> origFiles;
                    for ( auto p: CurrentRawDataFilenameList ) {
                        origFiles.push_back(p.substr(0,p.find_last_of(".")));
                    }
#ifdef ORIG
                    auto unsuccessfullyCalibFiles = origFiles.Except(successfullyCalibFiles).ToList();
#endif
                    std::vector<std::string> unsuccessfullyCalibFiles;
                    for ( auto p: origFiles ) {
                        bool found = false;
                        for (auto q: successfullyCalibFiles ) {
                            if ( q == p ) {
                                found = true;
                                break;
                            }
                        }
                        if ( !found ) {
                            unsuccessfullyCalibFiles.push_back(p);
                        }
                    }
#ifdef ORIG
                    auto unsuccessfullyCalibFilePaths = CurrentRawDataFilenameList.Where([&] (std::any p)   {
                            std::find(unsuccessfullyCalibFiles.begin(), unsuccessfullyCalibFiles.end(),
                                      Path::GetFileNameWithoutExtension(p)) != unsuccessfullyCalibFiles.end();
                        }).ToList();
#endif
                    std::vector<std::string> unsuccessfullyCalibFilePaths;
                    for ( auto p: CurrentRawDataFilenameList ) {
                        std::string p_dash = p.substr(0, p.find_last_of("."));
                        if ( std::find(unsuccessfullyCalibFiles.begin(), unsuccessfullyCalibFiles.end(),
                                       p_dash) != unsuccessfullyCalibFiles.end() ){
                            unsuccessfullyCalibFilePaths.push_back(p);
                        }
                    }
                    
                    CurrentRawDataFilenameList = myTaskResults->NewSpectra;

                    CurrentRawDataFilenameList.insert(CurrentRawDataFilenameList.end(), unsuccessfullyCalibFilePaths.begin(),
                                                      unsuccessfullyCalibFilePaths.end());
                } 
                
                NewSpectras(myTaskResults->NewSpectra);
            }
            if (!myTaskResults->NewFileSpecificTomls.empty() )
            {
                NewFileSpecificToml(myTaskResults->NewFileSpecificTomls);
            }

            std::string sst= "\r\n\r\n\r\n\r\n";
            allResultsText->appendLine( sst + myTaskResults->ToString());
        }

        clock_t end = clock();
        
        auto resultsFileName = OutputFolder + "/allResults.txt";
        std::ofstream file(resultsFileName);
        file << "MetaMorpheus: version " << GlobalVariables::getMetaMorpheusVersion() << std::endl;
        file <<"Total time: " << (end - begin)/CLOCKS_PER_SEC << std::endl;
        file << allResultsText->toString() << std::endl;
        
        std::vector<std::string> tmpvec;
        StringEventArgs tempVar(resultsFileName, tmpvec);
        if ( FinishedWritingAllResultsFileHandler != nullptr ) {
            FinishedWritingAllResultsFileHandler->Invoke(tempVar);
        }
        FinishedAllTasks(OutputFolder);
        
        delete allResultsText;
    }
    
    void EverythingRunnerEngine::Warn(const std::string &v)
    {
        std::vector<std::string> tmpvec;
        StringEventArgs tempVar(v, tmpvec);
        if ( WarnHandler != nullptr ) {
            WarnHandler->Invoke(tempVar);
        }
    }
    
    void EverythingRunnerEngine::StartingAllTasks()
    {
        std::vector<std::string> tmpvec;
        std::string v;
        StringEventArgs tempVar(v, tmpvec);
        if ( StartingAllTasksEngineHandler != nullptr ) {
            StartingAllTasksEngineHandler->Invoke(tempVar);
        }
    }
    
    void EverythingRunnerEngine::FinishedAllTasks(const std::string &rootOutputDir)
    {
        std::vector<std::string> tmpvec;
        StringEventArgs tempVar(rootOutputDir, tmpvec);
        if ( FinishedAllTasksEngineHandler != nullptr ) {
            FinishedAllTasksEngineHandler->Invoke(tempVar);
        }
    }
    
    void EverythingRunnerEngine::NewSpectras(std::vector<std::string> &newSpectra)
    {
        StringListEventArgs tempVar(newSpectra);
        if ( NewSpectrasHandler != nullptr ) {
            NewSpectrasHandler->Invoke(tempVar);
        }
    }
    
    void EverythingRunnerEngine::NewFileSpecificToml(std::vector<std::string> &newFileSpecificTomls)
    {
        StringListEventArgs tempVar(newFileSpecificTomls);
        if ( NewFileSpecificTomlHandler != nullptr ) {
            NewFileSpecificTomlHandler->Invoke(tempVar);
        }
    }
    
    void EverythingRunnerEngine::NewDBs(std::vector<DbForTask*> &newDatabases)
    {
        XmlForTaskListEventArgs tempVar(newDatabases);
        if ( NewDbsHandler != nullptr ) {
            NewDbsHandler->Invoke(tempVar);
        }
    }
}
