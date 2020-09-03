#include "EverythingRunnerEngine.h"
#include "MetaMorpheusTask.h"
#include "DbForTask.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/StringListEventArgs.h"
#include "EventHandler.h"
#include "EventArgs/XmlForTaskListEventArgs.h"
#include "CalibrationTask/CalibrationTask.h"

#include <filesystem>
#include <ctime>
#include <stdio.h>
#include <sys/time.h>


static double timediff (struct timeval t1, struct timeval t2)
{
    double elapsedtime;
    elapsedtime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedtime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

    return elapsedtime/1000;                            //ms to sec
}


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
        
    }
    
    void EverythingRunnerEngine::Run()
    {
        time_t timer;
        time(&timer);
        struct tm *tmi_start = localtime(&timer);

        struct timeval t1, t1e;
        gettimeofday (&t1, NULL);
        
        if (CurrentRawDataFilenameList.empty())
        {
            MetaMorpheusTask::Warn("No spectra files selected");
            return;
        }
        
        std::string startTimeForAllFilenames = std::to_string(tmi_start->tm_year) + "-" +
            std::to_string(tmi_start->tm_mon) + "-" +  std::to_string(tmi_start->tm_mday) + "-" +
            std::to_string(tmi_start->tm_hour) + "-" + std::to_string(tmi_start->tm_min) + "-" +
            std::to_string(tmi_start->tm_sec);
        
        OutputFolder = StringHelper::replace(OutputFolder, "$DATETIME", startTimeForAllFilenames);
        StringBuilder *allResultsText = new StringBuilder();
        
        for (int i = 0; i < (int)TaskList.size(); i++)
        {
            if (CurrentRawDataFilenameList.empty())
            {
                MetaMorpheusTask::Warn("Cannot proceed. No spectra files selected.");
                delete allResultsText;
                return;
            }
            if (CurrentXmlDbFilenameList.empty())
            {
                MetaMorpheusTask::Warn("Cannot proceed. No protein database files selected.");
                delete allResultsText;
                return;
            }
            auto ok = TaskList[i];
            
            auto outputFolderForThisTask = OutputFolder + "/" + std::get<0>(ok);
            if (!std::filesystem::exists(outputFolderForThisTask))
            {
                std::filesystem::create_directory(outputFolderForThisTask);
            }
            
            // Actual task running code
            auto myTaskResults = std::get<1>(ok)->RunTask(outputFolderForThisTask,
                                                          CurrentXmlDbFilenameList,
                                                          CurrentRawDataFilenameList,
                                                          std::get<0>(ok));
            
            if (!myTaskResults->NewDatabases.empty() )
            {
                CurrentXmlDbFilenameList = myTaskResults->NewDatabases;
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
                    std::vector<std::string> successfullyCalibFiles;
                    for (auto p: myTaskResults->NewSpectra ) {
                        std::string fname = p.substr(0, p.find_last_of("."));
                        successfullyCalibFiles.push_back( StringHelper::replace(fname, CalibrationTask::CalibSuffix, ""));
                    }

                    std::vector<std::string> origFiles;
                    for ( auto p: CurrentRawDataFilenameList ) {
                        origFiles.push_back(p.substr(0,p.find_last_of(".")));
                    }

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
            }

            std::string sst= "\r\n\r\n\r\n\r\n";
            allResultsText->appendLine( sst + myTaskResults->ToString());
        }

        gettimeofday (&t1e, NULL);
        double time = timediff (t1, t1e);
        int time_hour = (int) (time/3600);
        int time_min = (int) ((time - (time_hour*3600))/60);
        int time_sec = time - ((time_min * 60) + (time_hour *3600));
        char timestr[64];
        sprintf(timestr, "%02d:%02d:%02d", time_hour, time_min, time_sec);
        
        auto resultsFileName = OutputFolder + "/allResults.txt";
        std::ofstream file(resultsFileName);
        file << "MetaMorpheus: version " << GlobalVariables::getMetaMorpheusVersion() << std::endl;
        file <<"Total time: " << timestr << " " << time << std::endl;
        file << allResultsText->toString() << std::endl;
        
        //FinishedAllTasks(OutputFolder);
        
        delete allResultsText;
    }
    
}
