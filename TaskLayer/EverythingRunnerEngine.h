#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_filesystem.h"
#include "EventHandler.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/StringListEventArgs.h"
#include "EventArgs/XmlForTaskListEventArgs.h"

#include "MetaMorpheusTask.h"
#include "DbForTask.h"

using namespace EngineLayer;

namespace TaskLayer
{
    class EverythingRunnerEngine
    {
    private:
        const std::vector<std::tuple<std::string, MetaMorpheusTask*>> TaskList;
        std::string OutputFolder;
        std::vector<std::string> CurrentRawDataFilenameList;
        std::vector<DbForTask*> CurrentXmlDbFilenameList;
        
    public:
        EverythingRunnerEngine(std::vector<std::tuple<std::string, MetaMorpheusTask*>> &taskList,
                               std::vector<std::string> &startingRawFilenameList,
                               std::vector<DbForTask*> &startingXmlDbFilenameList,
                               const std::string &outputFolder);
        
        static EventHandler<StringEventArgs> *StartingAllTasksEngineHandler;
        static EventHandler<StringEventArgs> *FinishedAllTasksEngineHandler;
        static EventHandler<XmlForTaskListEventArgs> *NewDbsHandler;
        static EventHandler<StringListEventArgs> *NewSpectrasHandler;
        static EventHandler<StringListEventArgs> *NewFileSpecificTomlHandler;
        static EventHandler<StringEventArgs> *WarnHandler;
        static EventHandler<StringEventArgs> *FinishedWritingAllResultsFileHandler;
            
        void Run();
        
    private:
        void Warn(const std::string &v);
        
        void StartingAllTasks();
        
        void FinishedAllTasks(const std::string &rootOutputDir);
        
        void NewSpectras(std::vector<std::string> &newSpectra);
        
        void NewFileSpecificToml(std::vector<std::string> &newFileSpecificTomls);
        
        void NewDBs(std::vector<DbForTask*> &newDatabases);
    };
}
