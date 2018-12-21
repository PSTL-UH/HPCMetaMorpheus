#include "Program.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineFinishedEventArgs.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/EventArgs/SingleTaskEventArgs.h"
#include "../EngineLayer/EventArgs/SingleFileEventArgs.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
#include "../TaskLayer/XLSearchTask/TaskLayer.XLSearchTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"

#include <getopt.h>

using namespace EngineLayer;
//using namespace Fclp; Has been replaced with getopt_long
//using namespace Nett; 
//using namespace Proteomics;
using namespace TaskLayer;


int main( int argc, char *argv[] )
{

    Program::Main ( argc, argv);
    return 0;
}


namespace MetaMorpheusCommandLine
{

    static void print_usage ( void )
    {
        std::wcout << "Usage: MetaMorpheusC++ -t <Task> -m <MetaTask> -o <OutputFolder> "
            "-s <Spectra> -d <Databases>" << std::endl;
        
        return;
    }
    
    static void print_config ( std::vector<std::string> taskname,
                               std::vector<std::string> metataskname,
                               std::vector<std::string> outputfolder,
                               std::vector<std::string> spectra,
                               std::vector<std::string> dbases )
        
    {
        std::wcout << "Number of tasks: " << taskname.size() << std::endl;
        for ( int i=0; i< taskname.size(); i ++ ) {
            std::wcout << "  Taskname " << i << " : " << taskname.at(i).c_str()
                       << std::endl;
        }
        
        std::wcout << "Number of Metatasks: " << metataskname.size() << std::endl;
        for ( int i=0; i< metataskname.size(); i ++ ) {
            std::wcout << "  Metataskname " << i << " : " << taskname.at(i).c_str()
                       << std::endl;
        }
        
        if ( outputfolder.size() == 0 ) {
            std::wcout << "Outputfolder: Not set " << std::endl;    
        }
        else {
            std::wcout << "Outputfolder: " << outputfolder.at(0).c_str() << std::endl;    
        }
        
        if ( spectra.size() == 0 ) {
            std::wcout << "Spectra: Not set" << std::endl;    
        }
        else {
            std::wcout << "Spectra: " << spectra.at(0).c_str() << std::endl;    
        }
        std::wcout << "Number of Databases: " << dbases.size() << std::endl;
        for ( int i=0; i< dbases.size(); i ++ ) {
            std::wcout << "  Databasename " << i << " : " << dbases.at(i).c_str()
                       << std::endl;
        }
        
        return;
    }
    
    
    bool Program::InProgress = false;
    System::CodeDom::Compiler::IndentedTextWriter *Program::MyWriter = new System::CodeDom::Compiler::IndentedTextWriter(Console::Out, L"\t");

    void Program::Main( int argc, char *argv[] )
    {
        std::wcout << L"Welcome to MetaMorpheus" << std::endl;
        std::wcout << GlobalVariables::getMetaMorpheusVersion() << std::endl;
        
        int opt=0;
        static struct option long_options[] = {
            {"tasks", required_argument, 0, 't'},
            {'outputFolder', required_argument, 0, 'o'},
            {'meta-task', required_argument, 0, 'm'},
            {'spectra', required_argument, 0, 's'},
            {'databases', required_argument, 0, 'd'},
            {0, 0, 0, 0}
        };
        
        int long_index=0;
        std::vector<std::string> taskname;
        std::vector<std::string> metataskname;
        std::vector<std::string> outputfolder;
        std::vector<std::string> spectra;
        std::vector<std::string> dbases;
        while ((opt = getopt_long (argc, argv, "t:o:m:s:d:",
                               long_options, &long_index)) != -1 ) {
            switch (opt) {
                case 't':
                    taskname.push_back (std::string(optarg));
                    break;
                case 'o':
                    outputfolder.push_back (std::string(optarg));
                    break;
                case 'm':
                    metataskname.push_back (std::string(optarg));
                    break;
                case 's':
                    spectra.push_back (std::string(optarg));
                    break;
                case 'd':
                    dbases.push_back (std::string(optarg));
                    break;
                default:
                    print_usage();
                    exit ( -1);
            }
        }

        print_config(taskname, metataskname, outputfolder, spectra, dbases);

        if (metataskname.size() != 0 || taskname.size() != 0) {
            if (!result->HasErrors)  {
                MetaMorpheusEngine::WarnHandler->addListener(L"WarnHandler", [&] (std::any sender, StringEventArgs* e) {WarnHandler(sender, e);});
                MetaMorpheusEngine::OutProgressHandler->addListener(L"MyEngine_outProgressHandler", [&] (std::any sender, ProgressEventArgs* e) {MyEngine_outProgressHandler(sender, e);});
                MetaMorpheusEngine::StartingSingleEngineHander->addListener(L"MyEngine_startingSingleEngineHander", [&] (std::any sender, SingleEngineEventArgs* e) {MyEngine_startingSingleEngineHander(sender, e);});
                MetaMorpheusEngine::FinishedSingleEngineHandler->addListener(L"MyEngine_finishedSingleEngineHandler", [&] (std::any sender, SingleEngineFinishedEventArgs* e) {MyEngine_finishedSingleEngineHandler(sender, e);});
                
                MetaMorpheusTask::WarnHandler->addListener(L"WarnHandler", [&] (std::any sender, StringEventArgs* e) {WarnHandler(sender, e);});
                MetaMorpheusTask::LogHandler->addListener(L"LogHandler", [&] (std::any sender, StringEventArgs* e) {LogHandler(sender, e);});
                MetaMorpheusTask::StartingSingleTaskHander->addListener(L"MyTaskEngine_startingSingleTaskHander", [&] (std::any sender, SingleTaskEventArgs* e) {MyTaskEngine_startingSingleTaskHander(sender, e);});
                MetaMorpheusTask::FinishedSingleTaskHandler->addListener(L"MyTaskEngine_finishedSingleTaskHandler", [&] (std::any sender, SingleTaskEventArgs* e) {MyTaskEngine_finishedSingleTaskHandler(sender, e);});
                MetaMorpheusTask::FinishedWritingFileHandler->addListener(L"MyTaskEngine_finishedWritingFileHandler", [&] (std::any sender, SingleFileEventArgs* e) {MyTaskEngine_finishedWritingFileHandler(sender, e);});
                
                for (auto db =begin(dbases); db != end(dbases); ++db ) {
                    if (Path::GetExtension(*db) != L".fasta")   {
                        GlobalVariables::AddMods(UsefulProteomicsDatabases::ProteinDbLoader::GetPtmListFromProteinXml(*db).OfType<Modification*>(), true);

                        // print any error messages reading the mods to the console
                        for (auto error : GlobalVariables::ErrorsReadingMods) {
                            std::wcout << error << std::endl;
                        }
                        GlobalVariables::ErrorsReadingMods.clear();
                    }
                }
                
                std::vector<(std::wstring, MetaMorpheusTask)*> taskList;

                for (auto i= begin(taskname); i != end (taskname) ; ++i) {
                    std::vector<string> filePath = *i;
                    auto uhum = Toml::ReadFile(filePath, MetaMorpheusTask::tomlConfig);

                    if (uhum->Get<std::wstring>(L"TaskType") == L"Search") {
                        auto ye1 = Toml::ReadFile<SearchTask*>(filePath, MetaMorpheusTask::tomlConfig);
                        taskList.push_back((L"Task" + std::to_wstring(i + 1) + L"SearchTask", ye1));
                    }
                    else if (uhum->Get<std::wstring>(L"TaskType") == L"Calibrate") {
                        auto ye2 = Toml::ReadFile<CalibrationTask*>(filePath, MetaMorpheusTask::tomlConfig);
                        taskList.push_back((L"Task" + std::to_wstring(i + 1) + L"CalibrationTask", ye2));
                    }
                    else if (uhum->Get<std::wstring>(L"TaskType") == L"Gptmd")  {
                        auto ye3 = Toml::ReadFile<GptmdTask*>(filePath, MetaMorpheusTask::tomlConfig);
                        taskList.push_back((L"Task" + std::to_wstring(i + 1) + L"GptmdTask", ye3));                        
                    }
                    else if (uhum->Get<std::wstring>(L"TaskType") == L"XLSearch")  {
                        auto ye4 = Toml::ReadFile<XLSearchTask*>(filePath, MetaMorpheusTask::tomlConfig);
                        taskList.push_back((L"Task" + std::to_wstring(i + 1) + L"XLSearchTask", ye4));
                    }
                    else {
                        std::wcout << uhum->Get<std::wstring>(L"TaskType") << L" is not a known task type! Skipping." << std::endl;
                    }
                }

                for (auto i = begin(metataskname); i != end(metataskname); ++i) {
                    auto filePath = *i;
                    auto uhum = Toml::ReadFile(filePath, MetaMorpheusTask::tomlConfig);

                    if (uhum->Get<std::wstring>(L"TaskType") == L"Search") {
                        std::wcout << L"Search tasks are individual tasks. Please use -t for task instead of -m. Skipping." << std::endl;
                    }
                    else if (uhum->Get<std::wstring>(L"TaskType") == L"Calibrate") {
                        std::wcout << L"Calibrate tasks are individual tasks. Please use -t for task instead of -m. Skipping." << std::endl;
                    }
                    else if (uhum->Get<std::wstring>(L"TaskType") == L"Gptmd")  {
                        std::wcout << L"Gptmd tasks are individual tasks. Please use -t for task instead of -m. Skipping." << std::endl;
                    }
                    else if (uhum->Get<std::wstring>(L"TaskType") == L"XLSearch")  {
                        std::wcout << L"XLSearch tasks are individual tasks. Please use -t for task instead of -m. Skipping." << std::endl;
                    }
                    else {
                        std::wcout << uhum->Get<std::wstring>(L"TaskType") << L" is not a known task type! Skipping." << std::endl;
                    }
                }

                std::vector<std::wstring> startingRawFilenameList = p->Object.Spectra->Select([&] (std::any b)
                                                                                              {
                                                                                                  FileSystem::getFullPath(b);
                                                                                              }).ToList();
                std::vector<DbForTask*> startingXmlDbFilenameList = p->Object.Databases->Select([&] (std::any b)
                                                                                                {
                                                                                                    new DbForTask(FileSystem::getFullPath(b), IsContaminant(b));
                                                                                                }).ToList();

                std::wstring outputFolder = p->Object.OutputFolder;
                if (outputFolder == L"")
                {
                    auto pathOfFirstSpectraFile = FileSystem::getDirectoryName(startingRawFilenameList.front());
                    outputFolder = FileSystem::combine(pathOfFirstSpectraFile, LR"($DATETIME)");
                }
                
                EverythingRunnerEngine *a = new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList, outputFolder);
                
                try
                {
                    a->Run();
                }
                catch (const std::runtime_error &e)
                {
                    while (e.InnerException != nullptr)
                    {
                        e = e.InnerException;
                    }
                    auto message = L"Run failed, Exception: " + e.what();
                    std::wcout << message << std::endl;
                }
                
                delete a;
            }
            else
            {
                std::wcout << L"Error Text:" << result->ErrorText << std::endl;
            }
        }
        else
        {
            std::wcout << L"Error Text: No toml file was specified. Use -t for tasks or -m for meta-tasks." << std::endl;
        }
        
    }

    void Program::WriteMultiLineIndented(const std::wstring &toWrite)
    {
        std::vector<std::wstring> tokens = Regex::Split(toWrite, LR"(\r?\n|\r)");
        for (auto str : tokens)
        {
            MyWriter->WriteLine(str);
        }
    }
    
    bool Program::IsContaminant(const std::wstring &b)
    {
        if (StringHelper::toUpper(b).find(StringHelper::toUpper(L"contaminant")) != std::wstring::npos || StringHelper::toUpper(b).find(L"CRAP") != std::wstring::npos)
        {
            return true;
        }
        return false;
    }
    
    void Program::MyTaskEngine_startingSingleTaskHander(std::any sender, SingleTaskEventArgs *e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented(L"Starting task: " + e->getDisplayName());
        MyWriter->Indent = MyWriter->Indent + 1;
    }
    
    void Program::MyTaskEngine_finishedWritingFileHandler(std::any sender, SingleFileEventArgs *e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented(L"Finished writing file: " + e->getWrittenFile());
    }
    
    void Program::MyTaskEngine_finishedSingleTaskHandler(std::any sender, SingleTaskEventArgs *e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        MyWriter->Indent = MyWriter->Indent - 1;
        WriteMultiLineIndented(L"Finished task: " + e->getDisplayName());
    }
    
    void Program::MyEngine_startingSingleEngineHander(std::any sender, SingleEngineEventArgs *e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented(L"Starting engine: " + e->getMyEngine()->GetType()->Name + L" " + e->getMyEngine()->GetId());
        MyWriter->Indent = MyWriter->Indent + 1;
    }
    
    void Program::MyEngine_finishedSingleEngineHandler(std::any sender, SingleEngineFinishedEventArgs *e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented(L"Engine results: " + e);
        MyWriter->Indent = MyWriter->Indent - 1;
        WriteMultiLineIndented(L"Finished engine: " + e->MyResults->getMyEngine()->GetType()->Name + L" " + e->MyResults->getMyEngine()->GetId());
    }
    
    void Program::MyEngine_outProgressHandler(std::any sender, ProgressEventArgs *e)
    {
        MyWriter->Write(std::to_wstring(e->NewProgress) + L" ");
        InProgress = true;
    }
    
    void Program::WarnHandler(std::any sender, StringEventArgs *e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented(L"WARN: " + e->getS());
    }
    
    void Program::LogHandler(std::any sender, StringEventArgs *e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented(L"Log: " + e->getS());
    }
    
    std::vector<std::wstring> Program::ApplicationArguments::getTasks() const
    {
        return privateTasks;
    }
    
    void Program::ApplicationArguments::setTasks(const std::vector<std::wstring> &value)
    {
        privateTasks = value;
    }
    
    std::vector<std::wstring> Program::ApplicationArguments::getDatabases() const
    {
        return privateDatabases;
    }
    
    void Program::ApplicationArguments::setDatabases(const std::vector<std::wstring> &value)
    {
        privateDatabases = value;
    }
    
    std::vector<std::wstring> Program::ApplicationArguments::getSpectra() const
    {
        return privateSpectra;
    }
    
    void Program::ApplicationArguments::setSpectra(const std::vector<std::wstring> &value)
    {
        privateSpectra = value;
    }
    
    std::vector<std::wstring> Program::ApplicationArguments::getMetaTasks() const
    {
        return privateMetaTasks;
    }
    
    void Program::ApplicationArguments::setMetaTasks(const std::vector<std::wstring> &value)
    {
        privateMetaTasks = value;
    }
    
    std::wstring Program::ApplicationArguments::getOutputFolder() const
    {
        return privateOutputFolder;
    }
    
    void Program::ApplicationArguments::setOutputFolder(const std::wstring &value)
    {
        privateOutputFolder = value;
    }
}
