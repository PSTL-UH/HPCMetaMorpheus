#include "Program.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/EventArgs/SingleFileEventArgs.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineFinishedEventArgs.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/EventArgs/SingleTaskEventArgs.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
//#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
#include "../TaskLayer/XLSearchTask/TaskLayer.XLSearchTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"

#include <getopt.h>
#include <experimental/filesystem>

using namespace EngineLayer;
using namespace TaskLayer;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;


int main( int argc, char *argv[] )
{

    MetaMorpheusCommandLine::Program::Main ( argc, argv);
    return 0;
}


namespace MetaMorpheusCommandLine
{

    static void print_usage ( void )
    {
        std::cout << "Usage: MetaMorpheusC++ -t <Task> -m <MetaTask> -o <OutputFolder> "
            "-s <Spectra> -d <Databases>" << std::endl;
        
        return;
    }
    
    static void print_config ( std::vector<std::string> taskname,
                               std::vector<std::string> metataskname,
                               std::vector<std::string> outputfolder,
                               std::vector<std::string> spectra,
                               std::vector<std::string> dbases )
        
    {
        std::cout << "Number of tasks: " << taskname.size() << std::endl;
        for ( int i=0; i< (int)taskname.size(); i ++ ) {
            std::cout << "  Taskname " << i << " : " << taskname.at(i).c_str()
                       << std::endl;
        }
        
        std::cout << "Number of Metatasks: " << metataskname.size() << std::endl;
        for ( int i=0; i< (int)metataskname.size(); i ++ ) {
            std::cout << "  Metataskname " << i << " : " << taskname.at(i).c_str()
                       << std::endl;
        }
        
        if ( outputfolder.size() == 0 ) {
            std::cout << "Outputfolder: Not set " << std::endl;    
        }
        else {
            std::cout << "Outputfolder: " << outputfolder.at(0).c_str() << std::endl;    
        }
        
        if ( spectra.size() == 0 ) {
            std::cout << "Spectra: Not set" << std::endl;    
        }
        else {
            std::cout << "Spectra: " << spectra.at(0).c_str() << std::endl;    
        }
        std::cout << "Number of Databases: " << dbases.size() << std::endl;
        for ( int i=0; i< (int)dbases.size(); i ++ ) {
            std::cout << "  Databasename " << i << " : " << dbases.at(i).c_str()
                       << std::endl;
        }
        
        return;
    }
    
    
    bool Program::InProgress = false;

    void Program::Main( int argc, char *argv[] )
    {
        std::cout << "Welcome to MetaMorpheus" << std::endl;
        std::cout << GlobalVariables::getMetaMorpheusVersion() << std::endl;

        //IndentedTextWriter *Program::MyWriter = new IndentedTextWriter(Console::Out, "\t");
        MyWriter = new IndentedTextWriter();

        int opt=0;
        static struct option long_options[] = {
            {"tasks", required_argument, 0, 't'},
            {"outputFolder", required_argument, 0, 'o'},
            {"meta-task", required_argument, 0, 'm'},
            {"spectra", required_argument, 0, 's'},
            {"databases", required_argument, 0, 'd'},
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
                MetaMorpheusEngine::WarnHandler->addListener("WarnHandler", WarnHandler );
                MetaMorpheusEngine::OutProgressHandler->addListener("MyEngine_outProgressHandler",
                                                                    MyEngine_outProgressHandler );
                MetaMorpheusEngine::StartingSingleEngineHandler->addListener("MyEngine_startingSingleEngineHandler",
                                                                            MyEngine_startingSingleEngineHandler);
                MetaMorpheusEngine::FinishedSingleEngineHandler->addListener("MyEngine_finishedSingleEngineHandler",
                                                                             MyEngine_finishedSingleEngineHandler);
                
                MetaMorpheusTask::LogHandler->addListener("LogHandler", LogHandler);
                MetaMorpheusTask::StartingSingleTaskHandler->addListener("MyTaskEngine_startingSingleTaskHandler", 
                                                                        MyTaskEngine_startingSingleTaskHandler);
                MetaMorpheusTask::FinishedSingleTaskHandler->addListener("MyTaskEngine_finishedSingleTaskHandler", 
                                                                         MyTaskEngine_finishedSingleTaskHandler);
                MetaMorpheusTask::FinishedWritingFileHandler->addListener("MyTaskEngine_finishedWritingFileHandler",
                                                                          MyTaskEngine_finishedWritingFileHandler);
                for (auto db =dbases.begin(); db != dbases.end(); ++db ) {
                    if (Path::GetExtension(*db) != ".fasta")   {
                        GlobalVariables::AddMods(UsefulProteomicsDatabases::ProteinDbLoader::GetPtmListFromProteinXml(*db).OfType<Modification*>(), true);

                        // print any error messages reading the mods to the console
                        for (auto error : GlobalVariables::ErrorsReadingMods) {
                            std::cout << error << std::endl;
                        }
                        GlobalVariables::ErrorsReadingMods.clear();
                    }
                }
                
                std::vector<std::tuple<std::string, MetaMorpheusTask*>> taskList;

                for (auto i= taskname.begin(); i != taskname.end() ; ++i) {
                    std::vector<std::string> filePath = *i;
                    auto uhum = Toml::ReadFile(filePath, MetaMorpheusTask::tomlConfig);

                    if (uhum->Get<std::string>("TaskType") == "Search") {
                        auto ye1 = Toml::ReadFile<SearchTask*>(filePath, MetaMorpheusTask::tomlConfig);
                        taskList.push_back(std::make_tuple("Task" + std::to_string(i + 1) + "SearchTask", ye1));
                    }
                    else if (uhum->Get<std::string>("TaskType") == "Calibrate") {
                        auto ye2 = Toml::ReadFile<CalibrationTask*>(filePath, MetaMorpheusTask::tomlConfig);
                        taskList.push_back(std::make_tuple("Task" + std::to_string(i + 1) + "CalibrationTask", ye2));
                    }
                    else if (uhum->Get<std::string>("TaskType") == "Gptmd")  {
                        auto ye3 = Toml::ReadFile<GptmdTask*>(filePath, MetaMorpheusTask::tomlConfig);
                        taskList.push_back(std::make_tuple("Task" + std::to_string(i + 1) + "GptmdTask", ye3));                        
                    }
                    else if (uhum->Get<std::string>("TaskType") == "XLSearch")  {
                        auto ye4 = Toml::ReadFile<XLSearchTask*>(filePath, MetaMorpheusTask::tomlConfig);
                        taskList.push_back(std::make_tuple("Task" + std::to_string(i + 1) + "XLSearchTask", ye4));
                    }
                    else {
                        std::cout << uhum->Get<std::string>("TaskType") << " is not a known task type! Skipping." << std::endl;
                    }
                }

                for (auto i = metataskname.begin(); i != metataskname.end(); ++i) {
                    auto filePath = *i;
                    auto uhum = Toml::ReadFile(filePath, MetaMorpheusTask::tomlConfig);

                    if (uhum->Get<std::string>("TaskType") == "Search") {
                        std::cout << "Search tasks are individual tasks. Please use -t for task instead of -m. Skipping." <<
                            std::endl;
                    }
                    else if (uhum->Get<std::string>("TaskType") == "Calibrate") {
                        std::cout << "Calibrate tasks are individual tasks. Please use -t for task instead of -m. Skipping." <<
                            std::endl;
                    }
                    else if (uhum->Get<std::string>("TaskType") == "Gptmd")  {
                        std::cout << "Gptmd tasks are individual tasks. Please use -t for task instead of -m. Skipping." <<
                            std::endl;
                    }
                    else if (uhum->Get<std::string>("TaskType") == "XLSearch")  {
                        std::cout << "XLSearch tasks are individual tasks. Please use -t for task instead of -m. Skipping." <<
                            std::endl;
                    }
                    else {
                        std::cout << uhum->Get<std::string>("TaskType") << " is not a known task type! Skipping." <<
                            std::endl;
                    }
                }

#ifdef ORIG
                std::vector<std::string> startingRawFilenameList = p->Object.Spectra->Select([&] (std::any b)  {
                        FileSystem::getFullPath(b);
                    }).ToList();
#endif
                std::vector<std::string> startingRawFilenameList;
                for ( auto b: p->Object().Spectra() ) {
                    startingRawFilenameList.push_back(std::experimental::filesystem::absolute(b));
                }

#ifdef ORIG
                std::vector<DbForTask*> startingXmlDbFilenameList = p->Object.Databases->Select([&] (std::any b) {
                        new DbForTask(FileSystem::getFullPath(b), IsContaminant(b));
                    }).ToList();
#endif
                std::vector<DbForTask*> startingXmlDbFilenameList;
                for ( auto b: p->Object().Databases()) {
                    auto newdb = new DbForTask(std::experimental::filesystem::absolute(b), IsContaminant(b));
                    startingXmlDbFilenameList.push_back(newdb);       
                }

                std::string outputFolder = p->Object.OutputFolder;
                if (outputFolder == "")
                {
                    auto pathOfFirstSpectraFile = FileSystem::getDirectoryName(startingRawFilenameList.front());
                    outputFolder = pathOfFirstSpectraFile+ "/($DATETIME)";
                }
                
                EverythingRunnerEngine *a = new EverythingRunnerEngine(taskList, startingRawFilenameList,
                                                                       startingXmlDbFilenameList, outputFolder);
                
                try
                {
                    a->Run();
                }
                catch (const std::runtime_error &e)
                {
                    auto message = "Run failed, Exception: " + e.what();
                    std::cout << message << std::endl;
                }
                
                delete a;
            }
            else
            {
                std::cout << "Error Text:" << result->ErrorText << std::endl;
            }
        }
        else
        {
            std::cout << "Error Text: No toml file was specified. Use -t for tasks or -m for meta-tasks." <<
                std::endl;
        }
        
    }

    void Program::WriteMultiLineIndented(const std::string &toWrite)
    {
        std::vector<std::string> tokens = Regex::Split(toWrite, R"(\r?\n|\r)");
        for (auto str : tokens)
        {
            MyWriter->WriteLine(str);
        }
    }
    
    bool Program::IsContaminant(const std::string &b)
    {
        if (StringHelper::toUpper(b).find(StringHelper::toUpper("contaminant")) != std::string::npos ||
            StringHelper::toUpper(b).find("CRAP") != std::string::npos)
        {
            return true;
        }
        return false;
    }
    
    void Program::MyTaskEngine_startingSingleTaskHandler( SingleTaskEventArgs e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented("Starting task: " + e.getDisplayName());
        MyWriter->Indent = MyWriter->Indent + 1;
    }
    
    void Program::MyTaskEngine_finishedWritingFileHandler( SingleFileEventArgs e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented("Finished writing file: " + e.getWrittenFile());
    }
    
    void Program::MyTaskEngine_finishedSingleTaskHandler( SingleTaskEventArgs e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        MyWriter->Indent = MyWriter->Indent - 1;
        WriteMultiLineIndented("Finished task: " + e.getDisplayName());
    }
    
    void Program::MyEngine_startingSingleEngineHandler( SingleEngineEventArgs e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented("Starting engine: " + e.getMyEngine()->GetType()->Name + " " +
                               e.getMyEngine()->GetId());
        MyWriter->Indent = MyWriter->Indent + 1;
    }
    
    void Program::MyEngine_finishedSingleEngineHandler( SingleEngineFinishedEventArgs e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented("Engine results: " + e);
        MyWriter->Indent = MyWriter->Indent - 1;
        WriteMultiLineIndented("Finished engine: " + e.MyResults->getMyEngine()->GetType()->Name + " " +
                               e.MyResults->getMyEngine()->GetId());
    }
    
    void Program::MyEngine_outProgressHandler( ProgressEventArgs e)
    {
        MyWriter->Write(std::to_string(e.NewProgress) + " ");
        InProgress = true;
    }

    void Program::WarnHandler(StringEventArgs e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented("WARN: " + e.getS());
    }
    
    void Program::LogHandler( StringEventArgs e)
    {
        if (InProgress)
        {
            MyWriter->WriteLine();
        }
        InProgress = false;
        WriteMultiLineIndented("Log: " + e.getS());
    }
    
    std::vector<std::string> Program::ApplicationArguments::getTasks() const
    {
        return privateTasks;
    }
    
    void Program::ApplicationArguments::setTasks(const std::vector<std::string> &value)
    {
        privateTasks = value;
    }
    
    std::vector<std::string> Program::ApplicationArguments::getDatabases() const
    {
        return privateDatabases;
    }
    
    void Program::ApplicationArguments::setDatabases(const std::vector<std::string> &value)
    {
        privateDatabases = value;
    }
    
    std::vector<std::string> Program::ApplicationArguments::getSpectra() const
    {
        return privateSpectra;
    }
    
    void Program::ApplicationArguments::setSpectra(const std::vector<std::string> &value)
    {
        privateSpectra = value;
    }
    
    std::vector<std::string> Program::ApplicationArguments::getMetaTasks() const
    {
        return privateMetaTasks;
    }
    
    void Program::ApplicationArguments::setMetaTasks(const std::vector<std::string> &value)
    {
        privateMetaTasks = value;
    }
    
    std::string Program::ApplicationArguments::getOutputFolder() const
    {
        return privateOutputFolder;
    }
    
    void Program::ApplicationArguments::setOutputFolder(const std::string &value)
    {
        privateOutputFolder = value;
    }
}
