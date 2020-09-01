#include <vector>
#include <string>
#include <getopt.h>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <ctime>

#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
//#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
#include "../TaskLayer/XLSearchTask/XLSearchTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"

#include "stringhelper.h"

using namespace EngineLayer;
using namespace TaskLayer;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;

static void print_usage ( void )
{
    std::cout << "Usage: MetaMorpheusC++ -t <Task> -m <MetaTask> -o <OutputFolder> "
        "-s <Spectra> -d <Databases> -v <verbositLevel> "<< std::endl;
    
    return;
}

static void print_config ( std::vector<std::string> &taskname,
                           std::vector<std::string> &metataskname,
                           std::string &outputfolder,
                           std::vector<std::string> &spectra,
                           std::vector<std::string> &dbases )
    
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
        
        if ( outputfolder.length() == 0 ) {
            std::cout << "Outputfolder: Not set " << std::endl;    
        }
        else {
            std::cout << "Outputfolder: " << outputfolder << std::endl;    
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

static bool IsContaminant(const std::string &b)
{
    if (StringHelper::toUpper(b).find(StringHelper::toUpper("contaminant")) != std::string::npos ||
        StringHelper::toUpper(b).find("CRAP") != std::string::npos)
    {
        return true;
    }
    return false;
}

int main( int argc, char *argv[] )
{
    std::vector<std::string> Tasks;
    std::vector<std::string> Databases;
    std::vector<std::string> Spectra;
    std::vector<std::string> MetaTasks;
    std::string OutputFolder;
    
    std::cout << "Welcome to HPCMetaMorpheus" << std::endl;
    std::cout << GlobalVariables::getMetaMorpheusVersion() << std::endl;

    int verbosity=1; // Setting the default to provide Status updates
    int opt=0;
    static struct option long_options[] = {
        {"tasks", required_argument, 0, 't'},
        {"outputFolder", required_argument, 0, 'o'},
        {"meta-task", required_argument, 0, 'm'},
        {"spectra", required_argument, 0, 's'},
        {"databases", required_argument, 0, 'd'},
        {"verbosity", required_argument, 0, 'v'},
        {0, 0, 0, 0}
    };
    
    int long_index=0;
    while ((opt = getopt_long (argc, argv, "t:o:m:s:d:v:",
                               long_options, &long_index)) != -1 ) {
        switch (opt) {
            case 't':
                Tasks.push_back (std::string(optarg));
                break;
            case 'o':
                OutputFolder = std::string(optarg) ;
                break;
            case 'm':
                MetaTasks.push_back (std::string(optarg));
                break;
            case 's':
                Spectra.push_back (std::string(optarg));
                break;
            case 'd':
                Databases.push_back (std::string(optarg));
                break;
            case 'v':
                verbosity = std::stoi(optarg);
            default:
                print_usage();
                exit ( -1);
        }
    }
    
    print_config(Tasks, MetaTasks, OutputFolder, Spectra, Databases);
    
    if ( MetaTasks.size() != 0 || Tasks.size() != 0) {
        for (auto db = Databases.begin(); db != Databases.end(); ++db ) {
            if ( (*db).substr((*db).find_last_of(".")) != ".fasta")   {
                //GlobalVariables::AddMods(UsefulProteomicsDatabases::ProteinDbLoader::GetPtmListFromProteinXml(*db).OfType<Modification*>(), true);
                auto mods = UsefulProteomicsDatabases::ProteinDbLoader::GetPtmListFromProteinXml(*db);
                GlobalVariables::AddMods (mods, true );
            }
        }
        
        std::vector<std::tuple<std::string, MetaMorpheusTask*>> taskList;

        int i = 0;
        for (auto task : Tasks ) {
            std::string filePath = task;
            //read toml file
            Toml trw;
            toml::Value uhum = trw.tomlReadFile(filePath);
            
            if (uhum.find("TaskType")->as<std::string>() == "Search") {
                auto ye1 = new SearchTask(filePath);
                ye1->setVerbose(verbosity);
                taskList.push_back(std::make_tuple("Task" + std::to_string(i + 1) + "SearchTask", ye1));
            }
            else if (uhum.find("TaskType")->as<std::string>() == "Calibrate") {
                auto ye2 = new CalibrationTask(filePath);
                ye2->setVerbose(verbosity);
                taskList.push_back(std::make_tuple("Task" + std::to_string(i + 1) + "CalibrationTask", ye2));
            }
            else if (uhum.find("TaskType")->as<std::string>() == "Gptmd")  {
                //auto ye3 = new Gptmd(filePath);
                //ye3->setVerbose(verbosity);
                //taskList.push_back(std::make_tuple("Task" + std::to_string(i + 1) + "GptmdTask", ye3));
                std::cout << "Gptmd tasks are currently not support by HPCMetaMorpheus \n";
            }
            else if (uhum.find("TaskType")->as<std::string>() == "XLSearch")  {
                auto ye4 = new XLSearchTask(filePath);
                ye4->setVerbose(verbosity);
                taskList.push_back(std::make_tuple("Task" + std::to_string(i + 1) + "XLSearchTask", ye4));
            }
            else {
                std::cout << uhum.find("TaskType")->as<std::string>() << " is not a known task type! Skipping." << std::endl;
            }
            i++;
        }
        
        for (auto task : MetaTasks ) {
            auto filePath = task;
            
            std::cout << "MetaTasks are not currently supported by HPCMetaMorpheusTask\n";
        }

        std::vector<std::string> startingRawFilenameList;
        for ( auto b: Spectra) {
            startingRawFilenameList.push_back(std::filesystem::absolute(b));
        }
        
        std::vector<DbForTask*> startingXmlDbFilenameList;
        for ( auto b: Databases ) {
            auto newdb = new DbForTask(std::filesystem::absolute(b), IsContaminant(b));
            startingXmlDbFilenameList.push_back(newdb);       
        }
        
        std::string outputFolder = OutputFolder;
        if (outputFolder == "")
        {
            char dates[100];
            time_t curr_time;
            tm *curr_tm;

            time(&curr_time);
            curr_tm = localtime(&curr_time);
            strftime(dates, 100, "%Y-%m-%d-%H-%M-%S", curr_tm);
            std::string date_string(dates);
            
            auto pathOfFirstSpectraFile = FileSystem::getDirectoryName(startingRawFilenameList.front());
            outputFolder = pathOfFirstSpectraFile + "/" + date_string;
        }

        if (!std::filesystem::exists(outputFolder) ) {
            std::filesystem::create_directory(outputFolder);
        }
        
        EverythingRunnerEngine *a = new EverythingRunnerEngine(taskList, startingRawFilenameList,
                                                               startingXmlDbFilenameList, outputFolder);
        
        try
        {
            a->Run();
        }
        catch (const std::runtime_error &e)
        {
            std::cout << "Run failed, Exception: " << e.what() << std::endl;
        }
            
        delete a;
    }

}
