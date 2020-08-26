#include "TestToml.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
//#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
#include "../TaskLayer/XLSearchTask/XLSearchTask.h"
#include "../TaskLayer/FileSpecificParameters.h"
#include "../EngineLayer/CommonParameters.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include "Assert.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace TaskLayer;


int main ( int argc, char **argv )
{
    int i=0;

    std::cout << i << ". PeriodicTableLoader" << std::endl;
    const std::string elfile="elements.dat";
    const std::string &elr=elfile;
    UsefulProteomicsDatabases::PeriodicTableLoader::Load (elr);

    std::cout << ++i << ". TestToml::TestTomlFunction" << std::endl;
    Test::TestToml::TestTomlFunction();

#ifdef LATER
    std::cout << ++i << ". TestToml::TestTomlForSpecficFiles" << std::endl;
    Test::TestToml::TestTomlForSpecficFiles();
#endif
    
    return 0;
}
    
namespace Test
{

    void TestToml::TestTomlFunction()
    {
        std::string testdir=std::filesystem::current_path().string();        
        
        std::cout << "  SearchTask" << std::endl;
        SearchTask *searchTask = new SearchTask();
        std::vector<std::tuple<std::string, std::string>> v1 = {std::make_tuple("e", "f"), std::make_tuple("g", "h")};
        std::vector<std::tuple<std::string, std::string>> v2 = {std::make_tuple("a", "b"), std::make_tuple("c", "d")};
        auto tempVar = new CommonParameters ("", DissociationType::HCD, true, true, 3, 12, true, false, 1, 5, 200, 0.01,
                                  false, true, false, false, new PpmTolerance(666), nullptr, nullptr, -1, nullptr,
                                  &v1, &v2);
        searchTask->setCommonParameters(tempVar);
        // EDGAR: set searchtype to MODERN, since classic search is not yet working properly,
        //        and disbale protein quantification
        searchTask->getSearchParameters()->setSearchType(SearchType::Modern);
        searchTask->getSearchParameters()->setDoQuantification(false);
        
        std::ofstream tomlFile;
        std::string filename  = "SearchTask.toml";
        searchTask->writeTomlConfig( filename, tomlFile);

        std::string outputFolder = testdir + "/TestConsistency";
        std::filesystem::create_directory(outputFolder);

        SearchTask* searchTaskLoaded = new SearchTask("SearchTask.toml");


        Assert::AreEqual(searchTask->getCommonParameters()->getDeconvolutionMassTolerance()->getValue(),
                         searchTaskLoaded->getCommonParameters()->getDeconvolutionMassTolerance()->getValue());
        Assert::AreEqual(searchTask->getCommonParameters()->getProductMassTolerance()->getValue(),
                         searchTaskLoaded->getCommonParameters()->getProductMassTolerance()->getValue());
        Assert::AreEqual(searchTask->getCommonParameters()->getPrecursorMassTolerance()->getValue(),
                         searchTaskLoaded->getCommonParameters()->getPrecursorMassTolerance()->getValue());
        
        Assert::AreEqual(searchTask->getCommonParameters()->getListOfModsFixed()->size(),
                         searchTaskLoaded->getCommonParameters()->getListOfModsFixed()->size());
        Assert::AreEqual(std::get<0>(searchTask->getCommonParameters()->getListOfModsFixed()->front()),
                         std::get<0>(searchTaskLoaded->getCommonParameters()->getListOfModsFixed()->front()));
        Assert::AreEqual(std::get<1>(searchTask->getCommonParameters()->getListOfModsFixed()->front()),
                         std::get<1>(searchTaskLoaded->getCommonParameters()->getListOfModsFixed()->front()));
        
        Assert::AreEqual(searchTask->getCommonParameters()->getListOfModsVariable()->size(),
                         searchTaskLoaded->getCommonParameters()->getListOfModsVariable()->size());
        Assert::IsTrue(searchTask->getSearchParameters()->getMassDiffAcceptorType() ==
                       searchTaskLoaded->getSearchParameters()->getMassDiffAcceptorType());
        Assert::AreEqual(searchTask->getSearchParameters()->getCustomMdac(),
                         searchTaskLoaded->getSearchParameters()->getCustomMdac());
        
#ifdef NOT_NOW
        std::string myFile = testdir + "/TestData/PrunedDbSpectra.mzml";
        std::string myDatabase = testdir + "/TestData/DbForPrunedDb.fasta";

        std::vector<std::tuple<std::string, MetaMorpheusTask*> > tv1 = {std::make_tuple("Search", searchTask)};
        std::vector<std::string> tv2 = {myFile};
        std::vector<DbForTask*> tv3 = {new DbForTask(myDatabase, false)};
        auto engine = new EverythingRunnerEngine(tv1, tv2, tv3, outputFolder);
        engine->Run();

        std::vector<std::tuple<std::string, MetaMorpheusTask*>> tv4 = {std::make_tuple("SearchTOML", searchTaskLoaded)};
        auto engineToml = new EverythingRunnerEngine(tv4, tv2, tv3, outputFolder);
        engineToml->Run();
        
        //auto results = File::ReadAllLines(outputFolder + "/Search/AllPSMs.psmtsv");
        std::vector<std::string> results;
        std::ifstream if1(outputFolder + "/Search/AllPSMs.psmtsv");
        if ( if1.is_open() ) {
            std::string line;
            while ( getline(if1, line ) ){
                results.push_back(line);
            }
        }
        if1.close();
        
        //auto resultsToml = File::ReadAllLines(outputFolder + "/SearchTOML/AllPSMs.psmtsv");
        std::vector<std::string> resultsToml;
        std::ifstream if2(outputFolder + "/SearchTOML/AllPSMs.psmtsv");
        if ( if2.is_open() ) {
            std::string line;
            while ( getline(if2, line ) ){
                resultsToml.push_back(line);
            }
        }
        if2.close();
        
        Assert::SequenceEqual(results, resultsToml);
#endif
        
        std::cout << "  CalibrationTask" << std::endl;
        CalibrationTask *calibrationTask = new CalibrationTask();
        filename = "CalibrationTask.toml";
        calibrationTask->writeTomlConfig(filename, tomlFile);
        auto calibrationTaskLoaded = new CalibrationTask("CalibrationTask.toml");

#ifdef LATER
        // GptmdTask will be done later
        std::cout << "  GptmdTask" << std::endl;
        GptmdTask *gptmdTask = new GptmdTask();
        Toml::WriteFile(gptmdTask, "GptmdTask.toml", MetaMorpheusTask::tomlConfig);
        GptmdTask* gptmdTaskLoaded = Toml::ReadFile<GptmdTask*>("GptmdTask.toml", MetaMorpheusTask::tomlConfig);

        std::vector<std::tuple<std::string, MetaMorpheusTask*>> tv5 = {std::make_tuple("GPTMD", gptmdTask)};
        std::vector<std::string> tv6 = {myFile};
        std::vector<DbForTask*> tv7 = {new DbForTask(myDatabase, false)};
        auto gptmdEngine = new EverythingRunnerEngine(tv5, tv6, tv7, outputFolder);
        gptmdEngine->Run();

        std::vector<std::tuple<std::string, MetaMorpheusTask*>> tv8 = {make_tuple("GPTMDTOML", gptmdTaskLoaded)};
        auto gptmdEngineToml = new EverythingRunnerEngine(tv8, tv6, tv7, outputFolder);
        gptmdEngineToml->Run();
        
        //auto gptmdResults = File::ReadAllLines(outputFolder + "/GPTMD/GPTMD_Candidates.psmtsv");
        std::vector<std::string> gptmdResults;
        std::ifstream if3(outputFolder + "/GPTMD/GPTMD_Candidates.psmtsv");
        if ( if3.is_open() ) {
            std::string line;
            while ( getline(if3, line ) ){
                gptmdResults.push_back(line);
            }
        }
        if3.close();

        //auto gptmdResultsToml = File::ReadAllLines(outputFolder + "/GPTMDTOML/GPTMD_Candidates.psmtsv");
        std::vector<std::string> gptmdResultsToml;
        std::ifstream if4(outputFolder + "/GPTMDTOML/GPTMD_Candidates.psmtsv");
        if ( if4.is_open() ) {
            std::string line;
            while ( getline(if4, line ) ){
                gptmdResultsToml.push_back(line);
            }
        }
        if4.close();
        Assert::SequenceEqual(gptmdResults, gptmdResultsToml);
#endif
        
        std::cout << "  XLSearchTask" << std::endl;
        XLSearchTask *xLSearchTask = new XLSearchTask();
        filename = "XLSearchTask.toml";
        xLSearchTask->writeTomlConfig(filename, tomlFile);

        XLSearchTask* xLSearchTaskLoaded = new XLSearchTask("XLSearchTask.toml");
        
        std::string myFileXl = testdir + "/XlTestData/BSA_DSSO_ETchD6010.mgf";
        std::string myDatabaseXl = testdir + "/XlTestData/BSA.fasta";

        std::vector<std::tuple<std::string, MetaMorpheusTask*>> tv9 = {std::make_tuple("XLSearch", xLSearchTask)};
        std::vector<std::string> tv10 = {myFileXl};
        std::vector<DbForTask*> tv11 = {new DbForTask(myDatabaseXl, false)};
        auto xlEngine = new EverythingRunnerEngine( tv9, tv10, tv11, outputFolder);
        xlEngine->Run();

        std::vector<std::tuple<std::string, MetaMorpheusTask*>> tv12 = {std::make_tuple("XLSearchTOML", xLSearchTaskLoaded)};
        auto xlEngineToml = new EverythingRunnerEngine( tv12, tv10, tv11, outputFolder);
        xlEngineToml->Run();
        
        //auto xlResults = File::ReadAllLines(outputFolder + "/XLSearch/XL_Intralinks.tsv");
        std::vector<std::string> xlResults;
        std::ifstream if5(outputFolder + "/XLSearch/XL_Intralinks.tsv");
        if ( if5.is_open() ) {
            std::string line;
            while ( getline(if5, line ) ){
                xlResults.push_back(line);
            }
        }
        if5.close();

        //auto xlResultsToml = File::ReadAllLines(outputFolder + "/XLSearchTOML/XL_Intralinks.tsv");
        std::vector<std::string> xlResultsToml;
        std::ifstream if6(outputFolder + "/XLSearchTOML/XL_Intralinks.tsv");
        if ( if6.is_open() ) {
            std::string line;
            while ( getline(if6, line ) ){
                xlResultsToml.push_back(line);
            }
        }
        if6.close();
        
        Assert::SequenceEqual(xlResults, xlResultsToml);
        
        //std::filesystem::remove_all(outputFolder);
        //std::filesystem::remove(testdir + "/GptmdTask.toml");
        //std::filesystem::remove(testdir + "/XLSearchTask.toml");
        //std::filesystem::remove(testdir + "/SearchTask.toml");
        //std::filesystem::remove(testdir + "/CalibrationTask.toml");
        
        delete xlEngineToml;
        delete xlEngine;
        delete xLSearchTask;
        //delete gptmdEngineToml;
        //delete gptmdEngine;
        //delete gptmdTask;
        delete calibrationTask;
        //delete engineToml;
        //delete engine;
        delete searchTask;
    }

#ifdef LATER
    void TestToml::TestTomlForSpecficFiles()
    {

        std::string testdir=std::experimental::filesystem::current_path().string();
        
        auto fileSpecificToml = Toml::ReadFile(testdir+ "/testFileSpecfic.toml", MetaMorpheusTask::tomlConfig);
        auto tomlSettingsList = fileSpecificToml->ToDictionary([&] (std::any p)
                                                               {
                                                                   p::Key;
                                                               });
        std::string s1 = "Asp-N";
        Assert::AreEqual(tomlSettingsList["Protease"].Value->Get<std::string>(), s1);
        Assert::IsFalse(tomlSettingsList.find("maxMissedCleavages") != tomlSettingsList.end());
        Assert::IsFalse(tomlSettingsList.find("InitiatorMethionineBehavior") != tomlSettingsList.end());
        
        FileSpecificParameters *f = new FileSpecificParameters(fileSpecificToml);
        
        Assert::AreEqual(s1, f->getProtease()->getName());
        Assert::IsFalse(f->getMaxMissedCleavages().has_value());
        
        CommonParameters *tempVar = new CommonParameters();
        CommonParameters *c = MetaMorpheusTask::SetAllFileSpecificCommonParams(tempVar, f);
        
        Assert::AreEqual(s1, c->getDigestionParams()->getProtease()->getName());
        Assert::AreEqual(2, c->getDigestionParams()->getMaxMissedCleavages());
        
        delete f;
    }
#endif

}
