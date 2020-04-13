#include "TestToml.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
#include "../TaskLayer/XLSearchTask/TaskLayer.XLSearchTask.h"
#include "../TaskLayer/FileSpecificParameters.h"
#include "../EngineLayer/CommonParameters.h"

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
    
    std::cout << ++i << ". TestToml::TestTomlForSpecficFiles" << std::endl;
    Test::TestToml::TestTomlForSpecficFiles();
    
    return 0;
}
    
namespace Test
{

	void TestToml::TestTomlFunction()
	{
		SearchTask *searchTask = new SearchTask();
		CommonParameters tempVar(, , , = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, new PpmTolerance(666), , , , , {("e", "f"), ("g", "h")}, new std::vector<(std::string, std::string)*> {("a", "b"), ("c", "d")});
		searchTask->setCommonParameters(&tempVar);
		Toml::WriteFile(searchTask, "SearchTask.toml", MetaMorpheusTask::tomlConfig);
		auto searchTaskLoaded = Toml::ReadFile<SearchTask*>("SearchTask.toml", MetaMorpheusTask::tomlConfig);

		Assert::AreEqual(searchTask->getCommonParameters()->getDeconvolutionMassTolerance()->ToString(), searchTaskLoaded->CommonParameters.DeconvolutionMassTolerance.ToString());
		Assert::AreEqual(searchTask->getCommonParameters()->getProductMassTolerance()->ToString(), searchTaskLoaded->CommonParameters.ProductMassTolerance.ToString());
		Assert::AreEqual(searchTask->getCommonParameters()->getPrecursorMassTolerance()->ToString(), searchTaskLoaded->CommonParameters.PrecursorMassTolerance.ToString());

		Assert::AreEqual(searchTask->getCommonParameters()->ListOfModsFixed->Count(), searchTaskLoaded->CommonParameters.ListOfModsFixed->Count());
		Assert::AreEqual(searchTask->getCommonParameters()->ListOfModsFixed.First().Item1, searchTaskLoaded->CommonParameters.ListOfModsFixed.First().Item1);
		Assert::AreEqual(searchTask->getCommonParameters()->ListOfModsFixed.First().Item2, searchTaskLoaded->CommonParameters.ListOfModsFixed.First().Item2);

		Assert::AreEqual(searchTask->getCommonParameters()->ListOfModsVariable->Count(), searchTaskLoaded->CommonParameters.ListOfModsVariable->Count());

		Assert::AreEqual(searchTask->getSearchParameters()->getMassDiffAcceptorType(), searchTaskLoaded->SearchParameters.MassDiffAcceptorType);
		Assert::AreEqual(searchTask->getSearchParameters()->getCustomMdac(), searchTaskLoaded->SearchParameters.CustomMdac);

		std::string outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestConsistency)");
		std::string myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\PrunedDbSpectra.mzml)");
		std::string myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\DbForPrunedDb.fasta)");

		auto engine = new EverythingRunnerEngine(std::vector<(std::string, MetaMorpheusTask)*> {("Search", searchTask)}, std::vector<std::string> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		engine->Run();
		auto engineToml = new EverythingRunnerEngine(std::vector<(std::string, MetaMorpheusTask)*> {("SearchTOML", searchTaskLoaded)}, std::vector<std::string> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		engineToml->Run();

		auto results = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(Search\AllPSMs.psmtsv)"));
		auto resultsToml = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(SearchTOML\AllPSMs.psmtsv)"));
		Assert::That(results.SequenceEqual(resultsToml));

		CalibrationTask *calibrationTask = new CalibrationTask();
		Toml::WriteFile(calibrationTask, "CalibrationTask.toml", MetaMorpheusTask::tomlConfig);
		auto calibrationTaskLoaded = Toml::ReadFile<CalibrationTask*>("CalibrationTask.toml", MetaMorpheusTask::tomlConfig);

		GptmdTask *gptmdTask = new GptmdTask();
		Toml::WriteFile(gptmdTask, "GptmdTask.toml", MetaMorpheusTask::tomlConfig);
		auto gptmdTaskLoaded = Toml::ReadFile<GptmdTask*>("GptmdTask.toml", MetaMorpheusTask::tomlConfig);

		auto gptmdEngine = new EverythingRunnerEngine(std::vector<(std::string, MetaMorpheusTask)*> {("GPTMD", gptmdTask)}, std::vector<std::string> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		gptmdEngine->Run();
		auto gptmdEngineToml = new EverythingRunnerEngine(std::vector<(std::string, MetaMorpheusTask)*> {("GPTMDTOML", gptmdTaskLoaded)}, std::vector<std::string> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		gptmdEngineToml->Run();

		auto gptmdResults = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(GPTMD\GPTMD_Candidates.psmtsv)"));
		auto gptmdResultsToml = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(GPTMDTOML\GPTMD_Candidates.psmtsv)"));

		Assert::That(gptmdResults.SequenceEqual(gptmdResultsToml));

		XLSearchTask *xLSearchTask = new XLSearchTask();
		Toml::WriteFile(xLSearchTask, "XLSearchTask.toml", MetaMorpheusTask::tomlConfig);
		auto xLSearchTaskLoaded = Toml::ReadFile<XLSearchTask*>("XLSearchTask.toml", MetaMorpheusTask::tomlConfig);

		std::string myFileXl = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XlTestData\BSA_DSSO_ETchD6010.mgf)");
		std::string myDatabaseXl = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XlTestData\BSA.fasta)");

		auto xlEngine = new EverythingRunnerEngine(std::vector<(std::string, MetaMorpheusTask)*> {("XLSearch", xLSearchTask)}, std::vector<std::string> {myFileXl}, std::vector<DbForTask*> {new DbForTask(myDatabaseXl, false)}, outputFolder);
		xlEngine->Run();
		auto xlEngineToml = new EverythingRunnerEngine(std::vector<(std::string, MetaMorpheusTask)*> {("XLSearchTOML", xLSearchTaskLoaded)}, std::vector<std::string> {myFileXl}, std::vector<DbForTask*> {new DbForTask(myDatabaseXl, false)}, outputFolder);
		xlEngineToml->Run();

		auto xlResults = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(XLSearch\XL_Intralinks.tsv)"));
		auto xlResultsToml = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(XLSearchTOML\XL_Intralinks.tsv)"));

		Assert::That(xlResults.SequenceEqual(xlResultsToml));
		Directory::Delete(outputFolder, true);
		File::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(GptmdTask.toml)"));
		File::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XLSearchTask.toml)"));
		File::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(SearchTask.toml)"));
		File::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(CalibrationTask.toml)"));

		delete xlEngineToml;
		delete xlEngine;
                delete xLSearchTask;
		delete gptmdEngineToml;
		delete gptmdEngine;
                delete gptmdTask;
                delete calibrationTask;
		delete engineToml;
		delete engine;
                delete searchTask;
	}

	void TestToml::TestTomlForSpecficFiles()
	{
		auto fileSpecificToml = Toml::ReadFile(FileSystem::combine(TestContext::CurrentContext->TestDirectory, "testFileSpecfic.toml"), MetaMorpheusTask::tomlConfig);
		auto tomlSettingsList = fileSpecificToml->ToDictionary([&] (std::any p)
		{
			p::Key;
		});
		Assert::AreEqual(tomlSettingsList["Protease"].Value->Get<std::string>(), "Asp-N");
		Assert::IsFalse(tomlSettingsList.find("maxMissedCleavages") != tomlSettingsList.end());
		Assert::IsFalse(tomlSettingsList.find("InitiatorMethionineBehavior") != tomlSettingsList.end());

		FileSpecificParameters *f = new FileSpecificParameters(fileSpecificToml);

		Assert::AreEqual("Asp-N", f->getProtease()->Name);
		Assert::IsNull(f->getMaxMissedCleavages());

		CommonParameters tempVar();
		CommonParameters *c = MetaMorpheusTask::SetAllFileSpecificCommonParams(&tempVar, f);

		Assert::AreEqual("Asp-N", c->getDigestionParams()->Protease->Name);
		Assert::AreEqual(2, c->getDigestionParams()->MaxMissedCleavages);

                delete f;
	}
}
