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

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Nett;
using namespace NUnit::Framework;
using namespace TaskLayer;

namespace Test
{

	void TestToml::TestTomlFunction()
	{
		SearchTask *searchTask = new SearchTask();
		CommonParameters tempVar(, , , = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, new PpmTolerance(666), , , , , {(L"e", L"f"), (L"g", L"h")}, new std::vector<(std::wstring, std::wstring)*> {(L"a", L"b"), (L"c", L"d")});
		searchTask->setCommonParameters(&tempVar);
		Toml::WriteFile(searchTask, L"SearchTask.toml", MetaMorpheusTask::tomlConfig);
		auto searchTaskLoaded = Toml::ReadFile<SearchTask*>(L"SearchTask.toml", MetaMorpheusTask::tomlConfig);

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		Assert::AreEqual(searchTask->getCommonParameters()->getDeconvolutionMassTolerance()->ToString(), searchTaskLoaded->CommonParameters.DeconvolutionMassTolerance.ToString());
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		Assert::AreEqual(searchTask->getCommonParameters()->getProductMassTolerance()->ToString(), searchTaskLoaded->CommonParameters.ProductMassTolerance.ToString());
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		Assert::AreEqual(searchTask->getCommonParameters()->getPrecursorMassTolerance()->ToString(), searchTaskLoaded->CommonParameters.PrecursorMassTolerance.ToString());

		Assert::AreEqual(searchTask->getCommonParameters()->ListOfModsFixed->Count(), searchTaskLoaded->CommonParameters.ListOfModsFixed->Count());
		Assert::AreEqual(searchTask->getCommonParameters()->ListOfModsFixed.First().Item1, searchTaskLoaded->CommonParameters.ListOfModsFixed.First().Item1);
		Assert::AreEqual(searchTask->getCommonParameters()->ListOfModsFixed.First().Item2, searchTaskLoaded->CommonParameters.ListOfModsFixed.First().Item2);

		Assert::AreEqual(searchTask->getCommonParameters()->ListOfModsVariable->Count(), searchTaskLoaded->CommonParameters.ListOfModsVariable->Count());

		Assert::AreEqual(searchTask->getSearchParameters()->getMassDiffAcceptorType(), searchTaskLoaded->SearchParameters.MassDiffAcceptorType);
		Assert::AreEqual(searchTask->getSearchParameters()->getCustomMdac(), searchTaskLoaded->SearchParameters.CustomMdac);

		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestConsistency)");
		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\PrunedDbSpectra.mzml)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\DbForPrunedDb.fasta)");

		auto engine = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"Search", searchTask)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		engine->Run();
		auto engineToml = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"SearchTOML", searchTaskLoaded)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		engineToml->Run();

		auto results = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(Search\AllPSMs.psmtsv)"));
		auto resultsToml = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(SearchTOML\AllPSMs.psmtsv)"));
		Assert::That(results.SequenceEqual(resultsToml));

		CalibrationTask *calibrationTask = new CalibrationTask();
		Toml::WriteFile(calibrationTask, L"CalibrationTask.toml", MetaMorpheusTask::tomlConfig);
		auto calibrationTaskLoaded = Toml::ReadFile<CalibrationTask*>(L"CalibrationTask.toml", MetaMorpheusTask::tomlConfig);

		GptmdTask *gptmdTask = new GptmdTask();
		Toml::WriteFile(gptmdTask, L"GptmdTask.toml", MetaMorpheusTask::tomlConfig);
		auto gptmdTaskLoaded = Toml::ReadFile<GptmdTask*>(L"GptmdTask.toml", MetaMorpheusTask::tomlConfig);

		auto gptmdEngine = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"GPTMD", gptmdTask)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		gptmdEngine->Run();
		auto gptmdEngineToml = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"GPTMDTOML", gptmdTaskLoaded)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		gptmdEngineToml->Run();

		auto gptmdResults = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(GPTMD\GPTMD_Candidates.psmtsv)"));
		auto gptmdResultsToml = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(GPTMDTOML\GPTMD_Candidates.psmtsv)"));

		Assert::That(gptmdResults.SequenceEqual(gptmdResultsToml));

		XLSearchTask *xLSearchTask = new XLSearchTask();
		Toml::WriteFile(xLSearchTask, L"XLSearchTask.toml", MetaMorpheusTask::tomlConfig);
		auto xLSearchTaskLoaded = Toml::ReadFile<XLSearchTask*>(L"XLSearchTask.toml", MetaMorpheusTask::tomlConfig);

		std::wstring myFileXl = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XlTestData\BSA_DSSO_ETchD6010.mgf)");
		std::wstring myDatabaseXl = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XlTestData\BSA.fasta)");

		auto xlEngine = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"XLSearch", xLSearchTask)}, std::vector<std::wstring> {myFileXl}, std::vector<DbForTask*> {new DbForTask(myDatabaseXl, false)}, outputFolder);
		xlEngine->Run();
		auto xlEngineToml = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"XLSearchTOML", xLSearchTaskLoaded)}, std::vector<std::wstring> {myFileXl}, std::vector<DbForTask*> {new DbForTask(myDatabaseXl, false)}, outputFolder);
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
//C# TO C++ CONVERTER TODO TASK: A 'delete xLSearchTask' statement was not added since xLSearchTask was passed to a method or constructor. Handle memory management manually.
		delete gptmdEngineToml;
		delete gptmdEngine;
//C# TO C++ CONVERTER TODO TASK: A 'delete gptmdTask' statement was not added since gptmdTask was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete calibrationTask' statement was not added since calibrationTask was passed to a method or constructor. Handle memory management manually.
		delete engineToml;
		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTask' statement was not added since searchTask was passed to a method or constructor. Handle memory management manually.
	}

	void TestToml::TestTomlForSpecficFiles()
	{
		auto fileSpecificToml = Toml::ReadFile(FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"testFileSpecfic.toml"), MetaMorpheusTask::tomlConfig);
		auto tomlSettingsList = fileSpecificToml->ToDictionary([&] (std::any p)
		{
			p::Key;
		});
		Assert::AreEqual(tomlSettingsList[L"Protease"].Value->Get<std::wstring>(), L"Asp-N");
		Assert::IsFalse(tomlSettingsList.find(L"maxMissedCleavages") != tomlSettingsList.end());
		Assert::IsFalse(tomlSettingsList.find(L"InitiatorMethionineBehavior") != tomlSettingsList.end());

		FileSpecificParameters *f = new FileSpecificParameters(fileSpecificToml);

		Assert::AreEqual(L"Asp-N", f->getProtease()->Name);
		Assert::IsNull(f->getMaxMissedCleavages());

		CommonParameters tempVar();
		CommonParameters *c = MetaMorpheusTask::SetAllFileSpecificCommonParams(&tempVar, f);

		Assert::AreEqual(L"Asp-N", c->getDigestionParams()->Protease->Name);
		Assert::AreEqual(2, c->getDigestionParams()->MaxMissedCleavages);

//C# TO C++ CONVERTER TODO TASK: A 'delete f' statement was not added since f was passed to a method or constructor. Handle memory management manually.
	}
}
