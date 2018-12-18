#include "MsDataFileTest.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"

using namespace NUnit::Framework;
using namespace TaskLayer;

namespace Test
{

	void MsDataFileTest::Setup()
	{
		Environment::CurrentDirectory = TestContext::CurrentContext->TestDirectory;
	}

	void MsDataFileTest::TestLoadAndRunMgf()
	{
		//The purpose of this test is to ensure that mgfs can be run without crashing. 
		//Whenever a new feature is added that may require things an mgf does not have, 
		//there should be a check that prevents mgfs from using that feature.
		std::wstring mgfName = LR"(TestData\ok.mgf)";
		std::wstring xmlName = LR"(TestData\okk.xml)";
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestLoadAndRunMgf)");

		SearchTask *task1 = new SearchTask();
		SearchParameters tempVar();
		task1->setSearchParameters(&tempVar);
		task1->getSearchParameters()->setDoParsimony(true);
		task1->getSearchParameters()->setDoQuantification(true);
		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"task1", task1)};
		//run!

		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mgfName}, std::vector<DbForTask*> {new DbForTask(xmlName, false)}, outputFolder);
		engine->Run();
		//Just don't crash! There should also be at least one psm at 1% FDR, but can't check for that.
		Directory::Delete(outputFolder, true);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete task1' statement was not added since task1 was passed to a method or constructor. Handle memory management manually.
	}
}
