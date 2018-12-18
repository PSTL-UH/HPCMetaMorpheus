#include "SlicedTest.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"
#include "SetUpTests.h"

using namespace Nett;
using namespace NUnit::Framework;
using namespace TaskLayer;

namespace Test
{

	void SlicedTest::SlicedTest1()
	{
		auto task = Toml::ReadFile<SearchTask*>(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(SlicedSearchTaskConfig.toml)"), MetaMorpheusTask::tomlConfig);

		DbForTask *db = new DbForTask(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(sliced-db.fasta)"), false);
		std::wstring raw = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(sliced-raw.mzML)");
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestSlicedTest1)");
		EverythingRunnerEngine *a = new EverythingRunnerEngine({(L"Task", task)}, {raw}, {db}, outputFolder);

		a->Run();

		auto thisTaskOutputFolder = MySetUpClass::outputFolder;

		auto peaks = FileSystem::combine(thisTaskOutputFolder, L"Task", L"AllQuantifiedPeaks.tsv");

		Assert::AreEqual(2, File::ReadLines(peaks).size()());

		auto psms = FileSystem::combine(thisTaskOutputFolder, L"Task", L"AllPSMs.psmtsv");

		Assert::AreEqual(3, File::ReadLines(psms).size()());
		auto protGroups = FileSystem::combine(thisTaskOutputFolder, L"Task", L"AllProteinGroups.tsv");

		Assert::AreEqual(2, File::ReadLines(protGroups).size()());
		Directory::Delete(outputFolder, true);

		delete a;
//C# TO C++ CONVERTER TODO TASK: A 'delete db' statement was not added since db was passed to a method or constructor. Handle memory management manually.
	}

	void SlicedTest::FaFormatTest()
	{
		auto task = Toml::ReadFile<SearchTask*>(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(SlicedSearchTaskConfig.toml)"), MetaMorpheusTask::tomlConfig);

		DbForTask *db = new DbForTask(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(sliced-db.fa)"), false);
		std::wstring raw = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(sliced-raw.mzML)");
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(FaFormatTest)");
		EverythingRunnerEngine *a = new EverythingRunnerEngine({(L"Task", task)}, {raw}, {db}, outputFolder);

		a->Run();

		auto thisTaskOutputFolder = MySetUpClass::outputFolder;

		auto peaks = FileSystem::combine(thisTaskOutputFolder, L"Task", L"AllQuantifiedPeaks.tsv");

		Assert::AreEqual(2, File::ReadLines(peaks).size()());

		auto psms = FileSystem::combine(thisTaskOutputFolder, L"Task", L"AllPSMs.psmtsv");

		Assert::AreEqual(3, File::ReadLines(psms).size()());
		auto protGroups = FileSystem::combine(thisTaskOutputFolder, L"Task", L"AllProteinGroups.tsv");

		Assert::AreEqual(2, File::ReadLines(protGroups).size()());
		Directory::Delete(outputFolder, true);

		delete a;
//C# TO C++ CONVERTER TODO TASK: A 'delete db' statement was not added since db was passed to a method or constructor. Handle memory management manually.
	}
}
