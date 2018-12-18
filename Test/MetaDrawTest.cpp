#include "MetaDrawTest.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../EngineLayer/MetaDraw/MetaDrawPsm.h"
#include "../EngineLayer/MetaDraw/TsvResultReader.h"
#include "../TaskLayer/XLSearchTask/TaskLayer.XLSearchTask.h"

using namespace EngineLayer;
using namespace NUnit::Framework;
using namespace TaskLayer;

namespace Test
{

	void MetaDrawTest::TestMetaDrawReadPsmFile()
	{
		SearchTask *searchTask = new SearchTask();

		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\PrunedDbSpectra.mzml)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\DbForPrunedDb.fasta)");
		std::wstring folderPath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestMetaDrawReadPsmFile)");

		DbForTask *db = new DbForTask(myDatabase, false);
		FileSystem::createDirectory(folderPath);

		searchTask->RunTask(folderPath, {db}, {myFile}, L"metadraw");
		std::wstring psmFile = Directory::GetFiles(folderPath).First([&] (std::any f)
		{
			f->Contains(L"AllPSMs.psmtsv");
		});

		std::vector<string> warnings;
		std::vector<MetaDrawPsm*> parsedPsms = TsvResultReader::ReadTsv(psmFile, warnings);

		Assert::AreEqual(11, parsedPsms.size());
		Assert::AreEqual(0, warnings->Count);

		Directory::Delete(folderPath, true);

//C# TO C++ CONVERTER TODO TASK: A 'delete db' statement was not added since db was passed to a method or constructor. Handle memory management manually.
		delete searchTask;
	}

	void MetaDrawTest::TestMetaDrawReadCrossPsmFile()
	{
		XLSearchTask *searchTask = new XLSearchTask();
		searchTask->getXlSearchParameters()->setCrosslinkerType(EngineLayer::CrosslinkSearch::CrosslinkerType::DSS);

		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XlTestData\BSA_DSS_23747.mzML)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XlTestData\BSA.fasta)");
		std::wstring folderPath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestMetaDrawReadPsmFile)");

		DbForTask *db = new DbForTask(myDatabase, false);
		FileSystem::createDirectory(folderPath);

		searchTask->RunTask(folderPath, {db}, {myFile}, L"metadraw");

		std::wstring psmFile = Directory::GetFiles(folderPath).First([&] (std::any f)
		{
			f->Contains(L"XL_Intralinks.tsv");
		});

		std::vector<string> warnings;
		std::vector<MetaDrawPsm*> parsedPsms = TsvResultReader::ReadTsv(psmFile, warnings);

		Directory::Delete(folderPath, true);

//C# TO C++ CONVERTER TODO TASK: A 'delete db' statement was not added since db was passed to a method or constructor. Handle memory management manually.
		delete searchTask;
	}
}
