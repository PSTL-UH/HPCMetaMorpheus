#include "SearchTaskTest.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/SearchTask/MassDiffAcceptorType.h"
#include "../EngineLayer/MetaMorpheusException.h"
#include "../EngineLayer/ProteinScoringAndFdr/FdrCategory.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"

using namespace EngineLayer;
using namespace NUnit::Framework;
using namespace TaskLayer;
using namespace Proteomics::ProteolyticDigestion;
using namespace Proteomics::Fragmentation;

namespace Test
{

	void SearchTaskTest::MassDiffAceptorTest()
	{
		SearchTask *searchTask = new SearchTask();
		auto result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), searchTask->getSearchParameters()->getMassDiffAcceptorType(), searchTask->getSearchParameters()->getCustomMdac());
		Assert::That(result->getFileNameAddition() == L"1mm");

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::TwoMM, searchTask->getSearchParameters()->getCustomMdac());
		Assert::That(result->getFileNameAddition() == L"2mm");

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::ThreeMM, searchTask->getSearchParameters()->getCustomMdac());
		Assert::That(result->getFileNameAddition() == L"3mm");

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::ModOpen, searchTask->getSearchParameters()->getCustomMdac());
		Assert::That(result->getFileNameAddition() == L"-187andUp");

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Open, searchTask->getSearchParameters()->getCustomMdac());
		Assert::That(result->getFileNameAddition() == L"OpenSearch");

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Custom, L"custom ppmAroundZero 4");
		Assert::That(result->getFileNameAddition() == L"4ppmAroundZero");

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Exact, searchTask->getSearchParameters()->getCustomMdac());
		Assert::That(result->getFileNameAddition() == L"5ppmAroundZero");

		delete searchTask;
	}

	void SearchTaskTest::ParseSearchModeTest()
	{
		SearchTask *searchTask = new SearchTask();
		auto result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Custom, L"TestCustom dot 5 ppm 0,1.0029,2.0052");
		Assert::That(result->NumNotches == 3);

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Custom, L"TestCustom dot 5 da 0,1.0029,2.0052");
		Assert::That(result->NumNotches == 3);

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Custom, L"TestCustom interval [0,5];[0,5]");
		Assert::That(result->NumNotches == 1);

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Custom, L"TestCustom OpenSearch 5");
		Assert::That(result->getFileNameAddition() == L"OpenSearch");

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Custom, L"TestCustom daltonsAroundZero 5");
		Assert::That(result->getFileNameAddition() == L"5daltonsAroundZero");

		result = SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Custom, L"TestCustom ppmAroundZero 5");
		Assert::That(result->getFileNameAddition() == L"5ppmAroundZero");

		Assert::That([&] ()
		{
			SearchTask::GetMassDiffAcceptor(searchTask->getCommonParameters()->getPrecursorMassTolerance(), MassDiffAcceptorType::Custom, L"TestCustom Test 5");
		}, Throws::TypeOf<MetaMorpheusException*>());

		delete searchTask;
	}

	void SearchTaskTest::SemiSpecificFullAndSmallMatches()
	{
		SearchTask *searchTask = new SearchTask();
		SearchParameters tempVar();
		searchTask->setSearchParameters(&tempVar);
		searchTask->getSearchParameters()->setSearchType(SearchType::NonSpecific);
		searchTask->getSearchParameters()->setLocalFdrCategories(std::vector<FdrCategory> {FdrCategory::FullySpecific, FdrCategory::SemiSpecific});
		CommonParameters tempVar2(, , , = true, = 3, = 12, = true, true, , 11, , , , , , , , , , , new DigestionParams(minPeptideLength: 7, searchModeType: CleavageSpecificity::Semi, fragmentationTerminus: FragmentationTerminus::N));
		searchTask->setCommonParameters(&tempVar2);

		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\tinySemi.mgf)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\semiTest.fasta)");
		DbForTask *db = new DbForTask(myDatabase, false);

		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"TestSemiSpecificSmall", searchTask)};

		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, Environment::CurrentDirectory);
		engine->Run();

		std::wstring outputPath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestSemiSpecificSmall\AllPSMs.psmtsv)");
		auto output = File::ReadAllLines(outputPath);
		Assert::IsTrue(output.size() == 3);

		delete engine;
		delete db;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTask' statement was not added since searchTask was passed to a method or constructor. Handle memory management manually.
	}

	void SearchTaskTest::SemiSpecificTest()
	{
		std::vector<FragmentationTerminus*> terminiToTest = {FragmentationTerminus::N, FragmentationTerminus::C};
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestSemiSpecific)");
		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\PrunedDbSpectra.mzml)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\DbForPrunedDb.fasta)");
		for (auto fragTerm : terminiToTest)
		{
			SearchTask *searchTask = new SearchTask();
			SearchParameters tempVar();
			searchTask->setSearchParameters(&tempVar);
			searchTask->getSearchParameters()->setSearchType(SearchType::NonSpecific);
			searchTask->getSearchParameters()->setLocalFdrCategories(std::vector<FdrCategory> {FdrCategory::FullySpecific, FdrCategory::SemiSpecific});
			CommonParameters tempVar2(, , , = true, = 3, = 12, = true, true, , 4, , , , , , , , , , , new DigestionParams(searchModeType: CleavageSpecificity::Semi, fragmentationTerminus: fragTerm));
			searchTask->setCommonParameters(&tempVar2);

			DbForTask *db = new DbForTask(myDatabase, false);

			std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"TestSemiSpecific", searchTask)};

			auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
			engine->Run();

			std::wstring outputPath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestSemiSpecific\TestSemiSpecific\AllPSMs.psmtsv)");
			auto output = File::ReadAllLines(outputPath);
			Assert::That(output.size() == 13); //if N is only producing 11 lines, then the c is not being searched with it. //If only 12 lines, maybe missed mono issue

			delete engine;
			delete db;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTask' statement was not added since searchTask was passed to a method or constructor. Handle memory management manually.
		}
		Directory::Delete(outputFolder, true);
	}

	void SearchTaskTest::PostSearchNormalizeTest()
	{
		SearchTask *searchTask = new SearchTask();
		SearchParameters tempVar();
		searchTask->setSearchParameters(&tempVar);
		searchTask->getSearchParameters()->setNormalize(true);

		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\PrunedDbSpectra.mzml)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\DbForPrunedDb.fasta)");
		std::wstring folderPath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestNormalization)");
		std::wstring filePath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\ExperimentalDesign.tsv)");
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(filePath))
		{
			StreamWriter output = StreamWriter(filePath);
			output.WriteLine(L"FileName\tCondition\tBiorep\tFraction\tTechrep");
			output.WriteLine(std::wstring(L"PrunedDbSpectra") + L"\t" + L"condition" + L"\t" + L"1" + L"\t" + L"1" + L"\t" + L"1");
		}
		DbForTask *db = new DbForTask(myDatabase, false);
		FileSystem::createDirectory(folderPath);

		searchTask->RunTask(folderPath, {db}, {myFile}, L"normal");

		File::Delete(filePath);

		Assert::That([&] ()
		{
			searchTask->RunTask(folderPath, {db}, {myFile}, L"normal");
		}, Throws::TypeOf<MetaMorpheusException*>());
		Directory::Delete(folderPath, true);

//C# TO C++ CONVERTER TODO TASK: A 'delete db' statement was not added since db was passed to a method or constructor. Handle memory management manually.
		delete searchTask;
	}

	void SearchTaskTest::ProteinGroupsNoParsimonyTest()
	{
		SearchTask *searchTask = new SearchTask();
		SearchParameters tempVar();
		searchTask->setSearchParameters(&tempVar);
		searchTask->getSearchParameters()->setDoParsimony(false);

		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\PrunedDbSpectra.mzml)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\DbForPrunedDb.fasta)");
		std::wstring folderPath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestProteinGroupsNoParsimony)");

		DbForTask *db = new DbForTask(myDatabase, false);
		FileSystem::createDirectory(folderPath);

		searchTask->RunTask(folderPath, {db}, {myFile}, L"normal");
		Directory::Delete(folderPath, true);

//C# TO C++ CONVERTER TODO TASK: A 'delete db' statement was not added since db was passed to a method or constructor. Handle memory management manually.
		delete searchTask;
	}

	void SearchTaskTest::PrunedDbWithContaminantsTest()
	{
		SearchTask *searchTask = new SearchTask();
		SearchParameters tempVar();
		searchTask->setSearchParameters(&tempVar);
		searchTask->getSearchParameters()->setWritePrunedDatabase(true);

		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\PrunedDbSpectra.mzml)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\DbForPrunedDb.fasta)");
		std::wstring folderPath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestNormalization)");
		std::wstring filePath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\ExperimentalDesign.tsv)");

		// contaminant DB
		DbForTask *db = new DbForTask(myDatabase, true);
		FileSystem::createDirectory(folderPath);

		searchTask->RunTask(folderPath, {db}, {myFile}, L"normal");

		Assert::That(File::ReadAllLines(FileSystem::combine(folderPath, LR"(DbForPrunedDbproteinPruned.xml)")).length() > 0);
		Assert::That(File::ReadAllLines(FileSystem::combine(folderPath, LR"(DbForPrunedDbPruned.xml)")).length() > 0);
		Directory::Delete(folderPath, true);

//C# TO C++ CONVERTER TODO TASK: A 'delete db' statement was not added since db was passed to a method or constructor. Handle memory management manually.
		delete searchTask;
	}
}
