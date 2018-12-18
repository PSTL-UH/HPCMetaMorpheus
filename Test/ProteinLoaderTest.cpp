#include "ProteinLoaderTest.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/FileSpecificParameters.h"
#include "../TaskLayer/MyTaskResults.h"
#include "../EngineLayer/CommonParameters.h"

using namespace EngineLayer;
using namespace NUnit::Framework;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{

	void ProteinLoaderTest::TestProteinLoad()
	{
		ProteinLoaderTask tempVar(L"");
		(&tempVar)->Run(FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestData", L"gapdh.fasta"));
		ProteinLoaderTask tempVar2(L"");
		(&tempVar2)->Run(FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestData", L"gapdh.fa"));
		ProteinLoaderTask tempVar3(L"");
		(&tempVar3)->Run(FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestData", L"gapdh.fasta.gz"));
		ProteinLoaderTask tempVar4(L"");
		(&tempVar4)->Run(FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestData", L"gapdh.fa.gz"));
	}

	ProteinLoaderTest::ProteinLoaderTask::ProteinLoaderTask(const std::wstring &x) : ProteinLoaderTask()
	{
	}

	ProteinLoaderTest::ProteinLoaderTask::ProteinLoaderTask() : MetaMorpheusTask(MyTask::Search)
	{
	}

	void ProteinLoaderTest::ProteinLoaderTask::Run(const std::wstring &dbPath)
	{
		RunSpecific(L"", {new DbForTask(dbPath, false)}, nullptr, L"", nullptr);
	}

	MyTaskResults *ProteinLoaderTest::ProteinLoaderTask::RunSpecific(const std::wstring &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::wstring> &currentRawFileList, const std::wstring &taskId, std::vector<FileSpecificParameters*> &fileSettingsList)
	{
		EngineLayer::CommonParameters tempVar();
		LoadProteins(L"", dbFilenameList, true, DecoyType::None, std::vector<std::wstring>(), &tempVar);
		return nullptr;
	}
}
