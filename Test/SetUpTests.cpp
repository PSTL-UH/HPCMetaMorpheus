#include "SetUpTests.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"

using namespace EngineLayer;
using namespace NUnit::Framework;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{

std::wstring MySetUpClass::outputFolder = L"";
const std::wstring MySetUpClass::elementsLocation = LR"(elements.dat)";

	void MySetUpClass::Setup()
	{
		Environment::CurrentDirectory = TestContext::CurrentContext->TestDirectory;
		Loaders::LoadElements(FileSystem::combine(TestContext::CurrentContext->TestDirectory, elementsLocation));

		MetaMorpheusEngine::WarnHandler->addListener(L"WarnStatusHandler", [&] (std::any sender, StringEventArgs* e) {WarnStatusHandler(sender, e);});
		MetaMorpheusTask::WarnHandler->addListener(L"WarnStatusHandler", [&] (std::any sender, StringEventArgs* e) {WarnStatusHandler(sender, e);});

		EverythingRunnerEngine::FinishedAllTasksEngineHandler->addListener(L"SuccessfullyFinishedAllTasks", [&] (std::any sender, StringEventArgs* e) {SuccessfullyFinishedAllTasks(sender, e);});
	}

	void MySetUpClass::SuccessfullyFinishedAllTasks(std::any sender, StringEventArgs *rootOutputFolderPath)
	{
		outputFolder = rootOutputFolderPath->getS();
	}

	void MySetUpClass::WarnStatusHandler(std::any sender, StringEventArgs *e)
	{
		std::wcout << e->getS() << std::endl;
	}
}
