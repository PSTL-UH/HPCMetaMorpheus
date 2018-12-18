#include "Program.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineFinishedEventArgs.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/EventArgs/SingleTaskEventArgs.h"
#include "../EngineLayer/EventArgs/SingleFileEventArgs.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
#include "../TaskLayer/XLSearchTask/TaskLayer.XLSearchTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"

using namespace EngineLayer;
using namespace Fclp;
using namespace Nett;
using namespace Proteomics;
using namespace TaskLayer;

namespace MetaMorpheusCommandLine
{

bool Program::InProgress = false;
System::CodeDom::Compiler::IndentedTextWriter *Program::MyWriter = new System::CodeDom::Compiler::IndentedTextWriter(Console::Out, L"\t");

	void Program::Main(std::vector<std::wstring> &args)
	{
		std::wcout << L"Welcome to MetaMorpheus" << std::endl;
		std::wcout << GlobalVariables::getMetaMorpheusVersion() << std::endl;

		auto p = new FluentCommandLineParser<ApplicationArguments*>();

		p->Setup([&] (std::any arg)
		{
			arg::Tasks;
		}).As(L't', L"tasks").SetDefault(std::vector<std::wstring>());

		p->Setup([&] (std::any arg)
		{
			arg::OutputFolder;
		}).As(L'o', L"outputFolder").SetDefault(nullptr);

		p->Setup([&] (std::any arg)
		{
			arg::MetaTasks;
		}).As(L'm', L"meta-task").SetDefault(std::vector<std::wstring>());

		p->Setup([&] (std::any arg)
		{
			arg::Spectra;
		}).As(L's', L"spectra").Required();

		p->Setup([&] (std::any arg)
		{
			arg::Databases;
		}).As(L'd', L"databases").Required();

		auto result = p->Parse(args);

		if (p->Object.MetaTasks->Count != 0 || p->Object.Tasks->Count != 0)
		{
			if (!result->HasErrors)
			{
				MetaMorpheusEngine::WarnHandler->addListener(L"WarnHandler", [&] (std::any sender, StringEventArgs* e) {WarnHandler(sender, e);});
				MetaMorpheusEngine::OutProgressHandler->addListener(L"MyEngine_outProgressHandler", [&] (std::any sender, ProgressEventArgs* e) {MyEngine_outProgressHandler(sender, e);});
				MetaMorpheusEngine::StartingSingleEngineHander->addListener(L"MyEngine_startingSingleEngineHander", [&] (std::any sender, SingleEngineEventArgs* e) {MyEngine_startingSingleEngineHander(sender, e);});
				MetaMorpheusEngine::FinishedSingleEngineHandler->addListener(L"MyEngine_finishedSingleEngineHandler", [&] (std::any sender, SingleEngineFinishedEventArgs* e) {MyEngine_finishedSingleEngineHandler(sender, e);});

				MetaMorpheusTask::WarnHandler->addListener(L"WarnHandler", [&] (std::any sender, StringEventArgs* e) {WarnHandler(sender, e);});
				MetaMorpheusTask::LogHandler->addListener(L"LogHandler", [&] (std::any sender, StringEventArgs* e) {LogHandler(sender, e);});
				MetaMorpheusTask::StartingSingleTaskHander->addListener(L"MyTaskEngine_startingSingleTaskHander", [&] (std::any sender, SingleTaskEventArgs* e) {MyTaskEngine_startingSingleTaskHander(sender, e);});
				MetaMorpheusTask::FinishedSingleTaskHandler->addListener(L"MyTaskEngine_finishedSingleTaskHandler", [&] (std::any sender, SingleTaskEventArgs* e) {MyTaskEngine_finishedSingleTaskHandler(sender, e);});
				MetaMorpheusTask::FinishedWritingFileHandler->addListener(L"MyTaskEngine_finishedWritingFileHandler", [&] (std::any sender, SingleFileEventArgs* e) {MyTaskEngine_finishedWritingFileHandler(sender, e);});

				for (auto db : p->Object.Databases)
				{
					if (Path::GetExtension(db) != L".fasta")
					{
						GlobalVariables::AddMods(UsefulProteomicsDatabases::ProteinDbLoader::GetPtmListFromProteinXml(db).OfType<Modification*>(), true);

						// print any error messages reading the mods to the console
						for (auto error : GlobalVariables::ErrorsReadingMods)
						{
							std::wcout << error << std::endl;
						}
						GlobalVariables::ErrorsReadingMods.clear();
					}
				}

				std::vector<(std::wstring, MetaMorpheusTask)*> taskList;

				for (int i = 0; i < p->Object.Tasks->Count; i++)
				{
					auto filePath = p->Object.Tasks[i];

					auto uhum = Toml::ReadFile(filePath, MetaMorpheusTask::tomlConfig);

//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//					switch (uhum.Get<string>("TaskType"))
//ORIGINAL LINE: case "Search":
					if (uhum->Get<std::wstring>(L"TaskType") == L"Search")
					{
							auto ye1 = Toml::ReadFile<SearchTask*>(filePath, MetaMorpheusTask::tomlConfig);
							taskList.push_back((L"Task" + std::to_wstring(i + 1) + L"SearchTask", ye1));

					}
//ORIGINAL LINE: case "Calibrate":
					else if (uhum->Get<std::wstring>(L"TaskType") == L"Calibrate")
					{
							auto ye2 = Toml::ReadFile<CalibrationTask*>(filePath, MetaMorpheusTask::tomlConfig);
							taskList.push_back((L"Task" + std::to_wstring(i + 1) + L"CalibrationTask", ye2));

					}
//ORIGINAL LINE: case "Gptmd":
					else if (uhum->Get<std::wstring>(L"TaskType") == L"Gptmd")
					{
							auto ye3 = Toml::ReadFile<GptmdTask*>(filePath, MetaMorpheusTask::tomlConfig);
							taskList.push_back((L"Task" + std::to_wstring(i + 1) + L"GptmdTask", ye3));

					}
//ORIGINAL LINE: case "XLSearch":
					else if (uhum->Get<std::wstring>(L"TaskType") == L"XLSearch")
					{
							auto ye4 = Toml::ReadFile<XLSearchTask*>(filePath, MetaMorpheusTask::tomlConfig);
							taskList.push_back((L"Task" + std::to_wstring(i + 1) + L"XLSearchTask", ye4));

					}
					else
					{
							std::wcout << uhum->Get<std::wstring>(L"TaskType") << L" is not a known task type! Skipping." << std::endl;
					}
				}

				for (int i = 0; i < p->Object.MetaTasks->Count; i++)
				{
					auto filePath = p->Object.MetaTasks[i];
					auto uhum = Toml::ReadFile(filePath, MetaMorpheusTask::tomlConfig);
//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//					switch (uhum.Get<string>("TaskType"))
//ORIGINAL LINE: case "Search":
					if (uhum->Get<std::wstring>(L"TaskType") == L"Search")
					{
							std::wcout << L"Search tasks are individual tasks. Please use -t for task instead of -m. Skipping." << std::endl;

					}
//ORIGINAL LINE: case "Calibrate":
					else if (uhum->Get<std::wstring>(L"TaskType") == L"Calibrate")
					{
							std::wcout << L"Calibrate tasks are individual tasks. Please use -t for task instead of -m. Skipping." << std::endl;

					}
//ORIGINAL LINE: case "Gptmd":
					else if (uhum->Get<std::wstring>(L"TaskType") == L"Gptmd")
					{
							std::wcout << L"Gptmd tasks are individual tasks. Please use -t for task instead of -m. Skipping." << std::endl;

					}
//ORIGINAL LINE: case "XLSearch":
					else if (uhum->Get<std::wstring>(L"TaskType") == L"XLSearch")
					{
							std::wcout << L"XLSearch tasks are individual tasks. Please use -t for task instead of -m. Skipping." << std::endl;

					}
					else
					{
							std::wcout << uhum->Get<std::wstring>(L"TaskType") << L" is not a known task type! Skipping." << std::endl;
					}
				}

				std::vector<std::wstring> startingRawFilenameList = p->Object.Spectra->Select([&] (std::any b)
				{
					FileSystem::getFullPath(b);
				}).ToList();
				std::vector<DbForTask*> startingXmlDbFilenameList = p->Object.Databases->Select([&] (std::any b)
				{
					new DbForTask(FileSystem::getFullPath(b), IsContaminant(b));
				}).ToList();

				std::wstring outputFolder = p->Object.OutputFolder;
				if (outputFolder == L"")
				{
					auto pathOfFirstSpectraFile = FileSystem::getDirectoryName(startingRawFilenameList.front());
					outputFolder = FileSystem::combine(pathOfFirstSpectraFile, LR"($DATETIME)");
				}

				EverythingRunnerEngine *a = new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList, outputFolder);

				try
				{
					a->Run();
				}
				catch (const std::runtime_error &e)
				{
					while (e.InnerException != nullptr)
					{
						e = e.InnerException;
					}
					auto message = L"Run failed, Exception: " + e.what();
					std::wcout << message << std::endl;
				}

				delete a;
			}
			else
			{
				std::wcout << L"Error Text:" << result->ErrorText << std::endl;
			}
		}
		else
		{
			std::wcout << L"Error Text: No toml file was specified. Use -t for tasks or -m for meta-tasks." << std::endl;
		}

		delete p;
	}

	void Program::WriteMultiLineIndented(const std::wstring &toWrite)
	{
		std::vector<std::wstring> tokens = Regex::Split(toWrite, LR"(\r?\n|\r)");
		for (auto str : tokens)
		{
			MyWriter->WriteLine(str);
		}
	}

	bool Program::IsContaminant(const std::wstring &b)
	{
		if (StringHelper::toUpper(b).find(StringHelper::toUpper(L"contaminant")) != std::wstring::npos || StringHelper::toUpper(b).find(L"CRAP") != std::wstring::npos)
		{
			return true;
		}
		return false;
	}

	void Program::MyTaskEngine_startingSingleTaskHander(std::any sender, SingleTaskEventArgs *e)
	{
		if (InProgress)
		{
			MyWriter->WriteLine();
		}
		InProgress = false;
		WriteMultiLineIndented(L"Starting task: " + e->getDisplayName());
		MyWriter->Indent = MyWriter->Indent + 1;
	}

	void Program::MyTaskEngine_finishedWritingFileHandler(std::any sender, SingleFileEventArgs *e)
	{
		if (InProgress)
		{
			MyWriter->WriteLine();
		}
		InProgress = false;
		WriteMultiLineIndented(L"Finished writing file: " + e->getWrittenFile());
	}

	void Program::MyTaskEngine_finishedSingleTaskHandler(std::any sender, SingleTaskEventArgs *e)
	{
		if (InProgress)
		{
			MyWriter->WriteLine();
		}
		InProgress = false;
		MyWriter->Indent = MyWriter->Indent - 1;
		WriteMultiLineIndented(L"Finished task: " + e->getDisplayName());
	}

	void Program::MyEngine_startingSingleEngineHander(std::any sender, SingleEngineEventArgs *e)
	{
		if (InProgress)
		{
			MyWriter->WriteLine();
		}
		InProgress = false;
		WriteMultiLineIndented(L"Starting engine: " + e->getMyEngine()->GetType()->Name + L" " + e->getMyEngine()->GetId());
		MyWriter->Indent = MyWriter->Indent + 1;
	}

	void Program::MyEngine_finishedSingleEngineHandler(std::any sender, SingleEngineFinishedEventArgs *e)
	{
		if (InProgress)
		{
			MyWriter->WriteLine();
		}
		InProgress = false;
		WriteMultiLineIndented(L"Engine results: " + e);
		MyWriter->Indent = MyWriter->Indent - 1;
		WriteMultiLineIndented(L"Finished engine: " + e->MyResults->getMyEngine()->GetType()->Name + L" " + e->MyResults->getMyEngine()->GetId());
	}

	void Program::MyEngine_outProgressHandler(std::any sender, ProgressEventArgs *e)
	{
		MyWriter->Write(std::to_wstring(e->NewProgress) + L" ");
		InProgress = true;
	}

	void Program::WarnHandler(std::any sender, StringEventArgs *e)
	{
		if (InProgress)
		{
			MyWriter->WriteLine();
		}
		InProgress = false;
		WriteMultiLineIndented(L"WARN: " + e->getS());
	}

	void Program::LogHandler(std::any sender, StringEventArgs *e)
	{
		if (InProgress)
		{
			MyWriter->WriteLine();
		}
		InProgress = false;
		WriteMultiLineIndented(L"Log: " + e->getS());
	}

	std::vector<std::wstring> Program::ApplicationArguments::getTasks() const
	{
		return privateTasks;
	}

	void Program::ApplicationArguments::setTasks(const std::vector<std::wstring> &value)
	{
		privateTasks = value;
	}

	std::vector<std::wstring> Program::ApplicationArguments::getDatabases() const
	{
		return privateDatabases;
	}

	void Program::ApplicationArguments::setDatabases(const std::vector<std::wstring> &value)
	{
		privateDatabases = value;
	}

	std::vector<std::wstring> Program::ApplicationArguments::getSpectra() const
	{
		return privateSpectra;
	}

	void Program::ApplicationArguments::setSpectra(const std::vector<std::wstring> &value)
	{
		privateSpectra = value;
	}

	std::vector<std::wstring> Program::ApplicationArguments::getMetaTasks() const
	{
		return privateMetaTasks;
	}

	void Program::ApplicationArguments::setMetaTasks(const std::vector<std::wstring> &value)
	{
		privateMetaTasks = value;
	}

	std::wstring Program::ApplicationArguments::getOutputFolder() const
	{
		return privateOutputFolder;
	}

	void Program::ApplicationArguments::setOutputFolder(const std::wstring &value)
	{
		privateOutputFolder = value;
	}
}
