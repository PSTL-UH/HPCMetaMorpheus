#include "EverythingRunnerEngine.h"
#include "MetaMorpheusTask.h"
#include "DbForTask.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/StringListEventArgs.h"
#include "EventArgs/XmlForTaskListEventArgs.h"

using namespace EngineLayer;

namespace TaskLayer
{

	EverythingRunnerEngine::EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> &taskList, std::vector<std::wstring> &startingRawFilenameList, std::vector<DbForTask*> &startingXmlDbFilenameList, const std::wstring &outputFolder) : TaskList(taskList)
	{
		OutputFolder = StringHelper::trim(outputFolder, L""");

		CurrentRawDataFilenameList = startingRawFilenameList;
		CurrentXmlDbFilenameList = startingXmlDbFilenameList;
	}

	void EverythingRunnerEngine::Run()
	{
		StartingAllTasks();
		auto stopWatch = new Stopwatch();
		stopWatch->Start();

		if (!CurrentRawDataFilenameList.Any())
		{
			Warn(L"No spectra files selected");
			FinishedAllTasks(L"");

			delete stopWatch;
			return;
		}

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		auto startTimeForAllFilenames = DateTime::Now.ToString(L"yyyy-MM-dd-HH-mm-ss", CultureInfo::InvariantCulture);

		OutputFolder = StringHelper::replace(OutputFolder, L"$DATETIME", startTimeForAllFilenames);

		StringBuilder *allResultsText = new StringBuilder();

		for (int i = 0; i < TaskList.size(); i++)
		{
			if (!CurrentRawDataFilenameList.Any())
			{
				Warn(L"Cannot proceed. No spectra files selected.");
				FinishedAllTasks(OutputFolder);

				delete allResultsText;
				delete stopWatch;
				return;
			}
			if (!CurrentXmlDbFilenameList.Any())
			{
				Warn(L"Cannot proceed. No protein database files selected.");
				FinishedAllTasks(OutputFolder);

				delete allResultsText;
				delete stopWatch;
				return;
			}
			auto ok = TaskList[i];

			auto outputFolderForThisTask = FileSystem::combine(OutputFolder, ok->Item1);

			if (!FileSystem::directoryExists(outputFolderForThisTask))
			{
				FileSystem::createDirectory(outputFolderForThisTask);
			}

			// Actual task running code
			auto myTaskResults = ok->Item2->RunTask(outputFolderForThisTask, CurrentXmlDbFilenameList, CurrentRawDataFilenameList, ok->Item1);

			if (myTaskResults->NewDatabases != nullptr)
			{
				CurrentXmlDbFilenameList = myTaskResults->NewDatabases;
				NewDBs(myTaskResults->NewDatabases);
			}
			if (myTaskResults->NewSpectra != nullptr)
			{
				if (CurrentRawDataFilenameList.size() == myTaskResults->NewSpectra->Count)
				{
					CurrentRawDataFilenameList = myTaskResults->NewSpectra;
				}
				else
				{
					// at least one file was not successfully calibrated
					auto successfullyCalibFiles = myTaskResults->NewSpectra->Select([&] (std::any p)
					{
						StringHelper::replace(Path::GetFileNameWithoutExtension(p), CalibrationTask::CalibSuffix, L"");
					}).ToList();
					auto origFiles = CurrentRawDataFilenameList.Select([&] (std::any p)
					{
						Path::GetFileNameWithoutExtension(p);
					}).ToList();
					auto unsuccessfullyCalibFiles = origFiles.Except(successfullyCalibFiles).ToList();
					auto unsuccessfullyCalibFilePaths = CurrentRawDataFilenameList.Where([&] (std::any p)
					{
						std::find(unsuccessfullyCalibFiles.begin(), unsuccessfullyCalibFiles.end(), Path::GetFileNameWithoutExtension(p)) != unsuccessfullyCalibFiles.end());
					}).ToList();
					CurrentRawDataFilenameList = myTaskResults->NewSpectra;
					CurrentRawDataFilenameList.insert(CurrentRawDataFilenameList.end(), unsuccessfullyCalibFilePaths.begin(), unsuccessfullyCalibFilePaths.end());
				}

				NewSpectras(myTaskResults->NewSpectra);
			}
			if (myTaskResults->NewFileSpecificTomls != nullptr)
			{
				NewFileSpecificToml(myTaskResults->NewFileSpecificTomls);
			}
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			allResultsText->appendLine(L"\r\n" + L"\r\n" + L"\r\n" + L"\r\n" + myTaskResults->ToString());
		}
		stopWatch->Stop();
		auto resultsFileName = FileSystem::combine(OutputFolder, L"allResults.txt");
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter file = new StreamWriter(resultsFileName))
		{
			StreamWriter file = StreamWriter(resultsFileName);
			file.WriteLine(L"MetaMorpheus: version " + GlobalVariables::getMetaMorpheusVersion());
			file.WriteLine(L"Total time: " + stopWatch->Elapsed);
			file.Write(allResultsText->toString());
		}
		StringEventArgs tempVar(resultsFileName, nullptr);
		FinishedWritingAllResultsFileHandler +== nullptr ? nullptr : FinishedWritingAllResultsFileHandler::Invoke(this, &tempVar);
		FinishedAllTasks(OutputFolder);

		delete allResultsText;
		delete stopWatch;
	}

	void EverythingRunnerEngine::Warn(const std::wstring &v)
	{
		StringEventArgs tempVar(v, nullptr);
		WarnHandler +== nullptr ? nullptr : WarnHandler::Invoke(this, &tempVar);
	}

	void EverythingRunnerEngine::StartingAllTasks()
	{
		StartingAllTasksEngineHandler +== nullptr ? nullptr : StartingAllTasksEngineHandler::Invoke(this, EventArgs::Empty);
	}

	void EverythingRunnerEngine::FinishedAllTasks(const std::wstring &rootOutputDir)
	{
		StringEventArgs tempVar(rootOutputDir, nullptr);
		FinishedAllTasksEngineHandler +== nullptr ? nullptr : FinishedAllTasksEngineHandler::Invoke(this, &tempVar);
	}

	void EverythingRunnerEngine::NewSpectras(std::vector<std::wstring> &newSpectra)
	{
		StringListEventArgs tempVar(newSpectra);
		NewSpectrasHandler +== nullptr ? nullptr : NewSpectrasHandler::Invoke(this, &tempVar);
	}

	void EverythingRunnerEngine::NewFileSpecificToml(std::vector<std::wstring> &newFileSpecificTomls)
	{
		StringListEventArgs tempVar(newFileSpecificTomls);
		NewFileSpecificTomlHandler +== nullptr ? nullptr : NewFileSpecificTomlHandler::Invoke(this, &tempVar);
	}

	void EverythingRunnerEngine::NewDBs(std::vector<DbForTask*> &newDatabases)
	{
		XmlForTaskListEventArgs tempVar(newDatabases);
		NewDbsHandler +== nullptr ? nullptr : NewDbsHandler::Invoke(this, &tempVar);
	}
}
