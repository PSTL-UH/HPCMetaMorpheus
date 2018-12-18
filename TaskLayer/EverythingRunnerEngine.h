#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_event.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class MetaMorpheusTask; }
namespace TaskLayer { class DbForTask; }

using namespace EngineLayer;

namespace TaskLayer
{
	class EverythingRunnerEngine
	{
	private:
		const std::vector<(std::wstring, MetaMorpheusTask)*> TaskList;
		std::wstring OutputFolder;
		std::vector<std::wstring> CurrentRawDataFilenameList;
		std::vector<DbForTask*> CurrentXmlDbFilenameList;

	public:
		EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> &taskList, std::vector<std::wstring> &startingRawFilenameList, std::vector<DbForTask*> &startingXmlDbFilenameList, const std::wstring &outputFolder);

		static TangibleEvent<EventHandler<StringEventArgs>> *FinishedWritingAllResultsFileHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler> *StartingAllTasksEngineHandler = new TangibleEvent<EventHandler>();

		static TangibleEvent<EventHandler<StringEventArgs>> *FinishedAllTasksEngineHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler<XmlForTaskListEventArgs>> *NewDbsHandler = new TangibleEvent<EventHandler<XmlForTaskListEventArgs>>();

		static TangibleEvent<EventHandler<StringListEventArgs>> *NewSpectrasHandler = new TangibleEvent<EventHandler<StringListEventArgs>>();

		static TangibleEvent<EventHandler<StringListEventArgs>> *NewFileSpecificTomlHandler = new TangibleEvent<EventHandler<StringListEventArgs>>();

		static TangibleEvent<EventHandler<StringEventArgs>> *WarnHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		void Run();

	private:
		void Warn(const std::wstring &v);

		void StartingAllTasks();

		void FinishedAllTasks(const std::wstring &rootOutputDir);

		void NewSpectras(std::vector<std::wstring> &newSpectra);

		void NewFileSpecificToml(std::vector<std::wstring> &newFileSpecificTomls);

		void NewDBs(std::vector<DbForTask*> &newDatabases);
	};
}
