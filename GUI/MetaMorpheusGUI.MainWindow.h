#pragma once

#include <string>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <any>
#include <functional>
#include "stringhelper.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace MetaMorpheusGUI { class RawDataForDataGrid; }
namespace MetaMorpheusGUI { class ProteinDbForDataGrid; }
namespace MetaMorpheusGUI { class PreRunTask; }
namespace MetaMorpheusGUI { class InRunTask; }
namespace MetaMorpheusGUI { class GuiGlobalParams; }
namespace EngineLayer { class StringEventArgs; }
namespace TaskLayer { class XmlForTaskListEventArgs; }
namespace EngineLayer { class StringListEventArgs; }
namespace TaskLayer { class SingleTaskEventArgs; }
namespace TaskLayer { class MetaMorpheusTask; }
namespace EngineLayer { class ProgressEventArgs; }
namespace EngineLayer { class SingleFileEventArgs; }

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Nett;
using namespace Newtonsoft::Json::Linq;
using namespace Proteomics;
using namespace TaskLayer;

namespace MetaMorpheusGUI
{
	/// <summary>
	/// Interaction logic for MainWindow.xaml
	/// </summary>
	class MainWindow : public Window
	{
	private:
		static std::wstring privateNewestKnownVersion;

		ObservableCollection<RawDataForDataGrid*> *const SpectraFilesObservableCollection = new ObservableCollection<RawDataForDataGrid*>();
		ObservableCollection<ProteinDbForDataGrid*> *const ProteinDbObservableCollection = new ObservableCollection<ProteinDbForDataGrid*>();
		ObservableCollection<PreRunTask*> *const StaticTasksObservableCollection = new ObservableCollection<PreRunTask*>();
		ObservableCollection<RawDataForDataGrid*> *const SelectedRawFiles = new ObservableCollection<RawDataForDataGrid*>();
		ObservableCollection<InRunTask*> *DynamicTasksObservableCollection;
		bool WarnedAboutThermoAlready = false;

	public:
		virtual ~MainWindow()
		{
			delete SpectraFilesObservableCollection;
			delete ProteinDbObservableCollection;
			delete StaticTasksObservableCollection;
			delete SelectedRawFiles;
			delete DynamicTasksObservableCollection;
			delete GuiGlobalParams;
		}

		MainWindow();

	private:
		FlowDocument *YoutubeWikiNotification();

	public:
		static std::wstring getNewestKnownVersion();
		static void setNewestKnownVersion(const std::wstring &value);

		GuiGlobalParams *GuiGlobalParams = new GuiGlobalParams();

	private:
		static void GetVersionNumbersFromWeb();

		void MyWindow_Loaded(std::any sender, RoutedEventArgs *e);

		void EverythingRunnerEngine_FinishedWritingAllResultsFileHandler(std::any sender, StringEventArgs *e);

		void GuiWarnHandler(std::any sender, StringEventArgs *e);

		void MyTaskEngine_FinishedDataFileHandler(std::any sender, StringEventArgs *s);

		void MyTaskEngine_StartingDataFileHandler(std::any sender, StringEventArgs *s);

		void AddNewDB(std::any sender, XmlForTaskListEventArgs *e);

		void AddNewSpectra(std::any sender, StringListEventArgs *e);

		void AddNewFileSpecificToml(std::any sender, StringListEventArgs *e);

		void UpdateOutputFolderTextbox();

		void Po_startingSingleTaskHander(std::any sender, SingleTaskEventArgs *s);

		void Po_finishedSingleTaskHandler(std::any sender, SingleTaskEventArgs *s);

		void ClearSpectraFiles_Click(std::any sender, RoutedEventArgs *e);

		void OpenOutputFolder_Click(std::any sender, RoutedEventArgs *e);

		void AddProteinDatabase_Click(std::any sender, RoutedEventArgs *e);

		void AddSpectraFile_Click(std::any sender, RoutedEventArgs *e);

		void Window_Drop(std::any sender, DragEventArgs *e);

		void AddAFile(const std::wstring &draggedFilePath);

		void AddTaskToCollection(MetaMorpheusTask *ye);

		// handles double-clicking on a data grid row
		void Row_DoubleClick(std::any sender, MouseButtonEventArgs *e);

		void RunAllTasks_Click(std::any sender, RoutedEventArgs *e);

		void EverythingRunnerExceptionHandler(Task *obj);

		void ClearTasks_Click(std::any sender, RoutedEventArgs *e);

		void UpdateTaskGuiStuff();

		void UpdateSpectraFileGuiStuff();

		void AddSearchTaskButton_Click(std::any sender, RoutedEventArgs *e);

		void AddCalibrateTaskButton_Click(std::any sender, RoutedEventArgs *e);

		void AddGPTMDTaskButton_Click(std::any sender, RoutedEventArgs *e);

		void BtnAddCrosslinkSearch_Click(std::any sender, RoutedEventArgs *e);

		// deletes the selected task
		void DeleteSelectedTask(std::any sender, RoutedEventArgs *e);

		// move the task up or down in the GUI
		void MoveSelectedTask(std::any sender, RoutedEventArgs *e, bool moveTaskUp);

		// handles keyboard input in the main window
		void Window_KeyDown(std::any sender, KeyEventArgs *e);

		void NewCollectionHandler(std::any sender, StringEventArgs *s);

		void NewoutLabelStatus(std::any sender, StringEventArgs *s);

		// update progress bar for task/file
		void NewoutProgressBar(std::any sender, ProgressEventArgs *s);

		void NewRefreshBetweenTasks(std::any sender, EventArgs *e);

		void NewSuccessfullyStartingAllTasks(std::any sender, EventArgs *e);

		void NewSuccessfullyFinishedAllTasks(std::any sender, StringEventArgs *e);

		void NewSuccessfullyFinishedWritingFile(std::any sender, SingleFileEventArgs *v);

		void ClearXML_Click(std::any sender, RoutedEventArgs *e);

		void ResetTasksButton_Click(std::any sender, RoutedEventArgs *e);

		void TasksTreeView_MouseDoubleClick(std::any sender, MouseButtonEventArgs *e);

		void LoadTaskButton_Click(std::any sender, RoutedEventArgs *e);

		//run if fileSpecificParams are changed from GUI
		void UpdateFileSpecificParamsDisplay(std::vector<std::wstring> &tomlLocations);

		//run if data file has just been added with and checks for Existing fileSpecficParams
		void UpdateFileSpecificParamsDisplayJustAdded(const std::wstring &tomlLocation);

		void AddSelectedSpectra(std::any sender, RoutedEventArgs *e);

		void RemoveSelectedSpectra(std::any sender, RoutedEventArgs *e);

		void MenuItem_Click(std::any sender, RoutedEventArgs *e);

		void MenuItem_YouTube(std::any sender, RoutedEventArgs *e);

		void MenuItem_ProteomicsNewsBlog(std::any sender, RoutedEventArgs *e);

		void MenuItem_Click_1(std::any sender, RoutedEventArgs *e);

		bool DatabaseExists(ObservableCollection<ProteinDbForDataGrid*> *pDOC, ProteinDbForDataGrid *uuu);

		bool SpectraFileExists(ObservableCollection<RawDataForDataGrid*> *rDOC, RawDataForDataGrid *zzz);

		void CancelButton_Click(std::any sender, RoutedEventArgs *e);

		void ChangeFileParameters_Click(std::any sender, RoutedEventArgs *e);

		void BtnQuantSet_Click(std::any sender, RoutedEventArgs *e);

		void MenuItem_Click_2(std::any sender, RoutedEventArgs *e);

		void MenuItem_Click_4(std::any sender, RoutedEventArgs *e);

		void MenuItem_Click_5(std::any sender, RoutedEventArgs *e);

		void MenuItem_Twitter(std::any sender, RoutedEventArgs *e);

		void MenuItem_Slack(std::any sender, RoutedEventArgs *e);

		void MenuItem_Click_6(std::any sender, RoutedEventArgs *e);

		void MenuItem_Click_3(std::any sender, RoutedEventArgs *e);

		void MenuItem_Click_7(std::any sender, RoutedEventArgs *e);

		void MetaDrawMenu_Click(std::any sender, RoutedEventArgs *e);

		void PrintErrorsReadingMods();

		void AddContaminantXML_Click(std::any sender, RoutedEventArgs *e);

		void AddCustomMod_Click(std::any sender, RoutedEventArgs *e);

		// handle window closing
		void MainWindow_Closing(std::any sender, CancelEventArgs *e);
	};
}
