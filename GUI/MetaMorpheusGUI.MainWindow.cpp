#include "MetaMorpheusGUI.MainWindow.h"
#include "ForDisplayingInDataGrids/RawDataForDataGrid.h"
#include "ForDisplayingInDataGrids/ProteinDbForDataGrid.h"
#include "ForDisplayingInDataGrids/PreRunTask.h"
#include "Util/InRunTask.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../TaskLayer/EverythingRunnerEngine.h"
#include "../TaskLayer/EventArgs/XmlForTaskListEventArgs.h"
#include "../EngineLayer/EventArgs/StringListEventArgs.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/EventArgs/SingleTaskEventArgs.h"
#include "../EngineLayer/EventArgs/SingleFileEventArgs.h"
#include "../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../TaskLayer/MyFileManager.h"
#include "../TaskLayer/FileSpecificParameters.h"
#include "Util/SearchModifications.h"
#include "Util/GuiGlobalParams.h"
#include "../EngineLayer/MetaMorpheusException.h"
#include "MetaMorpheusGUI.MetaUpdater.h"
#include "Util/OutputFileForTreeView.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
#include "../TaskLayer/XLSearchTask/TaskLayer.XLSearchTask.h"
#include "MetaMorpheusGUI.SearchTaskWindow.h"
#include "MetaMorpheusGUI.CalibrateTaskWindow.h"
#include "MetaMorpheusGUI.GptmdTaskWindow.h"
#include "MetaMorpheusGUI.XLSearchTaskWindow.h"
#include "Util/ForTreeView.h"
#include "Util/CollectionForTreeView.h"
#include "MetaMorpheusGUI.GlobalSettingsWindow.h"
#include "MetaMorpheusGUI.FileSpecificParametersWindow.h"
#include "MetaMorpheusGUI.ExperimentalDesignWindow.h"
#include "MetaMorpheusGUI.MetaDraw.h"
#include "MetaMorpheusGUI.CustomModButtonWindow.h"
#include "MetaMorpheusGUI.CustomMsgBox.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Nett;
using namespace Newtonsoft::Json::Linq;
using namespace Proteomics;
using namespace TaskLayer;

namespace MetaMorpheusGUI
{

std::wstring MainWindow::privateNewestKnownVersion;

	MainWindow::MainWindow()
	{
		InitializeComponent();

		Title = L"MetaMorpheus: version " + GlobalVariables::getMetaMorpheusVersion();

		dataGridProteinDatabases->DataContext = ProteinDbObservableCollection;
		dataGridSpectraFiles->DataContext = SpectraFilesObservableCollection;
		tasksTreeView->DataContext = StaticTasksObservableCollection;

		EverythingRunnerEngine::NewDbsHandler->addListener(L"AddNewDB", [&] (std::any sender, XmlForTaskListEventArgs* e) {AddNewDB(sender, e);});
		EverythingRunnerEngine::NewSpectrasHandler->addListener(L"AddNewSpectra", [&] (std::any sender, StringListEventArgs* e) {AddNewSpectra(sender, e);});
		EverythingRunnerEngine::NewFileSpecificTomlHandler->addListener(L"AddNewFileSpecificToml", [&] (std::any sender, StringListEventArgs* e) {AddNewFileSpecificToml(sender, e);});
		EverythingRunnerEngine::StartingAllTasksEngineHandler->addListener(L"NewSuccessfullyStartingAllTasks", [&] (std::any sender, System::EventArgs e) {NewSuccessfullyStartingAllTasks(sender, e);});
		EverythingRunnerEngine::FinishedAllTasksEngineHandler->addListener(L"NewSuccessfullyFinishedAllTasks", [&] (std::any sender, StringEventArgs* e) {NewSuccessfullyFinishedAllTasks(sender, e);});
		EverythingRunnerEngine::WarnHandler->addListener(L"GuiWarnHandler", [&] (std::any sender, StringEventArgs* e) {GuiWarnHandler(sender, e);});
		EverythingRunnerEngine::FinishedWritingAllResultsFileHandler->addListener(L"EverythingRunnerEngine_FinishedWritingAllResultsFileHandler", [&] (std::any sender, StringEventArgs* e) {EverythingRunnerEngine_FinishedWritingAllResultsFileHandler(sender, e);});

		MetaMorpheusTask::StartingSingleTaskHander->addListener(L"Po_startingSingleTaskHander", [&] (std::any sender, SingleTaskEventArgs* e) {Po_startingSingleTaskHander(sender, e);});
		MetaMorpheusTask::FinishedSingleTaskHandler->addListener(L"Po_finishedSingleTaskHandler", [&] (std::any sender, SingleTaskEventArgs* e) {Po_finishedSingleTaskHandler(sender, e);});
		MetaMorpheusTask::FinishedWritingFileHandler->addListener(L"NewSuccessfullyFinishedWritingFile", [&] (std::any sender, SingleFileEventArgs* e) {NewSuccessfullyFinishedWritingFile(sender, e);});
		MetaMorpheusTask::StartingDataFileHandler->addListener(L"MyTaskEngine_StartingDataFileHandler", [&] (std::any sender, StringEventArgs* e) {MyTaskEngine_StartingDataFileHandler(sender, e);});
		MetaMorpheusTask::FinishedDataFileHandler->addListener(L"MyTaskEngine_FinishedDataFileHandler", [&] (std::any sender, StringEventArgs* e) {MyTaskEngine_FinishedDataFileHandler(sender, e);});
		MetaMorpheusTask::OutLabelStatusHandler->addListener(L"NewoutLabelStatus", [&] (std::any sender, StringEventArgs* e) {NewoutLabelStatus(sender, e);});
		MetaMorpheusTask::NewCollectionHandler->addListener(L"NewCollectionHandler", [&] (std::any sender, StringEventArgs* e) {NewCollectionHandler(sender, e);});
		MetaMorpheusTask::OutProgressHandler->addListener(L"NewoutProgressBar", [&] (std::any sender, ProgressEventArgs* e) {NewoutProgressBar(sender, e);});
		MetaMorpheusTask::WarnHandler->addListener(L"GuiWarnHandler", [&] (std::any sender, StringEventArgs* e) {GuiWarnHandler(sender, e);});

		MetaMorpheusEngine::OutProgressHandler->addListener(L"NewoutProgressBar", [&] (std::any sender, ProgressEventArgs* e) {NewoutProgressBar(sender, e);});
		MetaMorpheusEngine::OutLabelStatusHandler->addListener(L"NewoutLabelStatus", [&] (std::any sender, StringEventArgs* e) {NewoutLabelStatus(sender, e);});
		MetaMorpheusEngine::WarnHandler->addListener(L"GuiWarnHandler", [&] (std::any sender, StringEventArgs* e) {GuiWarnHandler(sender, e);});

		MyFileManager::WarnHandler->addListener(L"GuiWarnHandler", [&] (std::any sender, StringEventArgs* e) {GuiWarnHandler(sender, e);});

		UpdateSpectraFileGuiStuff();
		UpdateTaskGuiStuff();
		UpdateOutputFolderTextbox();
		FileSpecificParameters::ValidateFileSpecificVariableNames();
		SearchModifications::SetUpModSearchBoxes();

		// LOAD GUI SETTINGS

		if (FileSystem::fileExists(FileSystem::combine(GlobalVariables::getDataDir(), LR"(GUIsettings.toml)")))
		{
			GuiGlobalParams = Toml::ReadFile<GuiGlobalParams*>(FileSystem::combine(GlobalVariables::getDataDir(), LR"(GUIsettings.toml)"));
		}
		else
		{
			Toml::WriteFile(GuiGlobalParams, FileSystem::combine(GlobalVariables::getDataDir(), LR"(GUIsettings.toml)"), MetaMorpheusTask::tomlConfig);
			notificationsTextBox->Document = YoutubeWikiNotification();
		}

		if (GlobalVariables::getMetaMorpheusVersion().find(L"Not a release version") != std::wstring::npos)
		{
			GuiGlobalParams->setAskAboutUpdating(false);
		}

		try
		{
			GetVersionNumbersFromWeb();
		}
		catch (const std::runtime_error &e)
		{
			StringEventArgs tempVar(L"Could not get newest version from web: " + e.what(), nullptr);
			GuiWarnHandler(nullptr, &tempVar);
		}

		Application::Current->MainWindow.Closing += new CancelEventHandler(MainWindow_Closing);
	}

	FlowDocument *MainWindow::YoutubeWikiNotification()
	{

		FlowDocument *doc = notificationsTextBox::Document;
		Paragraph *p = new Paragraph();
		Run *run1 = new Run(L"Visit our ");
		Run *run2 = new Run(L"Wiki");
		Run *run3 = new Run(L" or ");
		Run *run4 = new Run(L"Youtube channel");
		Run *run5 = new Run(L" to check out what MetaMorpheus can do!" + L"\r\n");

		Hyperlink *wikiLink = new Hyperlink(run2);
		wikiLink->NavigateUri = new Uri(LR"(https://github.com/smith-chem-wisc/MetaMorpheus/wiki)");

		Hyperlink *youtubeLink = new Hyperlink(run4);
		youtubeLink->NavigateUri = new Uri(LR"(https://www.youtube.com/playlist?list=PLVk5tTSZ1aWlhNPh7jxPQ8pc0ElyzSUQb)");

		auto links = std::vector<Hyperlink*> {wikiLink, youtubeLink};

		p->Inlines->Add(run1);
		p->Inlines->Add(wikiLink);
		p->Inlines->Add(run3);
		p->Inlines->Add(youtubeLink);
		p->Inlines->Add(run5);

		for (auto link : links)
		{
			link->RequestNavigate += [&] (sender, e)
			{
				System::Diagnostics::Process::Start(e::Uri->ToString());
			};
		}

		doc->Blocks->Add(p);

//C# TO C++ CONVERTER TODO TASK: A 'delete youtubeLink' statement was not added since youtubeLink was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete wikiLink' statement was not added since wikiLink was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete run5' statement was not added since run5 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete run4' statement was not added since run4 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete run3' statement was not added since run3 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete run2' statement was not added since run2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete run1' statement was not added since run1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete p' statement was not added since p was passed to a method or constructor. Handle memory management manually.
		return doc;
	}

	std::wstring MainWindow::getNewestKnownVersion()
	{
		return privateNewestKnownVersion;
	}

	void MainWindow::setNewestKnownVersion(const std::wstring &value)
	{
		privateNewestKnownVersion = value;
	}

	void MainWindow::GetVersionNumbersFromWeb()
	{
		// Attempt to get current MetaMorpheus version
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var client = new HttpClient())
		{
			auto client = HttpClient();
			client.DefaultRequestHeaders->Add(L"User-Agent", L"Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; WOW64; Trident/6.0)");

//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var response = client.GetAsync("https://api.github.com/repos/smith-chem-wisc/MetaMorpheus/releases/latest").Result)
			{
				auto response = client.GetAsync(L"https://api.github.com/repos/smith-chem-wisc/MetaMorpheus/releases/latest").Result;
				auto json = response.Content::ReadAsStringAsync().Result;
				JObject *deserialized = JObject::Parse(json);
				auto assets = deserialized[L"assets"]->Select([&] (std::any b)
				{
					b[L"name"].ToString();
				}).ToList();
				if (!std::find(assets.begin(), assets.end(), L"MetaMorpheusInstaller.msi") != assets.end())
				{
					throw MetaMorpheusException(std::wstring(L"A new version of MetaMorpheus was detected, but the files haven't been") + L" uploaded yet. Try again in a few minutes.");
				}
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				setNewestKnownVersion(deserialized[L"tag_name"].ToString());
			}
		}
	}

	void MainWindow::MyWindow_Loaded(std::any sender, RoutedEventArgs *e)
	{
		if (getNewestKnownVersion() != L"" && GlobalVariables::getMetaMorpheusVersion() != getNewestKnownVersion() && GuiGlobalParams->getAskAboutUpdating())
		{
			try
			{
				MetaUpdater *newwind = new MetaUpdater();
				newwind->ShowDialog();

				delete newwind;
			}
			catch (const std::runtime_error &ex)
			{
				MessageBox::Show(ex.what());
			}
		}

		this->KeyDown += new KeyEventHandler(Window_KeyDown);

		// hide the "InProgress" column
		dataGridProteinDatabases::Columns::Where([&] (std::any p)
		{
			p::Header->Equals(L"InProgress");
		}).First()->Visibility = Visibility::Hidden;
		dataGridSpectraFiles::Columns::Where([&] (std::any p)
		{
			p::Header->Equals(L"InProgress");
		}).First()->Visibility = Visibility::Hidden;

		PrintErrorsReadingMods();
	}

	void MainWindow::EverythingRunnerEngine_FinishedWritingAllResultsFileHandler(std::any sender, StringEventArgs *e)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				EverythingRunnerEngine_FinishedWritingAllResultsFileHandler(sender, e);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			InRunTask tempVar2(L"All Task Results", nullptr);
			DynamicTasksObservableCollection->Add(&tempVar2);
			DynamicTasksObservableCollection->Last()->Progress = 100;
			OutputFileForTreeView tempVar3(e->getS(), Path::GetFileNameWithoutExtension(e->getS()));
			DynamicTasksObservableCollection->Last().Children->Add(&tempVar3);
		}
	}

	void MainWindow::GuiWarnHandler(std::any sender, StringEventArgs *e)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				GuiWarnHandler(sender, e);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			notificationsTextBox::AppendText(e->getS());
			notificationsTextBox::AppendText(L"\r\n");
		}
	}

	void MainWindow::MyTaskEngine_FinishedDataFileHandler(std::any sender, StringEventArgs *s)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				MyTaskEngine_FinishedDataFileHandler(sender, s);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			auto huh = SpectraFilesObservableCollection->First([&] (std::any b)
			{
				b::FilePath->Equals(s->getS());
			});
			huh->SetInProgress(false);

			dataGridSpectraFiles::Items->Refresh();
		}
	}

	void MainWindow::MyTaskEngine_StartingDataFileHandler(std::any sender, StringEventArgs *s)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				MyTaskEngine_StartingDataFileHandler(sender, s);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			auto huh = SpectraFilesObservableCollection->First([&] (std::any b)
			{
				b::FilePath->Equals(s->getS());
			});
			huh->SetInProgress(true);
			dataGridSpectraFiles::Items->Refresh();
		}
	}

	void MainWindow::AddNewDB(std::any sender, XmlForTaskListEventArgs *e)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				AddNewDB(sender, e);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			for (auto uu : ProteinDbObservableCollection)
			{
				uu->setUse(false);
			}
			for (auto uu : e->NewDatabases)
			{
				ProteinDbForDataGrid tempVar2(uu);
				ProteinDbObservableCollection->Add(&tempVar2);
			}
			dataGridProteinDatabases::Items->Refresh();
		}
	}

	void MainWindow::AddNewSpectra(std::any sender, StringListEventArgs *e)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				AddNewSpectra(sender, e);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			auto newFiles = e->getStringList().ToList();
			for (auto oldFile : SpectraFilesObservableCollection)
			{
				if (!std::find(newFiles.begin(), newFiles.end(), oldFile->getFilePath()) != newFiles.end())
				{
					oldFile->setUse(false);
				}
			}

			auto files = SpectraFilesObservableCollection->Select([&] (std::any p)
			{
				p::FilePath;
			}).ToList();
			for (auto newRawData : newFiles.Where([&] (std::any p)
			{
				!std::find(files.begin(), files.end(), p) != files.end();
			}))
			{
				RawDataForDataGrid tempVar2(newRawData);
				SpectraFilesObservableCollection->Add(&tempVar2);
			}

			UpdateOutputFolderTextbox();
		}
	}

	void MainWindow::AddNewFileSpecificToml(std::any sender, StringListEventArgs *e)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				AddNewFileSpecificToml(sender, e);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			for (auto path : e->getStringList())
			{
				UpdateFileSpecificParamsDisplayJustAdded(path);
			}
		}
	}

	void MainWindow::UpdateOutputFolderTextbox()
	{
		if (SpectraFilesObservableCollection->Any())
		{
			// if current output folder is blank and there is a spectra file, use the spectra file's path as the output path
			if (StringHelper::isEmptyOrWhiteSpace(OutputFolderTextBox->Text))
			{
				auto pathOfFirstSpectraFile = FileSystem::getDirectoryName(SpectraFilesObservableCollection->First().FilePath);
				OutputFolderTextBox->Text = FileSystem::combine(pathOfFirstSpectraFile, LR"($DATETIME)");
			}
			// else do nothing (do not override if there is a path already there; might clear user-defined path)
		}
		else
		{
			// no spectra files; clear the output folder from the GUI
			OutputFolderTextBox->Clear();
		}
	}

	void MainWindow::Po_startingSingleTaskHander(std::any sender, SingleTaskEventArgs *s)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				Po_startingSingleTaskHander(sender, s);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			auto theTask = DynamicTasksObservableCollection->First([&] (std::any b)
			{
				b::DisplayName->Equals(s->getDisplayName());
			});
			theTask->Status = L"Starting...";

			dataGridSpectraFiles::Items->Refresh();
			dataGridProteinDatabases::Items->Refresh();
		}
	}

	void MainWindow::Po_finishedSingleTaskHandler(std::any sender, SingleTaskEventArgs *s)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				Po_finishedSingleTaskHandler(sender, s);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			auto theTask = DynamicTasksObservableCollection->First([&] (std::any b)
			{
				b::DisplayName->Equals(s->getDisplayName());
			});
			theTask->IsIndeterminate = false;
			theTask->Progress = 100;
			theTask->Status = L"Done!";

			dataGridSpectraFiles::Items->Refresh();
			dataGridProteinDatabases::Items->Refresh();
		}
	}

	void MainWindow::ClearSpectraFiles_Click(std::any sender, RoutedEventArgs *e)
	{
		SpectraFilesObservableCollection->Clear();
		UpdateOutputFolderTextbox();
	}

	void MainWindow::OpenOutputFolder_Click(std::any sender, RoutedEventArgs *e)
	{
		std::wstring outputFolder = OutputFolderTextBox->Text;
		if (outputFolder.find(L"$DATETIME") != std::wstring::npos)
		{
			// the exact file path isn't known, so just open the parent directory
			outputFolder = Directory::GetParent(outputFolder)->FullName;
		}

		if (!FileSystem::directoryExists(outputFolder) && !outputFolder.empty())
		{
			// create the directory if it doesn't exist yet
			try
			{
				FileSystem::createDirectory(outputFolder);
			}
			catch (const std::runtime_error &ex)
			{
				StringEventArgs tempVar(L"Error opening directory: " + ex.what(), nullptr);
				GuiWarnHandler(nullptr, &tempVar);
			}
		}

		if (FileSystem::directoryExists(outputFolder))
		{
			// open the directory
			System::Diagnostics::ProcessStartInfo *tempVar2 = new System::Diagnostics::ProcessStartInfo();
			tempVar2->FileName = outputFolder;
			tempVar2->UseShellExecute = true;
			tempVar2->Verb = L"open";
			System::Diagnostics::Process::Start(tempVar2);
		}
		else
		{
			// this should only happen if the file path is empty or something unexpected happened
			StringEventArgs tempVar3(L"Output folder does not exist", nullptr);
			GuiWarnHandler(nullptr, &tempVar3);
		}
	}

	void MainWindow::AddProteinDatabase_Click(std::any sender, RoutedEventArgs *e)
	{
		Microsoft::Win32::OpenFileDialog *openPicker = new Microsoft::Win32::OpenFileDialog();
		openPicker->Filter = L"Database Files|*.xml;*.xml.gz;*.fasta;*.fa";
		openPicker->FilterIndex = 1;
		openPicker->RestoreDirectory = true;
		openPicker->Multiselect = true;
		if (openPicker->ShowDialog() == true)
		{
			for (auto filepath : openPicker->FileNames.OrderBy([&] (std::any p)
			{
			delete openPicker;
				return p;
			}))
			{
				AddAFile(filepath);
			}
		}
		dataGridProteinDatabases::Items->Refresh();

		delete openPicker;
	}

	void MainWindow::AddSpectraFile_Click(std::any sender, RoutedEventArgs *e)
	{
		Microsoft::Win32::OpenFileDialog *openFileDialog1 = new Microsoft::Win32::OpenFileDialog();
		openFileDialog1->Filter = L"Spectra Files(*.raw;*.mzML;*.mgf)|*.raw;*.mzML;*.mgf";
		openFileDialog1->FilterIndex = 1;
		openFileDialog1->RestoreDirectory = true;
		openFileDialog1->Multiselect = true;
		if (openFileDialog1->ShowDialog() == true)
		{
			for (auto rawDataFromSelected : openFileDialog1->FileNames.OrderBy([&] (std::any p)
			{
			delete openFileDialog1;
				return p;
			}))
			{
				AddAFile(rawDataFromSelected);
			}
		}
		dataGridSpectraFiles::Items->Refresh();

		delete openFileDialog1;
	}

	void MainWindow::Window_Drop(std::any sender, DragEventArgs *e)
	{
		if (LoadTaskButton::IsEnabled)
		{
			std::vector<std::wstring> files = (static_cast<std::vector<std::wstring>>(e->Data->GetData(DataFormats::FileDrop)))->OrderBy([&] (std::any p)
			{
				return p;
			})->ToArray();

			if (files.size() > 0)
			{
				for (auto draggedFilePath : files)
				{
					if (FileSystem::directoryExists(draggedFilePath))
					{
						for (auto file : Directory::EnumerateFiles(draggedFilePath, L"*.*", SearchOption::AllDirectories))
						{
							AddAFile(file);
						}
					}
					else
					{
						AddAFile(draggedFilePath);
					}
					dataGridSpectraFiles::CommitEdit(DataGridEditingUnit::Row, true);
					dataGridProteinDatabases::CommitEdit(DataGridEditingUnit::Row, true);
					dataGridSpectraFiles::Items->Refresh();
					dataGridProteinDatabases::Items->Refresh();
				}
			}
			UpdateTaskGuiStuff();
		}
	}

	void MainWindow::AddAFile(const std::wstring &draggedFilePath)
	{
		// this line is NOT used because .xml.gz (extensions with two dots) mess up with Path.GetExtension
		//var theExtension = Path.GetExtension(draggedFilePath).ToLowerInvariant();

		// we need to get the filename before parsing out the extension because if we assume that everything after the dot
		// is the extension and there are dots in the file path (i.e. in a folder name), this will mess up
		auto filename = FileSystem::getFileName(draggedFilePath);
//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
		auto theExtension = Path::GetExtension(filename).ToLowerInvariant();
		bool compressed = StringHelper::endsWith(theExtension, L"gz"); // allows for .bgz and .tgz, too which are used on occasion
//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
		theExtension = compressed ? Path::GetExtension(Path::GetFileNameWithoutExtension(filename)).ToLowerInvariant() : theExtension;

//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//		switch (theExtension)
//ORIGINAL LINE: case ".raw":
		if (theExtension == L".raw")
		{
				if (!WarnedAboutThermoAlready)
				{
					// check for MSFileReader and display a warning if the expected DLLs are not found
					auto versionCheckerResult = MyFileManager::ValidateThermoMsFileReaderVersion();

					if (versionCheckerResult->Equals(MyFileManager::ThermoMsFileReaderVersionCheck::IncorrectVersion))
					{
						StringEventArgs tempVar(L"Warning! Thermo MSFileReader is not version 3.0 SP2; a crash may result from searching this .raw file", nullptr);
						GuiWarnHandler(nullptr, &tempVar);
					}
					else if (versionCheckerResult->Equals(MyFileManager::ThermoMsFileReaderVersionCheck::DllsNotFound))
					{
						StringEventArgs tempVar2(L"Warning! Cannot find Thermo MSFileReader (v3.0 SP2 is preferred); a crash may result from searching this .raw file", nullptr);
						GuiWarnHandler(nullptr, &tempVar2);
					}
					else if (versionCheckerResult->Equals(MyFileManager::ThermoMsFileReaderVersionCheck::SomeDllsMissing))
					{
						StringEventArgs tempVar3(L"Warning! Found only some of the expected Thermo MSFileReader .dll files; a crash may result from searching this .raw file", nullptr);
						GuiWarnHandler(nullptr, &tempVar3);
					}

					// check for ManagedThermoHelperLayer.dll and display a warning if it's not found
					// this is one hacky way of checking if the user has C++ redistributable installed
					std::wstring assumedManagedThermoHelperLayerDllPath = FileSystem::combine(Environment::CurrentDirectory, L"ManagedThermoHelperLayer.dll");
					if (!FileSystem::fileExists(assumedManagedThermoHelperLayerDllPath))
					{
						StringEventArgs tempVar4(std::wstring(L"Warning! Cannot find Microsoft Visual C++ Redistributable; ") + L"a crash may result from searching this .raw file. If you have just installed the C++ redistributable, " + L"please uninstall and reinstall MetaMorpheus", nullptr);
						GuiWarnHandler(nullptr, &tempVar4);
					}
				}

				WarnedAboutThermoAlready = true;
				goto case L".mzml";

		}
//ORIGINAL LINE: case ".mgf":
		else if (theExtension == L".mgf")
		{
				StringEventArgs tempVar5(L".mgf files lack MS1 spectra, which are needed for quantification and searching for coisolated peptides. All other features of MetaMorpheus will function.", nullptr);
				GuiWarnHandler(nullptr, &tempVar5);
				goto case L".mzml";

		}
//ORIGINAL LINE: case ".mzml":
		else if (theExtension == L".mzml")
		{
				if (compressed) // not implemented yet
				{
					StringEventArgs tempVar6(L"Cannot read, try uncompressing: " + draggedFilePath, nullptr);
					GuiWarnHandler(nullptr, &tempVar6);
					break;
				}
				RawDataForDataGrid *zz = new RawDataForDataGrid(draggedFilePath);
				if (!SpectraFileExists(SpectraFilesObservableCollection, zz))
				{
					SpectraFilesObservableCollection->Add(zz);
				}
				UpdateFileSpecificParamsDisplayJustAdded(Path::ChangeExtension(draggedFilePath, L".toml"));
				UpdateOutputFolderTextbox();

		}
//ORIGINAL LINE: case ".xml":
		else if (theExtension == L".xml" || theExtension == L".fasta" || theExtension == L".fa")
		{
				ProteinDbForDataGrid *uu = new ProteinDbForDataGrid(draggedFilePath);
				if (!DatabaseExists(ProteinDbObservableCollection, uu))
				{
					ProteinDbObservableCollection->Add(uu);
					if (theExtension == L".xml")
					{
						try
						{
							GlobalVariables::AddMods(UsefulProteomicsDatabases::ProteinDbLoader::GetPtmListFromProteinXml(draggedFilePath).OfType<Modification*>(), true);

							PrintErrorsReadingMods();
						}
						catch (const std::runtime_error &ee)
						{
							MessageBox::Show(ee.what());
							StringEventArgs tempVar7(L"Cannot parse modification info from: " + draggedFilePath, nullptr);
							GuiWarnHandler(nullptr, &tempVar7);
							ProteinDbObservableCollection->Remove(uu);
						}
					}
				}

		}
//ORIGINAL LINE: case ".toml":
		else if (theExtension == L".toml")
		{
				TomlTable *tomlFile = nullptr;
				try
				{
					tomlFile = Toml::ReadFile(draggedFilePath, MetaMorpheusTask::tomlConfig);
				}
				catch (const std::runtime_error &e1)
				{
					StringEventArgs tempVar8(L"Cannot read toml: " + draggedFilePath, nullptr);
					GuiWarnHandler(nullptr, &tempVar8);
					break;
				}

				if (tomlFile->ContainsKey(L"TaskType"))
				{
					try
					{
//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//						switch (tomlFile.Get<string>("TaskType"))
//ORIGINAL LINE: case "Search":
						if (tomlFile->Get<std::wstring>(L"TaskType") == L"Search")
						{
								auto ye1 = Toml::ReadFile<SearchTask*>(draggedFilePath, MetaMorpheusTask::tomlConfig);
								AddTaskToCollection(ye1);

						}
//ORIGINAL LINE: case "Calibrate":
						else if (tomlFile->Get<std::wstring>(L"TaskType") == L"Calibrate")
						{
								auto ye2 = Toml::ReadFile<CalibrationTask*>(draggedFilePath, MetaMorpheusTask::tomlConfig);
								AddTaskToCollection(ye2);

						}
//ORIGINAL LINE: case "Gptmd":
						else if (tomlFile->Get<std::wstring>(L"TaskType") == L"Gptmd")
						{
								auto ye3 = Toml::ReadFile<GptmdTask*>(draggedFilePath, MetaMorpheusTask::tomlConfig);
								AddTaskToCollection(ye3);

						}
//ORIGINAL LINE: case "XLSearch":
						else if (tomlFile->Get<std::wstring>(L"TaskType") == L"XLSearch")
						{
								auto ye4 = Toml::ReadFile<XLSearchTask*>(draggedFilePath, MetaMorpheusTask::tomlConfig);
								AddTaskToCollection(ye4);
						}
					}
					catch (const std::runtime_error &e)
					{
						StringEventArgs tempVar9(L"Cannot read task toml: " + e.what(), nullptr);
						GuiWarnHandler(nullptr, &tempVar9);
					}
				}

		}
		else
		{
				StringEventArgs tempVar10(L"Unrecognized file type: " + theExtension, nullptr);
				GuiWarnHandler(nullptr, &tempVar10);

//C# TO C++ CONVERTER TODO TASK: A 'delete uu' statement was not added since uu was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete zz' statement was not added since zz was passed to a method or constructor. Handle memory management manually.
		}
	}

	void MainWindow::AddTaskToCollection(MetaMorpheusTask *ye)
	{
		PreRunTask *te = new PreRunTask(ye);
		StaticTasksObservableCollection->Add(te);
		StaticTasksObservableCollection->Last()->DisplayName = L"Task" + std::to_wstring(StaticTasksObservableCollection->find(te) + 1) + L"-" + ye->getCommonParameters()->getTaskDescriptor();

//C# TO C++ CONVERTER TODO TASK: A 'delete te' statement was not added since te was passed to a method or constructor. Handle memory management manually.
	}

	void MainWindow::Row_DoubleClick(std::any sender, MouseButtonEventArgs *e)
	{
		auto ye = dynamic_cast<DataGridCell*>(sender);

		// prevent opening protein DB or spectra files if a run is in progress
		if ((dynamic_cast<ProteinDbForDataGrid*>(ye->DataContext) != nullptr || dynamic_cast<RawDataForDataGrid*>(ye->DataContext) != nullptr) && !LoadTaskButton::IsEnabled)
		{
			return;
		}

		// open the file with the default process for this file format
//C# TO C++ CONVERTER TODO TASK: C++ has no equivalent to C# pattern variables in 'is' expressions:
//ORIGINAL LINE: if (ye.Content is TextBlock hm && hm != null && !string.IsNullOrEmpty(hm.Text))
		if (dynamic_cast<TextBlock*>(ye->Content) != nullptr hm && hm != nullptr && !hm->Text->empty())
		{
			try
			{
				System::Diagnostics::Process::Start(hm->Text);
			}
			catch (const std::runtime_error &e1)
			{
			}
		}
	}

	void MainWindow::RunAllTasks_Click(std::any sender, RoutedEventArgs *e)
	{
		GlobalVariables::setStopLoops(false);
		CancelButton->IsEnabled = true;

		// check for valid tasks/spectra files/protein databases
		if (!StaticTasksObservableCollection->Any())
		{
			StringEventArgs tempVar(L"You need to add at least one task!", nullptr);
			GuiWarnHandler(nullptr, &tempVar);
			return;
		}
		if (!SpectraFilesObservableCollection->Any())
		{
			StringEventArgs tempVar2(L"You need to add at least one spectra file!", nullptr);
			GuiWarnHandler(nullptr, &tempVar2);
			return;
		}
		if (!ProteinDbObservableCollection->Any())
		{
			StringEventArgs tempVar3(L"You need to add at least one protein database!", nullptr);
			GuiWarnHandler(nullptr, &tempVar3);
			return;
		}

		DynamicTasksObservableCollection = new ObservableCollection<InRunTask*>();

		for (int i = 0; i < StaticTasksObservableCollection->Count; i++)
		{
			InRunTask tempVar4(L"Task" + std::to_wstring(i + 1) + L"-" + StaticTasksObservableCollection[i]->metaMorpheusTask.getCommonParameters().getTaskDescriptor(), StaticTasksObservableCollection[i]->metaMorpheusTask);
			DynamicTasksObservableCollection->Add(&tempVar4);
		}
		tasksTreeView->DataContext = DynamicTasksObservableCollection;

		notificationsTextBox::Document::Blocks->Clear();

		// output folder
		if (OutputFolderTextBox->Text->empty())
		{
			auto pathOfFirstSpectraFile = FileSystem::getDirectoryName(SpectraFilesObservableCollection->First().FilePath);
			OutputFolderTextBox->Text = FileSystem::combine(pathOfFirstSpectraFile, LR"($DATETIME)");
		}

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		auto startTimeForAllFilenames = DateTime::Now.ToString(L"yyyy-MM-dd-HH-mm-ss", CultureInfo::InvariantCulture);
		std::wstring outputFolder = OutputFolderTextBox->Text->Replace(L"$DATETIME", startTimeForAllFilenames);
		OutputFolderTextBox->Text = outputFolder;

		// check that experimental design is defined if normalization is enabled
		// TODO: move all of this over to EverythingRunnerEngine
		auto searchTasks = StaticTasksObservableCollection->Where([&] (std::any p)
		{
			return p::metaMorpheusTask->TaskType == MyTask::Search;
		})->Select([&] (std::any p)
		{
			(SearchTask)p.metaMorpheusTask;
		});

		std::wstring pathToExperDesign = Directory::GetParent(SpectraFilesObservableCollection->First().FilePath)->FullName;
		pathToExperDesign = FileSystem::combine(pathToExperDesign, GlobalVariables::getExperimentalDesignFileName());

		for (auto searchTask : searchTasks->Where([&] (std::any p)
		{
			p::SearchParameters::Normalize;
		}))
		{
			if (!FileSystem::fileExists(pathToExperDesign))
			{
				MessageBox::Show(std::wstring(L"Experimental design must be defined for normalization!\n") + L"Click the \"Experimental Design\" button in the bottom left by the spectra files");
				return;
			}

			// check that experimental design is OK (spectra files may have been added after exper design was defined)
			// TODO: experimental design might still have flaws if user edited the file manually, need to check for this
			auto experDesign = File::ReadAllLines(pathToExperDesign).ToDictionary([&] (std::any p)
			{
				p->Split(L'\t')[0];
			}, [&] (std::any p)
			{
				return p;
			});
			auto filesToUse = std::unordered_set<std::wstring>(SpectraFilesObservableCollection->Select([&] (std::any p)
			{
				Path::GetFileNameWithoutExtension(p::FileName);
			}));
			auto experDesignFilesDefined = std::unordered_set<std::wstring>(experDesign.Keys);

			auto undefined = filesToUse.Except(experDesignFilesDefined);

			if (undefined->Any())
			{
				MessageBox::Show(L"Need to define experimental design parameters for file: " + undefined->First());
				return;
			}
		}
		BtnQuantSet->IsEnabled = false;

		// everything is OK to run
		EverythingRunnerEngine *a = new EverythingRunnerEngine(DynamicTasksObservableCollection->Select([&] (std::any b)
		{
			(b::DisplayName, b::Task);
		}).ToList(), SpectraFilesObservableCollection->Where([&] (std::any b)
		{
			b::Use;
		})->Select([&] (std::any b)
		{
			b::FilePath;
		}).ToList(), ProteinDbObservableCollection->Where([&] (std::any b)
		{
			b::Use;
		})->Select([&] (std::any b)
		{
			new DbForTask(b::FilePath, b::Contaminant);
		}).ToList(), outputFolder);

		auto t = new Task(a->Run);
		t->ContinueWith([&] (System::Threading::Tasks::Task* obj) {EverythingRunnerExceptionHandler(obj);}, TaskContinuationOptions::OnlyOnFaulted);
		t->Start();

		delete t;
		delete a;
	}

	void MainWindow::EverythingRunnerExceptionHandler(Task *obj)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				EverythingRunnerExceptionHandler(obj);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			std::runtime_error e = obj->Exception;
			while (e.InnerException != nullptr)
			{
				e = e.InnerException;
			}
			auto message = L"Run failed, Exception: " + e.what();
			auto messageBoxResult = System::Windows::MessageBox::Show(message + L"\n\nWould you like to report this crash?", L"Runtime Error", MessageBoxButton::YesNo);
			notificationsTextBox::AppendText(message + L"\r\n");
			std::runtime_error exception = e;
			//Find Output Folder
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			std::wstring outputFolder = e.Data[L"folder"].ToString();
			std::wstring tomlText = L"";
			if (FileSystem::directoryExists(outputFolder))
			{
				auto tomls = Directory::GetFiles(outputFolder, L"*.toml");
				//will only be 1 toml per task
				for (auto tomlFile : tomls)
				{
					tomlText += L"\n" + File::ReadAllText(tomlFile);
				}
				if (!tomls.Any())
				{
					tomlText = L"TOML not found";
				}
			}
			else
			{
				tomlText = L"Directory not found";
			}
			if (messageBoxResult == MessageBoxResult::Yes)
			{
				std::wstring body = exception.what() + L"%0D%0A" + exception.Data +
				   L"%0D%0A" + exception.StackTrace +
				   L"%0D%0A" + exception.Source +
				   L"%0D%0A %0D%0A %0D%0A %0D%0A SYSTEM INFO: %0D%0A " +
					SystemInfo::CompleteSystemInfo() +
				   L"%0D%0A%0D%0A MetaMorpheus: version " + GlobalVariables::getMetaMorpheusVersion()
				   + L"%0D%0A %0D%0A %0D%0A %0D%0A TOML: %0D%0A " + tomlText;
				body = StringHelper::replace(body, L'&', L' ');
				body = StringHelper::replace(body, L"\n", L"%0D%0A");
				body = StringHelper::replace(body, L"\r", L"%0D%0A");
				std::wstring mailto = StringHelper::formatSimple(L"mailto:{0}?Subject=MetaMorpheus. Issue:&Body={1}", L"mm_support@chem.wisc.edu", body);
				System::Diagnostics::Process::Start(mailto);
				std::wcout << body << std::endl;
			}
			ResetTasksButton->IsEnabled = true;
		}
	}

	void MainWindow::ClearTasks_Click(std::any sender, RoutedEventArgs *e)
	{
		StaticTasksObservableCollection->Clear();
		UpdateTaskGuiStuff();
	}

	void MainWindow::UpdateTaskGuiStuff()
	{
		if (StaticTasksObservableCollection->Count == 0)
		{
			RunTasksButton->IsEnabled = false;
			DeleteSelectedTaskButton->IsEnabled = false;
			ClearTasksButton->IsEnabled = false;
			ResetTasksButton->IsEnabled = false;
		}
		else
		{
			RunTasksButton->IsEnabled = true;
			DeleteSelectedTaskButton->IsEnabled = true;
			ClearTasksButton->IsEnabled = true;

			// this exists so that when a task is deleted, the remaining tasks are renamed to keep the task numbers correct
			for (int i = 0; i < StaticTasksObservableCollection->Count; i++)
			{
				std::wstring newName = L"Task" + std::to_wstring(i + 1) + L"-" + StaticTasksObservableCollection[i]->metaMorpheusTask.getCommonParameters().getTaskDescriptor();
				StaticTasksObservableCollection[i]->setDisplayName(newName);
			}
			tasksTreeView::Items->Refresh();
		}
	}

	void MainWindow::UpdateSpectraFileGuiStuff()
	{
		ChangeFileParameters->IsEnabled = SelectedRawFiles->Count > 0 && LoadTaskButton::IsEnabled;
	}

	void MainWindow::AddSearchTaskButton_Click(std::any sender, RoutedEventArgs *e)
	{
		auto dialog = new SearchTaskWindow();
		if (dialog->ShowDialog() == true)
		{
			AddTaskToCollection(dialog->getTheTask());
			UpdateTaskGuiStuff();
		}

		delete dialog;
	}

	void MainWindow::AddCalibrateTaskButton_Click(std::any sender, RoutedEventArgs *e)
	{
		auto dialog = new CalibrateTaskWindow();
		if (dialog->ShowDialog() == true)
		{
			AddTaskToCollection(dialog->getTheTask());
			UpdateTaskGuiStuff();
		}

		delete dialog;
	}

	void MainWindow::AddGPTMDTaskButton_Click(std::any sender, RoutedEventArgs *e)
	{
		auto dialog = new GptmdTaskWindow();
		if (dialog->ShowDialog() == true)
		{
			AddTaskToCollection(dialog->getTheTask());
			UpdateTaskGuiStuff();
		}

		delete dialog;
	}

	void MainWindow::BtnAddCrosslinkSearch_Click(std::any sender, RoutedEventArgs *e)
	{
		auto dialog = new XLSearchTaskWindow();
		if (dialog->ShowDialog() == true)
		{
			AddTaskToCollection(dialog->getTheTask());
			UpdateTaskGuiStuff();
		}

		delete dialog;
	}

	void MainWindow::DeleteSelectedTask(std::any sender, RoutedEventArgs *e)
	{
		auto selectedTask = static_cast<PreRunTask*>(tasksTreeView::SelectedItem);
		if (selectedTask != nullptr)
		{
			StaticTasksObservableCollection->Remove(selectedTask);
			UpdateTaskGuiStuff();
		}
	}

	void MainWindow::MoveSelectedTask(std::any sender, RoutedEventArgs *e, bool moveTaskUp)
	{
		auto selectedTask = static_cast<PreRunTask*>(tasksTreeView::SelectedItem);
		if (selectedTask == nullptr)
		{
			return;
		}

		int indexOfSelected = StaticTasksObservableCollection->find(selectedTask);
		int indexToMoveTo = indexOfSelected - 1;
		if (moveTaskUp)
		{
			indexToMoveTo = indexOfSelected + 1;
		}

		if (indexToMoveTo >= 0 && indexToMoveTo < StaticTasksObservableCollection->Count)
		{
			auto temp = StaticTasksObservableCollection[indexToMoveTo];
			StaticTasksObservableCollection[indexToMoveTo] = selectedTask;
			StaticTasksObservableCollection[indexOfSelected] = temp;

			UpdateTaskGuiStuff();

			auto item = tasksTreeView::ItemContainerGenerator::ContainerFromItem(selectedTask);
			(static_cast<TreeViewItem*>(item))->IsSelected = true;
		}
	}

	void MainWindow::Window_KeyDown(std::any sender, KeyEventArgs *e)
	{
		if (LoadTaskButton::IsEnabled)
		{
			// delete selected task
			if (e->Key == Key->Delete || e->Key == Key->Back)
			{
				DeleteSelectedTask(sender, e);
				e->Handled = true;
			}

			// move task up
			if (e->Key == Key->Add || e->Key == Key->OemPlus)
			{
				MoveSelectedTask(sender, e, true);
				e->Handled = true;
			}

			// move task down
			if (e->Key == Key->Subtract || e->Key == Key->OemMinus)
			{
				MoveSelectedTask(sender, e, false);
				e->Handled = true;
			}
		}
	}

	void MainWindow::NewCollectionHandler(std::any sender, StringEventArgs *s)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				NewCollectionHandler(sender, s);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			// Find the task or the collection!!!

			ForTreeView *theEntityOnWhichToUpdateLabel = DynamicTasksObservableCollection->First([&] (std::any b)
			{
				b::Id->Equals(s->NestedIDs[0]);
			});

			for (int i = 1; i < s->NestedIDs.size() - 1; i++)
			{
				auto hm = s->NestedIDs[i];
				try
				{
					theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel->getChildren()->First([&] (std::any b)
					{
						b::Id->Equals(hm);
					});
				}
				catch (...)
				{
					CollectionForTreeView tempVar2(hm, hm);
					theEntityOnWhichToUpdateLabel->getChildren()->Add(&tempVar2);
					theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel->getChildren()->First([&] (std::any b)
					{
						b::Id->Equals(hm);
					});
				}
			}

			CollectionForTreeView tempVar3(s->getS(), s->NestedIDs.back());
			theEntityOnWhichToUpdateLabel->getChildren()->Add(&tempVar3);
		}
	}

	void MainWindow::NewoutLabelStatus(std::any sender, StringEventArgs *s)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				NewoutLabelStatus(sender, s);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			// Find the task or the collection!!!

			ForTreeView *theEntityOnWhichToUpdateLabel = DynamicTasksObservableCollection->First([&] (std::any b)
			{
				b::Id->Equals(s->NestedIDs[0]);
			});

			for (auto hm : s->NestedIDs.Skip(1))
			{
				try
				{
					theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel->getChildren()->First([&] (std::any b)
					{
						b::Id->Equals(hm);
					});
				}
				catch (...)
				{
					CollectionForTreeView tempVar2(hm, hm);
					theEntityOnWhichToUpdateLabel->getChildren()->Add(&tempVar2);
					theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel->getChildren()->First([&] (std::any b)
					{
						b::Id->Equals(hm);
					});
				}
			}

			theEntityOnWhichToUpdateLabel->setStatus(s->getS());
			theEntityOnWhichToUpdateLabel->setIsIndeterminate(true);
		}
	}

	void MainWindow::NewoutProgressBar(std::any sender, ProgressEventArgs *s)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				NewoutProgressBar(sender, s);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			ForTreeView *theEntityOnWhichToUpdateLabel = DynamicTasksObservableCollection->First([&] (std::any b)
			{
				b::Id->Equals(s->NestedIDs[0]);
			});

			for (auto hm : s->NestedIDs.Skip(1))
			{
				try
				{
					theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel->getChildren()->First([&] (std::any b)
					{
						b::Id->Equals(hm);
					});
				}
				catch (...)
				{
					CollectionForTreeView tempVar2(hm, hm);
					theEntityOnWhichToUpdateLabel->getChildren()->Add(&tempVar2);
					theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel->getChildren()->First([&] (std::any b)
					{
						b::Id->Equals(hm);
					});
				}
			}

			theEntityOnWhichToUpdateLabel->setStatus(s->V);
			theEntityOnWhichToUpdateLabel->setIsIndeterminate(false);
			theEntityOnWhichToUpdateLabel->setProgress(s->NewProgress);
		}
	}

	void MainWindow::NewRefreshBetweenTasks(std::any sender, EventArgs *e)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				NewRefreshBetweenTasks(sender, e);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			dataGridSpectraFiles::Items->Refresh();
			dataGridProteinDatabases::Items->Refresh();
		}
	}

	void MainWindow::NewSuccessfullyStartingAllTasks(std::any sender, EventArgs *e)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				NewSuccessfullyStartingAllTasks(sender, e);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			dataGridSpectraFiles::Items->Refresh();

			ClearTasksButton->IsEnabled = false;
			DeleteSelectedTaskButton->IsEnabled = false;
			RunTasksButton->IsEnabled = false;
			LoadTaskButton->IsEnabled = false;

			addCalibrateTaskButton->IsEnabled = false;
			addGPTMDTaskButton->IsEnabled = false;
			addSearchTaskButton->IsEnabled = false;
			btnAddCrosslinkSearch->IsEnabled = false;

			AddXML->IsEnabled = false;
			ClearXML->IsEnabled = false;
			AddRaw->IsEnabled = false;
			ClearRaw->IsEnabled = false;

			OutputFolderTextBox->IsEnabled = false;

			dataGridSpectraFiles->IsReadOnly = true;
			dataGridProteinDatabases->IsReadOnly = true;
		}
	}

	void MainWindow::NewSuccessfullyFinishedAllTasks(std::any sender, StringEventArgs *e)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				NewSuccessfullyFinishedAllTasks(sender, e);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			ResetTasksButton->IsEnabled = true;

			dataGridSpectraFiles::Items->Refresh();
		}
	}

	void MainWindow::NewSuccessfullyFinishedWritingFile(std::any sender, SingleFileEventArgs *v)
	{
		if (!Dispatcher::CheckAccess())
		{
			std::function<void()> tempVar([&] ()
			{
				NewSuccessfullyFinishedWritingFile(sender, v);
			});
			Dispatcher::BeginInvoke(&tempVar);
		}
		else
		{
			ForTreeView *AddWrittenFileToThisOne = DynamicTasksObservableCollection->First([&] (std::any b)
			{
				b::Id->Equals(v->NestedIDs[0]);
			});

			for (auto hm : v->NestedIDs.Skip(1))
			{
				try
				{
					AddWrittenFileToThisOne = AddWrittenFileToThisOne->getChildren()->First([&] (std::any b)
					{
						b::Id->Equals(hm);
					});
				}
				catch (...)
				{
				}
			}
			OutputFileForTreeView tempVar2(v->getWrittenFile(), FileSystem::getFileName(v->getWrittenFile()));
			AddWrittenFileToThisOne->getChildren()->Add(&tempVar2);
		}
	}

	void MainWindow::ClearXML_Click(std::any sender, RoutedEventArgs *e)
	{
		ProteinDbObservableCollection->Clear();
	}

	void MainWindow::ResetTasksButton_Click(std::any sender, RoutedEventArgs *e)
	{
		tasksGroupBox->IsEnabled = true;
		ClearTasksButton->IsEnabled = true;
		DeleteSelectedTaskButton->IsEnabled = true;
		RunTasksButton->IsEnabled = true;
		addCalibrateTaskButton->IsEnabled = true;
		addGPTMDTaskButton->IsEnabled = true;
		addSearchTaskButton->IsEnabled = true;
		btnAddCrosslinkSearch->IsEnabled = true;
		ResetTasksButton->IsEnabled = false;
		OutputFolderTextBox->IsEnabled = true;

		dataGridSpectraFiles->IsReadOnly = false;
		dataGridProteinDatabases->IsReadOnly = false;

		AddXML->IsEnabled = true;
		ClearXML->IsEnabled = true;
		AddRaw->IsEnabled = true;
		ClearRaw->IsEnabled = true;
		BtnQuantSet->IsEnabled = true;

		LoadTaskButton->IsEnabled = true;

		tasksTreeView->DataContext = StaticTasksObservableCollection;
		UpdateSpectraFileGuiStuff();

		auto pathOfFirstSpectraFile = FileSystem::getDirectoryName(SpectraFilesObservableCollection->First().FilePath);
		OutputFolderTextBox->Text = FileSystem::combine(pathOfFirstSpectraFile, LR"($DATETIME)");
	}

	void MainWindow::TasksTreeView_MouseDoubleClick(std::any sender, MouseButtonEventArgs *e)
	{
		auto a = dynamic_cast<TreeView*>(sender);
//C# TO C++ CONVERTER TODO TASK: C++ has no equivalent to C# pattern variables in 'is' expressions:
//ORIGINAL LINE: if (a.SelectedItem is PreRunTask preRunTask)
		if (dynamic_cast<PreRunTask*>(a->SelectedItem) != nullptr preRunTask)
		{
			switch (preRunTask::metaMorpheusTask::TaskType)
			{
				case MyTask::Search:
				{

					auto searchDialog = new SearchTaskWindow(dynamic_cast<SearchTask*>(preRunTask::metaMorpheusTask));
					searchDialog->ShowDialog();
					preRunTask->DisplayName = L"Task" + std::to_wstring(StaticTasksObservableCollection->find(preRunTask) + 1) + L"-" + searchDialog->getTheTask()->getCommonParameters()->getTaskDescriptor();
					tasksTreeView::Items->Refresh();

					delete searchDialog;
					return;

				}
				case MyTask::Gptmd:
				{
					auto gptmddialog = new GptmdTaskWindow(dynamic_cast<GptmdTask*>(preRunTask::metaMorpheusTask));
					gptmddialog->ShowDialog();
					preRunTask->DisplayName = L"Task" + std::to_wstring(StaticTasksObservableCollection->find(preRunTask) + 1) + L"-" + gptmddialog->getTheTask()->getCommonParameters()->getTaskDescriptor();
					tasksTreeView::Items->Refresh();

					delete gptmddialog;
					return;

				}
				case MyTask::Calibrate:
				{
					auto calibratedialog = new CalibrateTaskWindow(dynamic_cast<CalibrationTask*>(preRunTask::metaMorpheusTask));
					calibratedialog->ShowDialog();
					preRunTask->DisplayName = L"Task" + std::to_wstring(StaticTasksObservableCollection->find(preRunTask) + 1) + L"-" + calibratedialog->getTheTask()->getCommonParameters()->getTaskDescriptor();
					tasksTreeView::Items->Refresh();

					delete calibratedialog;
					return;

				}
				case MyTask::XLSearch:
				{
					auto XLSearchdialog = new XLSearchTaskWindow(dynamic_cast<XLSearchTask*>(preRunTask::metaMorpheusTask));
					XLSearchdialog->ShowDialog();
					preRunTask->DisplayName = L"Task" + std::to_wstring(StaticTasksObservableCollection->find(preRunTask) + 1) + L"-" + XLSearchdialog->getTheTask()->getCommonParameters()->getTaskDescriptor();
					tasksTreeView::Items->Refresh();

					delete XLSearchdialog;
					return;
				}
			}
		}

//C# TO C++ CONVERTER TODO TASK: C++ has no equivalent to C# pattern variables in 'is' expressions:
//ORIGINAL LINE: if (a.SelectedItem is OutputFileForTreeView fileThing)
		if (dynamic_cast<OutputFileForTreeView*>(a->SelectedItem) != nullptr fileThing)
		{
			if (FileSystem::fileExists(fileThing::FullPath))
			{
				System::Diagnostics::Process::Start(fileThing::FullPath);
			}
			else
			{
				MessageBox::Show(L"File " + FileSystem::getFileName(fileThing::FullPath) + L" does not exist");
			}
		}
	}

	void MainWindow::LoadTaskButton_Click(std::any sender, RoutedEventArgs *e)
	{
		Microsoft::Win32::OpenFileDialog *openFileDialog1 = new Microsoft::Win32::OpenFileDialog();
		openFileDialog1->Filter = L"TOML files(*.toml)|*.toml";
		openFileDialog1->FilterIndex = 1;
		openFileDialog1->RestoreDirectory = true;
		openFileDialog1->Multiselect = true;
		if (openFileDialog1->ShowDialog() == true)
		{
			for (auto tomlFromSelected : openFileDialog1->FileNames.OrderBy([&] (std::any p)
			{
			delete openFileDialog1;
				return p;
			}))
			{
				AddAFile(tomlFromSelected);
			}
		}
		UpdateTaskGuiStuff();

		delete openFileDialog1;
	}

	void MainWindow::UpdateFileSpecificParamsDisplay(std::vector<std::wstring> &tomlLocations)
	{
		std::vector<std::wstring> fullPathofTomls = tomlLocations;

		for (auto file : SelectedRawFiles)
		{
			for (int j = 0; j < fullPathofTomls.Count(); j++)
			{
				if (Path::GetFileNameWithoutExtension(file->getFileName()) == Path::GetFileNameWithoutExtension(fullPathofTomls[j]))
				{
					if (FileSystem::fileExists(fullPathofTomls[j]))
					{
						TomlTable *fileSpecificSettings = Toml::ReadFile(fullPathofTomls[j], MetaMorpheusTask::tomlConfig);
						try
						{
							// parse to make sure toml is readable
							auto temp = new FileSpecificParameters(fileSpecificSettings);

							// toml is ok; display the file-specific settings in the gui
							file->SetParametersText(File::ReadAllText(fullPathofTomls[j]));

							delete temp;
						}
						catch (const MetaMorpheusException &e)
						{
							StringEventArgs tempVar(L"Problem parsing the file-specific toml " + FileSystem::getFileName(fullPathofTomls[j]) + L"; " + e->what() + L"; is the toml from an older version of MetaMorpheus?", nullptr);
							GuiWarnHandler(nullptr, &tempVar);
						}
					}
					else
					{
						file->SetParametersText(L"");
					}
				}
			}
		}
		UpdateSpectraFileGuiStuff();
		dataGridSpectraFiles::Items->Refresh();
	}

	void MainWindow::UpdateFileSpecificParamsDisplayJustAdded(const std::wstring &tomlLocation)
	{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		std::wstring assumedPathOfSpectraFileWithoutExtension = FileSystem::combine(Directory::GetParent(tomlLocation)->ToString(), Path::GetFileNameWithoutExtension(tomlLocation));

		for (int i = 0; i < SpectraFilesObservableCollection->Count; i++)
		{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			std::wstring thisFilesPathWihoutExtension = FileSystem::combine(Directory::GetParent(SpectraFilesObservableCollection[i]->getFilePath())->ToString(), Path::GetFileNameWithoutExtension(SpectraFilesObservableCollection[i]->getFilePath()));
			if (FileSystem::fileExists(tomlLocation) && assumedPathOfSpectraFileWithoutExtension == thisFilesPathWihoutExtension)
			{
				TomlTable *fileSpecificSettings = Toml::ReadFile(tomlLocation, MetaMorpheusTask::tomlConfig);
				try
				{
					// parse to make sure toml is readable
					auto temp = new FileSpecificParameters(fileSpecificSettings);

					// toml is ok; display the file-specific settings in the gui
					SpectraFilesObservableCollection[i]->SetParametersText(File::ReadAllText(tomlLocation));

					delete temp;
				}
				catch (const MetaMorpheusException &e)
				{
					StringEventArgs tempVar(L"Problem parsing the file-specific toml " + FileSystem::getFileName(tomlLocation) + L"; " + e->what() + L"; is the toml from an older version of MetaMorpheus?", nullptr);
					GuiWarnHandler(nullptr, &tempVar);
				}
			}
		}
		UpdateSpectraFileGuiStuff();
		dataGridSpectraFiles::Items->Refresh();
	}

	void MainWindow::AddSelectedSpectra(std::any sender, RoutedEventArgs *e)
	{
		DataGridRow *obj = std::any_cast<DataGridRow*>(sender);

		RawDataForDataGrid *ok = static_cast<RawDataForDataGrid*>(obj->DataContext);
		SelectedRawFiles->Add(ok);
		UpdateSpectraFileGuiStuff();
	}

	void MainWindow::RemoveSelectedSpectra(std::any sender, RoutedEventArgs *e)
	{
		DataGridRow *obj = std::any_cast<DataGridRow*>(sender);
		RawDataForDataGrid *ok = static_cast<RawDataForDataGrid*>(obj->DataContext);
		SelectedRawFiles->Remove(ok);
		UpdateSpectraFileGuiStuff();
	}

	void MainWindow::MenuItem_Click(std::any sender, RoutedEventArgs *e)
	{
		System::Diagnostics::Process::Start(LR"(https://github.com/smith-chem-wisc/MetaMorpheus/wiki)");
	}

	void MainWindow::MenuItem_YouTube(std::any sender, RoutedEventArgs *e)
	{
		System::Diagnostics::Process::Start(LR"(https://www.youtube.com/playlist?list=PLVk5tTSZ1aWlhNPh7jxPQ8pc0ElyzSUQb)");
	}

	void MainWindow::MenuItem_ProteomicsNewsBlog(std::any sender, RoutedEventArgs *e)
	{
		System::Diagnostics::Process::Start(LR"(https://proteomicsnews.blogspot.com/)");
	}

	void MainWindow::MenuItem_Click_1(std::any sender, RoutedEventArgs *e)
	{
		auto globalSettingsDialog = new GlobalSettingsWindow();
		globalSettingsDialog->ShowDialog();

		delete globalSettingsDialog;
	}

	bool MainWindow::DatabaseExists(ObservableCollection<ProteinDbForDataGrid*> *pDOC, ProteinDbForDataGrid *uuu)
	{
		for (auto pdoc : pDOC)
		{
			if (pdoc->getFilePath() == uuu->getFilePath())
			{
			return true;
			}
		}
		return false;
	}

	bool MainWindow::SpectraFileExists(ObservableCollection<RawDataForDataGrid*> *rDOC, RawDataForDataGrid *zzz)
	{
		for (auto rdoc : rDOC)
		{
			if (rdoc->getFileName() == zzz->getFileName())
			{
			return true;
			}
		}
		return false;
	}

	void MainWindow::CancelButton_Click(std::any sender, RoutedEventArgs *e)
	{
		std::wstring grammar = StaticTasksObservableCollection->Count <= 1 ? L"this task" : L"these tasks";
		if (MessageBox::Show(L"Are you sure you want to cancel " + grammar + L"?", L"Cancel Tasks", MessageBoxButton::OKCancel) == MessageBoxResult::OK)
		{
			GlobalVariables::setStopLoops(true);
			CancelButton->IsEnabled = false;
			notificationsTextBox::AppendText(L"Canceling...\n");
		}
	}

	void MainWindow::ChangeFileParameters_Click(std::any sender, RoutedEventArgs *e)
	{
		auto dialog = new FileSpecificParametersWindow(SelectedRawFiles);
		if (dialog->ShowDialog() == true)
		{
			auto tomlPathsForSelectedFiles = SelectedRawFiles->Select([&] (std::any p)
			{
			delete dialog;
				return FileSystem::combine(Directory::GetParent(p::FilePath)->ToString(), Path::GetFileNameWithoutExtension(p::FileName)) + L".toml";
			}).ToList();
			UpdateFileSpecificParamsDisplay(tomlPathsForSelectedFiles.ToArray());
		}

		delete dialog;
	}

	void MainWindow::BtnQuantSet_Click(std::any sender, RoutedEventArgs *e)
	{
		auto dialog = new ExperimentalDesignWindow(SpectraFilesObservableCollection);
		dialog->ShowDialog();

		delete dialog;
	}

	void MainWindow::MenuItem_Click_2(std::any sender, RoutedEventArgs *e)
	{
		try
		{
			GetVersionNumbersFromWeb();
		}
		catch (const std::runtime_error &ex)
		{
			StringEventArgs tempVar(L"Could not get newest MM version from web: " + ex.what(), nullptr);
			GuiWarnHandler(nullptr, &tempVar);
			return;
		}

		if (GlobalVariables::getMetaMorpheusVersion() == getNewestKnownVersion())
		{
			MessageBox::Show(L"You have the most updated version!");
		}
		else
		{
			try
			{
				MetaUpdater *newwind = new MetaUpdater();
				newwind->ShowDialog();

				delete newwind;
			}
			catch (const std::runtime_error &ex)
			{
				MessageBox::Show(ex.what());
			}
		}
	}

	void MainWindow::MenuItem_Click_4(std::any sender, RoutedEventArgs *e)
	{
		std::wstring mailto = StringHelper::formatSimple(L"mailto:{0}?Subject=MetaMorpheus. Issue:", L"mm_support@chem.wisc.edu");
		System::Diagnostics::Process::Start(mailto);
	}

	void MainWindow::MenuItem_Click_5(std::any sender, RoutedEventArgs *e)
	{
		System::Diagnostics::Process::Start(LR"(https://github.com/smith-chem-wisc/MetaMorpheus/issues/new)");
	}

	void MainWindow::MenuItem_Twitter(std::any sender, RoutedEventArgs *e)
	{
		System::Diagnostics::Process::Start(LR"(https://twitter.com/Smith_Chem_Wisc)");
	}

	void MainWindow::MenuItem_Slack(std::any sender, RoutedEventArgs *e)
	{
		System::Diagnostics::Process::Start(LR"(https://join.slack.com/t/smith-chem-public/shared_invite/enQtNDYzNTM5Mzg5NzY0LTRiYWQ5MzVmYmExZWIyMTcyZmNlODJjMWI0YjVhNGM2MmQ2NjE4ZDAzNmM4NWYxMDFhNTQyNDBiM2E0MWE0NGU)");
	}

	void MainWindow::MenuItem_Click_6(std::any sender, RoutedEventArgs *e)
	{
		System::Diagnostics::Process::Start(GlobalVariables::getDataDir());
	}

	void MainWindow::MenuItem_Click_3(std::any sender, RoutedEventArgs *e)
	{
		System::Diagnostics::Process::Start(FileSystem::combine(GlobalVariables::getDataDir(), LR"(settings.toml)"));
		Application::Current->Shutdown();
	}

	void MainWindow::MenuItem_Click_7(std::any sender, RoutedEventArgs *e)
	{
		System::Diagnostics::Process::Start(FileSystem::combine(GlobalVariables::getDataDir(), LR"(GUIsettings.toml)"));
		Application::Current->Shutdown();
	}

	void MainWindow::MetaDrawMenu_Click(std::any sender, RoutedEventArgs *e)
	{
		MetaDraw *metaDrawGui = new MetaDraw();
		metaDrawGui->Show();

		delete metaDrawGui;
	}

	void MainWindow::PrintErrorsReadingMods()
	{
		// print any error messages reading the mods to the notifications area
		for (auto error : GlobalVariables::ErrorsReadingMods)
		{
			StringEventArgs tempVar(error, nullptr);
			GuiWarnHandler(nullptr, &tempVar);
		}
		GlobalVariables::ErrorsReadingMods.clear();
	}

	void MainWindow::AddContaminantXML_Click(std::any sender, RoutedEventArgs *e)
	{
		std::vector<std::wstring> contaminantFiles = Directory::GetFiles(FileSystem::combine(GlobalVariables::getDataDir(), L"Contaminants"));
		for (auto contaminantFile : contaminantFiles)
		{
			AddAFile(contaminantFile);
		}
		dataGridProteinDatabases::Items->Refresh();
	}

	void MainWindow::AddCustomMod_Click(std::any sender, RoutedEventArgs *e)
	{
		auto dialog = new CustomModButtonWindow();
		dialog->ShowDialog();

		delete dialog;
	}

	void MainWindow::MainWindow_Closing(std::any sender, CancelEventArgs *e)
	{
		if (!GuiGlobalParams->getDisableCloseWindow() && !GlobalVariables::getMetaMorpheusVersion().find(L"DEBUG") != std::wstring::npos)
		{
			e->Cancel = true;
			auto exit = CustomMsgBox::Show(L"Exit MetaMorpheus", L"Are you sure you want to exit MetaMorpheus?", L"Yes", L"No", L"Yes and don't ask me again");

			if (exit == MessageBoxResult::Yes)
			{
				e->Cancel = false;
			}
			else if (exit == MessageBoxResult::OK)
			{
				GuiGlobalParams->setDisableCloseWindow(true);
				Toml::WriteFile(GuiGlobalParams, FileSystem::combine(GlobalVariables::getDataDir(), LR"(GUIsettings.toml)"), MetaMorpheusTask::tomlConfig);
				e->Cancel = false;
			}
		}
	}
}
