#include "MetaMorpheusGUI.MetaUpdater.h"
#include "MetaMorpheusGUI.MainWindow.h"
#include "../EngineLayer/GlobalVariables.h"

using namespace EngineLayer;
using namespace Newtonsoft::Json::Linq;

namespace MetaMorpheusGUI
{

	MetaUpdater::MetaUpdater()
	{
		InitializeComponent();
		lbl->Text = L"A newer version: " + MainWindow::getNewestKnownVersion() + L" is available!";
		ReleaseHandler();
	}

	void MetaUpdater::InstallerClicked(std::any sender, RoutedEventArgs *e)
	{
		DialogResult = true;
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var client = new WebClient())
		{
			auto client = WebClient();
			auto uri = new Uri(LR"(https://github.com/smith-chem-wisc/MetaMorpheus/releases/download/)" + MainWindow::getNewestKnownVersion() + LR"(/MetaMorpheusInstaller.msi)");

			try
			{
				auto tempDownloadLocation = FileSystem::combine(System::IO::Path::GetTempPath(), L"MetaMorpheusInstaller.msi");
				client.DownloadFile(uri, tempDownloadLocation);
				Process *p = new Process();
				p->StartInfo->FileName = tempDownloadLocation;
				Application::Current->Shutdown();
				p->Start();

				delete p;
			}
			catch (const std::runtime_error &ex)
			{
				MessageBox::Show(ex.what());
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete uri' statement was not added since uri was passed to a method or constructor. Handle memory management manually.
		}
	}

	void MetaUpdater::ReleaseHandler()
	{
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var client = new HttpClient())
		{
			auto client = HttpClient();
			client.DefaultRequestHeaders->Add(L"User-Agent", L"Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; WOW64; Trident/6.0)");

//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var response = client.GetAsync("https://api.github.com/repos/smith-chem-wisc/MetaMorpheus/releases").Result)
			{
				auto response = client.GetAsync(L"https://api.github.com/repos/smith-chem-wisc/MetaMorpheus/releases").Result;
				std::wstring json = response.Content::ReadAsStringAsync().Result;
				JArray *GitArray = JArray::Parse(json);
				auto currV = GetVersionNumber(GlobalVariables::getMetaMorpheusVersion());
				StringBuilder *allVersionsText = new StringBuilder();
				for (JObject *obj : GitArray->Children<JObject*>())
				{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					auto checkVersion = GetVersionNumber(obj->SelectToken(L"tag_name").ToString());
					if (checkVersion->Item1 < currV->Item1 || (checkVersion->Item1 == currV->Item1 && checkVersion->Item2 < currV->Item2) || (checkVersion->Item1 == currV->Item1 && checkVersion->Item2 == currV->Item2 && checkVersion->Item3 <= currV->Item3))
					{
						break;
					}
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					allVersionsText->appendLine(obj->SelectToken(L"tag_name").ToString());
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					allVersionsText->appendLine(obj->SelectToken(L"body").ToString());
					allVersionsText->appendLine();
				}
				releases->Text = allVersionsText->toString();

				delete allVersionsText;
			}
		}
	}

	void MetaUpdater::PortableClicked(std::any semder, RoutedEventArgs *e)
	{
		DialogResult = true;
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var client = new WebClient())
		{
			auto client = WebClient();
			auto uri = new Uri(LR"(https://github.com/smith-chem-wisc/MetaMorpheus/releases/download/)" + MainWindow::getNewestKnownVersion() + LR"(/MetaMorpheusGuiDotNetFrameworkAppveyor.zip)");

			try
			{
				auto tempDownloadLocation = FileSystem::combine(System::IO::Path::GetTempPath(), L"MetaMorpheusGuiDotNetFrameworkAppveyor.zip");
				client.DownloadFile(uri, tempDownloadLocation);
				Process *p = new Process();
				p->StartInfo->FileName = tempDownloadLocation;
				Application::Current->Shutdown();
				p->Start();

				delete p;
			}
			catch (const std::runtime_error &ex)
			{
				MessageBox::Show(ex.what());
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete uri' statement was not added since uri was passed to a method or constructor. Handle memory management manually.
		}
	}

	void MetaUpdater::NoClicked(std::any semder, RoutedEventArgs *e)
	{
		DialogResult = false;
	}
}
