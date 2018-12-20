#include "MyFileManager.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"

using namespace EngineLayer;
using namespace IO::MzML;
using namespace IO::Mgf;
#if defined(NETFRAMEWORK)
using namespace IO::Thermo;
#endif
using namespace MassSpectrometry;

namespace TaskLayer
{

const std::wstring MyFileManager::AssumedThermoMsFileReaderDllPath = LR"(C:\Program Files\Thermo\MSFileReader)";
const std::wstring MyFileManager::DesiredFileIoVersion = L"3.0";
const std::wstring MyFileManager::DesiredFregistryVersion = L"3.0";
const std::wstring MyFileManager::DesiredXRawFileVersion = L"3.0.29.0";

	MyFileManager::MyFileManager(bool disposeOfFileWhenDone) : DisposeOfFileWhenDone(disposeOfFileWhenDone)
	{
	}

	bool MyFileManager::SeeIfOpen(const std::wstring &path)
	{
		return (MyMsDataFiles.find(path) != MyMsDataFiles.end() && MyMsDataFiles[path] != nullptr);
	}

	MyFileManager::ThermoMsFileReaderVersionCheck MyFileManager::ValidateThermoMsFileReaderVersion()
	{
		std::wstring fileIoAssumedPath = FileSystem::combine(AssumedThermoMsFileReaderDllPath, L"Fileio_x64.dll");
		std::wstring fregistryAssumedPath = FileSystem::combine(AssumedThermoMsFileReaderDllPath, L"fregistry_x64.dll");
		std::wstring xRawFileAssumedPath = FileSystem::combine(AssumedThermoMsFileReaderDllPath, L"XRawfile2_x64.dll");

		if (FileSystem::fileExists(fileIoAssumedPath) && FileSystem::fileExists(fregistryAssumedPath) && FileSystem::fileExists(xRawFileAssumedPath))
		{
			std::wstring fileIoVersion = FileVersionInfo::GetVersionInfo(fileIoAssumedPath)->FileVersion;
			std::wstring fregistryVersion = FileVersionInfo::GetVersionInfo(fregistryAssumedPath)->FileVersion;
			std::wstring xRawFileVersion = FileVersionInfo::GetVersionInfo(xRawFileAssumedPath)->FileVersion;

			if (fileIoVersion == DesiredFileIoVersion && fregistryVersion == DesiredFregistryVersion && xRawFileVersion == DesiredXRawFileVersion)
			{
				return ThermoMsFileReaderVersionCheck::CorrectVersion;
			}
			else
			{
				return ThermoMsFileReaderVersionCheck::IncorrectVersion;
			}
		}
		else if (FileSystem::fileExists(fileIoAssumedPath) || FileSystem::fileExists(fregistryAssumedPath) || FileSystem::fileExists(xRawFileAssumedPath))
		{
			return ThermoMsFileReaderVersionCheck::SomeDllsMissing;
		}

		return ThermoMsFileReaderVersionCheck::DllsNotFound;
	}

	MsDataFile *MyFileManager::LoadFile(const std::wstring &origDataFile, std::optional<int> &topNpeaks, std::optional<double> &minRatio, bool trimMs1Peaks, bool trimMsMsPeaks, CommonParameters *commonParameters)
	{
		FilteringParams *filter = new FilteringParams(topNpeaks, minRatio, 1, trimMs1Peaks, trimMsMsPeaks);
		MsDataFile value;
		std::unordered_map<std::wstring, MsDataFile*>::const_iterator MyMsDataFiles_iterator = MyMsDataFiles.find(origDataFile);
		if (MyMsDataFiles_iterator != MyMsDataFiles.end() && value != nullptr)
		{
			value = MyMsDataFiles_iterator->second;

			delete filter;
			return value;
		}
		else
		{
			value = MyMsDataFiles_iterator->second;
		}

		{
		// By now know that need to load this file!!!
			std::lock_guard<std::mutex> lock(FileLoadingLock);
//C# TO C++ CONVERTER TODO TASK: The following .NET 'String.Equals' reference is not converted:
			if (Path::GetExtension(origDataFile).Equals(L".mzML", StringComparison::OrdinalIgnoreCase))
			{
				MyMsDataFiles[origDataFile] = Mzml::LoadAllStaticData(origDataFile, filter, commonParameters->getMaxThreadsToUsePerFile());
			}
//C# TO C++ CONVERTER TODO TASK: The following .NET 'String.Equals' reference is not converted:
			else if (Path::GetExtension(origDataFile).Equals(L".mgf", StringComparison::OrdinalIgnoreCase))
			{
				MyMsDataFiles[origDataFile] = Mgf::LoadAllStaticData(origDataFile, filter);
			}
			else
			{
	#if defined(NETFRAMEWORK)
				MyMsDataFiles[origDataFile] = ThermoStaticData::LoadAllStaticData(origDataFile, filter);
	#else
				Warn(L"No capability for reading " + origDataFile);
	#endif
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete filter' statement was not added since filter was passed to a method or constructor. Handle memory management manually.
			return MyMsDataFiles[origDataFile];
		}

//C# TO C++ CONVERTER TODO TASK: A 'delete filter' statement was not added since filter was passed to a method or constructor. Handle memory management manually.
	}

	void MyFileManager::DoneWithFile(const std::wstring &origDataFile)
	{
		if (DisposeOfFileWhenDone)
		{
			MyMsDataFiles[origDataFile] = nullptr;
		}
	}

	void MyFileManager::Warn(const std::wstring &v)
	{
		StringEventArgs tempVar(v, nullptr);
		WarnHandler +== nullptr ? nullptr : WarnHandler::Invoke(this, &tempVar);
	}
}
