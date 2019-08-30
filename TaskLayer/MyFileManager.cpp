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

const std::string MyFileManager::AssumedThermoMsFileReaderDllPath = LR"(C:\Program Files\Thermo\MSFileReader)";
const std::string MyFileManager::DesiredFileIoVersion = "3.0";
const std::string MyFileManager::DesiredFregistryVersion = "3.0";
const std::string MyFileManager::DesiredXRawFileVersion = "3.0.29.0";

	MyFileManager::MyFileManager(bool disposeOfFileWhenDone) : DisposeOfFileWhenDone(disposeOfFileWhenDone)
	{
	}

	bool MyFileManager::SeeIfOpen(const std::string &path)
	{
		return (MyMsDataFiles.find(path) != MyMsDataFiles.end() && MyMsDataFiles[path] != nullptr);
	}

	MyFileManager::ThermoMsFileReaderVersionCheck MyFileManager::ValidateThermoMsFileReaderVersion()
	{
		std::string fileIoAssumedPath = FileSystem::combine(AssumedThermoMsFileReaderDllPath, "Fileio_x64.dll");
		std::string fregistryAssumedPath = FileSystem::combine(AssumedThermoMsFileReaderDllPath, "fregistry_x64.dll");
		std::string xRawFileAssumedPath = FileSystem::combine(AssumedThermoMsFileReaderDllPath, "XRawfile2_x64.dll");

		if (FileSystem::fileExists(fileIoAssumedPath) && FileSystem::fileExists(fregistryAssumedPath) && FileSystem::fileExists(xRawFileAssumedPath))
		{
			std::string fileIoVersion = FileVersionInfo::GetVersionInfo(fileIoAssumedPath)->FileVersion;
			std::string fregistryVersion = FileVersionInfo::GetVersionInfo(fregistryAssumedPath)->FileVersion;
			std::string xRawFileVersion = FileVersionInfo::GetVersionInfo(xRawFileAssumedPath)->FileVersion;

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

	MsDataFile *MyFileManager::LoadFile(const std::string &origDataFile, std::optional<int> &topNpeaks, std::optional<double> &minRatio, bool trimMs1Peaks, bool trimMsMsPeaks, CommonParameters *commonParameters)
	{
		FilteringParams *filter = new FilteringParams(topNpeaks, minRatio, 1, trimMs1Peaks, trimMsMsPeaks);
		MsDataFile value;
		std::unordered_map<std::string, MsDataFile*>::const_iterator MyMsDataFiles_iterator = MyMsDataFiles.find(origDataFile);
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
			if (Path::GetExtension(origDataFile).Equals(".mzM", StringComparison::OrdinalIgnoreCase))
			{
				MyMsDataFiles[origDataFile] = Mzml::LoadAllStaticData(origDataFile, filter, commonParameters->getMaxThreadsToUsePerFile());
			}
//C# TO C++ CONVERTER TODO TASK: The following .NET 'String.Equals' reference is not converted:
			else if (Path::GetExtension(origDataFile).Equals(".mgf", StringComparison::OrdinalIgnoreCase))
			{
				MyMsDataFiles[origDataFile] = Mgf::LoadAllStaticData(origDataFile, filter);
			}
			else
			{
	#if defined(NETFRAMEWORK)
				MyMsDataFiles[origDataFile] = ThermoStaticData::LoadAllStaticData(origDataFile, filter);
	#else
				Warn("No capability for reading " + origDataFile);
	#endif
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete filter' statement was not added since filter was passed to a method or constructor. Handle memory management manually.
			return MyMsDataFiles[origDataFile];
		}

//C# TO C++ CONVERTER TODO TASK: A 'delete filter' statement was not added since filter was passed to a method or constructor. Handle memory management manually.
	}

	void MyFileManager::DoneWithFile(const std::string &origDataFile)
	{
		if (DisposeOfFileWhenDone)
		{
			MyMsDataFiles[origDataFile] = nullptr;
		}
	}

	void MyFileManager::Warn(const std::string &v)
	{
		StringEventArgs tempVar(v, nullptr);
		WarnHandler +== nullptr ? nullptr : WarnHandler::Invoke(this, &tempVar);
	}
}
