#pragma once

#include <string>
#include <unordered_map>
#include <optional>
#include <mutex>
#include "tangible_event.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class CommonParameters; }

using namespace EngineLayer;
using namespace IO::MzML;
using namespace IO::Mgf;

#if defined(NETFRAMEWORK)

using namespace IO::Thermo;

#else
#endif

using namespace MassSpectrometry;

namespace TaskLayer
{
	class MyFileManager
	{
	public:
		enum class ThermoMsFileReaderVersionCheck
		{
			DllsNotFound,
			IncorrectVersion,
			CorrectVersion,
			SomeDllsMissing
		};

	private:
		const bool DisposeOfFileWhenDone;
		const std::unordered_map<std::string, MsDataFile*> MyMsDataFiles = std::unordered_map<std::string, MsDataFile*>();
		std::mutex FileLoadingLock;
		static const std::string AssumedThermoMsFileReaderDllPath;
		static const std::string DesiredFileIoVersion;
		static const std::string DesiredFregistryVersion;
		static const std::string DesiredXRawFileVersion;

	public:
		MyFileManager(bool disposeOfFileWhenDone);

		static TangibleEvent<EventHandler<StringEventArgs>> *WarnHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		bool SeeIfOpen(const std::string &path);

		static ThermoMsFileReaderVersionCheck ValidateThermoMsFileReaderVersion();

		MsDataFile *LoadFile(const std::string &origDataFile, std::optional<int> &topNpeaks, std::optional<double> &minRatio, bool trimMs1Peaks, bool trimMsMsPeaks, CommonParameters *commonParameters);

		void DoneWithFile(const std::string &origDataFile);

	private:
		void Warn(const std::string &v);
	};
}
