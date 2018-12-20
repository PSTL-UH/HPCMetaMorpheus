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
		const std::unordered_map<std::wstring, MsDataFile*> MyMsDataFiles = std::unordered_map<std::wstring, MsDataFile*>();
		std::mutex FileLoadingLock;
		static const std::wstring AssumedThermoMsFileReaderDllPath;
		static const std::wstring DesiredFileIoVersion;
		static const std::wstring DesiredFregistryVersion;
		static const std::wstring DesiredXRawFileVersion;

	public:
		MyFileManager(bool disposeOfFileWhenDone);

		static TangibleEvent<EventHandler<StringEventArgs>> *WarnHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		bool SeeIfOpen(const std::wstring &path);

		static ThermoMsFileReaderVersionCheck ValidateThermoMsFileReaderVersion();

		MsDataFile *LoadFile(const std::wstring &origDataFile, std::optional<int> &topNpeaks, std::optional<double> &minRatio, bool trimMs1Peaks, bool trimMsMsPeaks, CommonParameters *commonParameters);

		void DoneWithFile(const std::wstring &origDataFile);

	private:
		void Warn(const std::wstring &v);
	};
}
