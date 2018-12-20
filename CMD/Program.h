#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <any>
#include "stringhelper.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class SingleTaskEventArgs; }
namespace EngineLayer { class SingleFileEventArgs; }
namespace EngineLayer { class SingleEngineEventArgs; }
namespace EngineLayer { class SingleEngineFinishedEventArgs; }
namespace EngineLayer { class ProgressEventArgs; }
namespace EngineLayer { class StringEventArgs; }

using namespace EngineLayer;
//using namespace Fclp;
//using namespace Nett;
//using namespace Proteomics;
using namespace TaskLayer;

namespace MetaMorpheusCommandLine
{
	class Program final
	{
	private:
		static bool InProgress;

		static System::CodeDom::Compiler::IndentedTextWriter *MyWriter;

		static void Main(std::vector<std::wstring> &args);

		static void WriteMultiLineIndented(const std::wstring &toWrite);

		static bool IsContaminant(const std::wstring &b);

		static void MyTaskEngine_startingSingleTaskHander(std::any sender, SingleTaskEventArgs *e);

		static void MyTaskEngine_finishedWritingFileHandler(std::any sender, SingleFileEventArgs *e);

		static void MyTaskEngine_finishedSingleTaskHandler(std::any sender, SingleTaskEventArgs *e);

		static void MyEngine_startingSingleEngineHander(std::any sender, SingleEngineEventArgs *e);

		static void MyEngine_finishedSingleEngineHandler(std::any sender, SingleEngineFinishedEventArgs *e);

		static void MyEngine_outProgressHandler(std::any sender, ProgressEventArgs *e);

		static void WarnHandler(std::any sender, StringEventArgs *e);

		static void LogHandler(std::any sender, StringEventArgs *e);

	public:
		class ApplicationArguments
		{
		private:
			std::vector<std::wstring> privateTasks;
			std::vector<std::wstring> privateDatabases;
			std::vector<std::wstring> privateSpectra;
			std::vector<std::wstring> privateMetaTasks;
			std::wstring privateOutputFolder;

		public:
			std::vector<std::wstring> getTasks() const;
			void setTasks(const std::vector<std::wstring> &value);
			std::vector<std::wstring> getDatabases() const;
			void setDatabases(const std::vector<std::wstring> &value);
			std::vector<std::wstring> getSpectra() const;
			void setSpectra(const std::vector<std::wstring> &value);
			std::vector<std::wstring> getMetaTasks() const;
			void setMetaTasks(const std::vector<std::wstring> &value);
			std::wstring getOutputFolder() const;
			void setOutputFolder(const std::wstring &value);
		};
	};
}
