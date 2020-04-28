#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <any>
#include "stringhelper.h"
#include "tangible_filesystem.h"

#include "../EngineLayer/EventArgs/SingleFileEventArgs.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineFinishedEventArgs.h"
#include "../TaskLayer/EventArgs/SingleTaskEventArgs.h"

using namespace EngineLayer;
using namespace TaskLayer;

//using namespace Proteomics;

namespace MetaMorpheusCommandLine
{
    class IndentedTextWriter {
    public:
        int Indent=0;
        void Write(std::string);
        void WriteLine(std::string);
        void WriteLine(void);
    };


    
    class Program final
	{
        public:
            static void Main( int argc, char *argv[]);

	private:
		static bool InProgress;

		//static System::CodeDom::Compiler::IndentedTextWriter *MyWriter;
                static IndentedTextWriter *MyWriter;

		static void WriteMultiLineIndented(const std::string &toWrite);

		static bool IsContaminant(const std::string &b);

		static void MyTaskEngine_startingSingleTaskHandler(SingleTaskEventArgs e);

		static void MyTaskEngine_finishedWritingFileHandler(SingleFileEventArgs e);

		static void MyTaskEngine_finishedSingleTaskHandler(SingleTaskEventArgs e);

		static void MyEngine_startingSingleEngineHandler(SingleEngineEventArgs e);

		static void MyEngine_finishedSingleEngineHandler(SingleEngineFinishedEventArgs e);

		static void MyEngine_outProgressHandler(ProgressEventArgs e);

		static void WarnHandler(StringEventArgs e);

		static void LogHandler(StringEventArgs e);

                static void print_config ( std::vector<std::string> taskname,
                                           std::vector<std::string> metataskname,
                                           std::vector<std::string> outputfolder,
                                           std::vector<std::string> spectra,
                                           std::vector<std::string> dbases );

                
	public:
		class ApplicationArguments
		{
		private:
			std::vector<std::string> privateTasks;
			std::vector<std::string> privateDatabases;
			std::vector<std::string> privateSpectra;
			std::vector<std::string> privateMetaTasks;
			std::string privateOutputFolder;

		public:
			std::vector<std::string> getTasks() const;
			void setTasks(const std::vector<std::string> &value);
			std::vector<std::string> getDatabases() const;
			void setDatabases(const std::vector<std::string> &value);
			std::vector<std::string> getSpectra() const;
			void setSpectra(const std::vector<std::string> &value);
			std::vector<std::string> getMetaTasks() const;
			void setMetaTasks(const std::vector<std::string> &value);
			std::string getOutputFolder() const;
			void setOutputFolder(const std::string &value);
		};
	};
}
