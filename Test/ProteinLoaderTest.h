#pragma once

#include "../TaskLayer/MetaMorpheusTask.h"
#include <string>
#include <vector>
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class DbForTask; }
namespace TaskLayer { class FileSpecificParameters; }
namespace TaskLayer { class MyTaskResults; }

using namespace EngineLayer;
using namespace NUnit::Framework;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public class ProteinLoaderTest
	class ProteinLoaderTest
	{
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public void TestProteinLoad()
		void TestProteinLoad();

	public:
		class ProteinLoaderTask : public MetaMorpheusTask
		{
		public:
			ProteinLoaderTask(const std::wstring &x);

		protected:
			ProteinLoaderTask();

		public:
			void Run(const std::wstring &dbPath);

		protected:
			MyTaskResults *RunSpecific(const std::wstring &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::wstring> &currentRawFileList, const std::wstring &taskId, std::vector<FileSpecificParameters*> &fileSettingsList) override;
		};
	};
}
