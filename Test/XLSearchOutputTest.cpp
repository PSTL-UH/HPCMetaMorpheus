#include "XLSearchOutputTest.h"
#include "../TaskLayer/XLSearchTask/TaskLayer.XLSearchTask.h"
#include "../TaskLayer/DbForTask.h"

using namespace NUnit::Framework;
using namespace TaskLayer;

namespace Test
{

	void XLSearchOutputTest::WriteTsvTest()
	{
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XlOutputTest1)");
		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XlTestData\BSA_DSS_23747.mzML)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(XlTestData\BSA.fasta)");

		FileSystem::createDirectory(outputFolder);

		XLSearchTask *xLSearch = new XLSearchTask();
		xLSearch->RunTask(outputFolder, {new DbForTask(myDatabase, false)}, {myFile}, L"test");

		auto resultsPath = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(XL_Interlinks.tsv)"));
		auto sections = StringHelper::split(resultsPath[1], L'\t');
		Assert::That(resultsPath.size() == 2);
		Assert::That(sections.size() == 45);
		Directory::Delete(outputFolder, true);

		delete xLSearch;
	}
}
