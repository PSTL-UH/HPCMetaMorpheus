#include "CalibrationTests.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
#include "../TaskLayer/DbForTask.h"

using namespace NUnit::Framework;
using namespace TaskLayer;

namespace Test
{

	void CalibrationTests::ExperimentalDesignCalibrationTest()
	{
		CalibrationTask *calibrationTask = new CalibrationTask();
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestCalibration)");
		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\SmallCalibratible_Yeast.mzML)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\smalldb.fasta)");
		std::wstring filePath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\ExperimentalDesign.tsv)");
		FileSystem::createDirectory(outputFolder);
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(filePath))
		{
			StreamWriter output = StreamWriter(filePath);
			output.WriteLine(L"FileName\tCondition\tBiorep\tFraction\tTechrep");
			output.WriteLine(std::wstring(L"SmallCalibratible_Yeast") + L"\t" + L"condition" + L"\t" + L"1" + L"\t" + L"1" + L"\t" + L"1");
		}

		calibrationTask->RunTask(outputFolder, {new DbForTask(myDatabase, false)}, {myFile}, L"test");
		auto expDesignPath = FileSystem::combine(outputFolder, LR"(ExperimentalDesign.tsv)");
		auto expDesign = File::ReadAllLines(expDesignPath);
		Assert::That(expDesign[1].find(L"SmallCalibratible_Yeast-calib") != std::wstring::npos);
		Assert::That(FileSystem::fileExists(FileSystem::combine(outputFolder, LR"(SmallCalibratible_Yeast-calib.mzML)")));
		Assert::That(FileSystem::fileExists(FileSystem::combine(outputFolder, LR"(SmallCalibratible_Yeast-calib.toml)")));
		auto lines = File::ReadAllLines(FileSystem::combine(outputFolder, LR"(SmallCalibratible_Yeast-calib.toml)"));
		auto tolerance = Regex::Match(lines[0], LR"(\d+\.\d*)")->Value;
		auto tolerance1 = Regex::Match(lines[1], LR"(\d+\.\d*)")->Value;
		double tol;
		Assert::That(double::TryParse(tolerance, tol) == true);
		double tol1;
		Assert::That(double::TryParse(tolerance1, tol1) == true);
		Assert::That(lines[0].find(L"PrecursorMassTolerance") != std::wstring::npos);
		Assert::That(lines[1].find(L"ProductMassTolerance") != std::wstring::npos);
		File::Delete(filePath);
		Directory::Delete(outputFolder, true);
		Directory::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(Task Settings)"), true);

		delete calibrationTask;
	}
}
