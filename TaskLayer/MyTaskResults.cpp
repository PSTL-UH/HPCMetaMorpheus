#include "MyTaskResults.h"
#include "DbForTask.h"
#include "MetaMorpheusTask.h"


namespace TaskLayer
{

	MyTaskResults::MyTaskResults(MetaMorpheusTask *s) : resultTexts(std::vector<std::wstring>())
	{
	}

	std::wstring MyTaskResults::ToString()
	{
		StringBuilder *sb = new StringBuilder();
		sb->appendLine(L"Time to run task: " + Time);
		sb->appendLine();
		sb->appendLine();
		sb->appendLine(L"--------------------------------------------------");
		if ((NewSpectra.size() > 0 && NewSpectra.Any()) || (NewDatabases.size() > 0 && NewDatabases.Any()))
		{
			sb->appendLine();
			sb->appendLine();
			sb->appendLine(L"New files:");
			if (NewSpectra.size() > 0 && NewSpectra.Any())
			{
				sb->appendLine(L"New spectra: ");
				sb->appendLine();
				sb->appendLine(std::wstring::Join(L"\r\n" + L"\t", NewSpectra));
			}
			if (NewDatabases.size() > 0 && NewDatabases.Any())
			{
				sb->appendLine(L"New databases: ");
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				sb->appendLine(std::wstring::Join(L"\r\n" + L"\t", NewDatabases.Select([&] (std::any b)
				{
					b::FilePath;
				})).ToString());
			}
			sb->appendLine();
			sb->appendLine();
			sb->appendLine(L"--------------------------------------------------");
		}
		sb->appendLine();
		sb->appendLine();
		sb->appendLine(niceText->toString());
		sb->appendLine();
		sb->appendLine();
		sb->appendLine(L"--------------------------------------------------");
		sb->appendLine();
		sb->appendLine();
		sb->appendLine(L"Engine Results:");
		sb->appendLine();
		for (auto ok : resultTexts)
		{
			sb->appendLine(ok);
			sb->appendLine();
		}

		delete sb;
		return sb->toString();
	}

	void MyTaskResults::AddResultText(const std::wstring &resultsText)
	{
		resultTexts.push_back(resultsText);
	}

	void MyTaskResults::AddNiceText(const std::wstring &niceTextString)
	{
		niceText->appendLine(niceTextString);
	}
}
