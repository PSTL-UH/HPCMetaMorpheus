#include "MyTaskResults.h"
#include "DbForTask.h"
#include "MetaMorpheusTask.h"


namespace TaskLayer
{

	MyTaskResults::MyTaskResults(MetaMorpheusTask *s) : resultTexts(std::vector<std::string>())
	{
	}

	std::string MyTaskResults::ToString()
	{
		StringBuilder *sb = new StringBuilder();
		sb->appendLine("Time to run task: " + Time);
		sb->appendLine();
		sb->appendLine();
		sb->appendLine("--------------------------------------------------");
		if ((NewSpectra.size() > 0 && NewSpectra.Any()) || (NewDatabases.size() > 0 && NewDatabases.Any()))
		{
			sb->appendLine();
			sb->appendLine();
			sb->appendLine("New files:");
			if (NewSpectra.size() > 0 && NewSpectra.Any())
			{
				sb->appendLine("New spectra: ");
				sb->appendLine();
				sb->appendLine(std::string::Join("\r\n" + "\t", NewSpectra));
			}
			if (NewDatabases.size() > 0 && NewDatabases.Any())
			{
				sb->appendLine("New databases: ");
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				sb->appendLine(std::string::Join("\r\n" + "\t", NewDatabases.Select([&] (std::any b)
				{
					b::FilePath;
				})).ToString());
			}
			sb->appendLine();
			sb->appendLine();
			sb->appendLine("--------------------------------------------------");
		}
		sb->appendLine();
		sb->appendLine();
		sb->appendLine(niceText->toString());
		sb->appendLine();
		sb->appendLine();
		sb->appendLine("--------------------------------------------------");
		sb->appendLine();
		sb->appendLine();
		sb->appendLine("Engine Results:");
		sb->appendLine();
		for (auto ok : resultTexts)
		{
			sb->appendLine(ok);
			sb->appendLine();
		}

		delete sb;
		return sb->toString();
	}

	void MyTaskResults::AddResultText(const std::string &resultsText)
	{
		resultTexts.push_back(resultsText);
	}

	void MyTaskResults::AddNiceText(const std::string &niceTextString)
	{
		niceText->appendLine(niceTextString);
	}
}
