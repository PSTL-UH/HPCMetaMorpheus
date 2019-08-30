#pragma once

#include <string>
#include <vector>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
//namespace TaskLayer { class DbForTask; }
#include "DbForTask.h"

//namespace TaskLayer { class MetaMorpheusTask; }
#include "MetaMorpheusTask.h"


namespace TaskLayer
{
	class MyTaskResults
	{
	public:
		std::vector<std::string> NewSpectra; // calibration writes new calibrated spectra
		std::vector<DbForTask*> NewDatabases; // gptmd writes new annotated databases
		std::vector<std::string> NewFileSpecificTomls; // calibration writes suggested ppm tolerances
		TimeSpan Time;

	private:
		const std::vector<std::string> resultTexts;

		StringBuilder *const niceText = new StringBuilder();

	public:
		virtual ~MyTaskResults()
		{
			delete niceText;
		}

		MyTaskResults(MetaMorpheusTask *s);

		std::string ToString();

		void AddResultText(const std::string &resultsText);

		void AddNiceText(const std::string &niceTextString);
	};
}
