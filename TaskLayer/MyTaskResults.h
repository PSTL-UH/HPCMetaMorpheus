#pragma once

#include <string>
#include <vector>
#include "stringbuilder.h"

#include "DbForTask.h"

// need to use forward declaration instead of including MetaMorpheusTask.h due to circular dependence
namespace TaskLayer { class MetaMorpheusTask; }



namespace TaskLayer
{
	class MyTaskResults
	{
	public:
		std::vector<std::string> NewSpectra; // calibration writes new calibrated spectra
		std::vector<DbForTask*> NewDatabases; // gptmd writes new annotated databases
		std::vector<std::string> NewFileSpecificTomls; // calibration writes suggested ppm tolerances
		//TimeSpan Time;
                double Time;
	private:
		std::vector<std::string> resultTexts;

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
