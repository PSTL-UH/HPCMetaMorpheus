#pragma once

#include <string>
#include <vector>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class DbForTask; }
namespace TaskLayer { class MetaMorpheusTask; }


namespace TaskLayer
{
	class MyTaskResults
	{
	public:
		std::vector<std::wstring> NewSpectra; // calibration writes new calibrated spectra
		std::vector<DbForTask*> NewDatabases; // gptmd writes new annotated databases
		std::vector<std::wstring> NewFileSpecificTomls; // calibration writes suggested ppm tolerances
		TimeSpan Time;

	private:
		const std::vector<std::wstring> resultTexts;

		StringBuilder *const niceText = new StringBuilder();

	public:
		virtual ~MyTaskResults()
		{
			delete niceText;
		}

		MyTaskResults(MetaMorpheusTask *s);

		std::wstring ToString() override;

		void AddResultText(const std::wstring &resultsText);

		void AddNiceText(const std::wstring &niceTextString);
	};
}
