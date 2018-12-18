#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <cctype>
#include "stringhelper.h"


namespace MetaMorpheusGUI
{
	/// <summary>
	/// Provides filters and error handling for GUI forms
	/// </summary>
	class GlobalGuiSettings
	{
		/// <summary>
		/// Checks the validity of each setting passed from the GUI task windows
		/// </summary>
	public:
		static bool CheckTaskSettingsValidity(const std::wstring &precursorMassTolerance, const std::wstring &productMassTolerance, const std::wstring &maxMissedCleavages, const std::wstring &maxModificationIsoforms, const std::wstring &minPeptideLength, const std::wstring &maxPeptideLength, const std::wstring &maxThreads, const std::wstring &minScore, const std::wstring &peakFindingTolerance, const std::wstring &histogramBinWidth, const std::wstring &deconMaxAssumedCharge, const std::wstring &numPeaks, const std::wstring &minRatio, const std::wstring &numberOfDatabaseSearches, const std::wstring &maxModsPerPeptide, const std::wstring &maxFragmentMass, const std::wstring &qValueFilter);

		/// <summary>
		/// Checks to see if the given text contains non-numerical characters (letters, etc.)
		/// </summary>
		static bool CheckIsNumber(const std::wstring &text);

//		#region Check Task Validity

		static std::wstring MaxValueConversion(const std::wstring &text);

		static bool CheckDeconvolutionMaxAssumedChargeState(const std::wstring &text);

		static bool CheckTopNPeaks(const std::wstring &text);

		static bool CheckMinRatio(const std::wstring &text);

		static bool CheckPrecursorMassTolerance(const std::wstring &text);

		static bool CheckProductMassTolerance(const std::wstring &text);

		static bool CheckNumberOfDatabasePartitions(const std::wstring &text);

		static bool CheckMaxMissedCleavages(const std::wstring &text);

		static bool CheckMaxModificationIsoForms(const std::wstring &text);

		static bool CheckPeptideLength(const std::wstring &min, const std::wstring &max);

		static bool CheckMaxModsPerPeptide(const std::wstring &text);

		static bool CheckMaxFragementMass(const std::wstring &text);

		static bool CheckMaxThreads(const std::wstring &text);

		static bool CheckMinScoreAllowed(const std::wstring &text);

		static bool CheckPeakFindingTolerance(const std::wstring &text);

		static bool CheckHistogramBinWidth(const std::wstring &text);

		static bool CheckQValueFilter(const std::wstring &text);

		static bool VariableModCheck(std::vector<(std::wstring, std::wstring)*> &listOfModsVariable);

//		#endregion
	};
}
