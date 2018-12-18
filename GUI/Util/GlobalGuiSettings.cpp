#include "GlobalGuiSettings.h"


namespace MetaMorpheusGUI
{

	bool GlobalGuiSettings::CheckTaskSettingsValidity(const std::wstring &precursorMassTolerance, const std::wstring &productMassTolerance, const std::wstring &maxMissedCleavages, const std::wstring &maxModificationIsoforms, const std::wstring &minPeptideLength, const std::wstring &maxPeptideLength, const std::wstring &maxThreads, const std::wstring &minScore, const std::wstring &peakFindingTolerance, const std::wstring &histogramBinWidth, const std::wstring &deconMaxAssumedCharge, const std::wstring &numPeaks, const std::wstring &minRatio, const std::wstring &numberOfDatabaseSearches, const std::wstring &maxModsPerPeptide, const std::wstring &maxFragmentMass, const std::wstring &qValueFilter)
	{
		maxMissedCleavages = MaxValueConversion(maxMissedCleavages);
		maxPeptideLength = MaxValueConversion(maxPeptideLength);

		std::vector<bool> results;
		results.push_back((CheckPrecursorMassTolerance(precursorMassTolerance)));
		results.push_back((CheckProductMassTolerance(productMassTolerance)));
		results.push_back((CheckMaxMissedCleavages(maxMissedCleavages)));
		results.push_back((CheckMaxModificationIsoForms(maxModificationIsoforms)));
		results.push_back((CheckPeptideLength(minPeptideLength, maxPeptideLength)));
		results.push_back((CheckMaxThreads(maxThreads)));
		results.push_back((CheckMinScoreAllowed(minScore)));
		results.push_back((CheckPeakFindingTolerance(peakFindingTolerance)));
		results.push_back((CheckHistogramBinWidth(histogramBinWidth)));
		results.push_back((CheckDeconvolutionMaxAssumedChargeState(deconMaxAssumedCharge)));
		results.push_back((CheckTopNPeaks(numPeaks)));
		results.push_back((CheckMinRatio(minRatio)));
		results.push_back((CheckNumberOfDatabasePartitions(numberOfDatabaseSearches)));
		results.push_back((CheckMaxModsPerPeptide(maxModsPerPeptide)));
		results.push_back((CheckMaxFragementMass(maxFragmentMass)));
		results.push_back((CheckQValueFilter(qValueFilter)));

		if (std::find(results.begin(), results.end(), false) != results.end())
		{
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckIsNumber(const std::wstring &text)
	{
		bool result = true;
		for (auto character : text)
		{
			if (!std::isdigit(character) && !(character == L'.') && !(character == L'-'))
			{
				result = false;
			}
		}
		return result;
	}

	std::wstring GlobalGuiSettings::MaxValueConversion(const std::wstring &text)
	{
		if (text.empty())
		{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			text = std::numeric_limits<int>::max().ToString();
		}
		return text;
	}

	bool GlobalGuiSettings::CheckDeconvolutionMaxAssumedChargeState(const std::wstring &text)
	{
		int deconMaxAssumedCharge;
		if (!int::TryParse(text, deconMaxAssumedCharge) || deconMaxAssumedCharge < 1)
		{
			MessageBox::Show(L"The maximum assumed charge state for deconvolution is invalid. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"') + L"\n Please enter a positive, non-zero number.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckTopNPeaks(const std::wstring &text)
	{
		if (text.length() == 0)
		{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			text = std::numeric_limits<int>::max().ToString();
		}

		int numPeaks;
		if (!int::TryParse(text, numPeaks) || numPeaks < 1)
		{
			MessageBox::Show(L"The Top N Peaks to be retained must be greater than zero. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"') + L"\n Please enter a positive number.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckMinRatio(const std::wstring &text)
	{
		double minRatio;
		if (!double::TryParse(text, minRatio) || minRatio < 0 || minRatio > 1)
		{
			MessageBox::Show(L"The minimum intensity ratio must be between zero and one. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"'));
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckPrecursorMassTolerance(const std::wstring &text)
	{

		double precursorMassTolerance;
		if (!double::TryParse(text, precursorMassTolerance) || precursorMassTolerance <= 0)
		{
			MessageBox::Show(L"The precursor mass tolerance is invalid. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"') + L"\n Please enter a positive number.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckProductMassTolerance(const std::wstring &text)
	{
		double productMassTolerance;
		if (!double::TryParse(text, productMassTolerance) || productMassTolerance <= 0)
		{
			MessageBox::Show(L"The product mass tolerance is invalid. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"') + L"\n Please enter a positive number.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckNumberOfDatabasePartitions(const std::wstring &text)
	{
		int numberOfDatabaseSearches;
		if (!int::TryParse(text, numberOfDatabaseSearches) || numberOfDatabaseSearches <= 0)
		{
			MessageBox::Show(L"The number of database partitions is invalid. At least one database is required for searching.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckMaxMissedCleavages(const std::wstring &text)
	{
		int maxMissedCleavages;
		if (!int::TryParse(text, maxMissedCleavages) || maxMissedCleavages < 0)
		{
			MessageBox::Show(L"The number of missed cleavages is invalid. Please enter an integer zero or greater.");
			return false;
		}

		return true;
	}

	bool GlobalGuiSettings::CheckMaxModificationIsoForms(const std::wstring &text)
	{
		int maxModificationIsoforms;
		if (!int::TryParse(text, maxModificationIsoforms) || maxModificationIsoforms < 1)
		{
			MessageBox::Show(L"The maximum number of modification isoforms is invalid. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"') + L"\n Please enter a positive, non-zero number.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckPeptideLength(const std::wstring &min, const std::wstring &max)
	{
		int minPeptideLength;
		if (!int::TryParse(min, minPeptideLength) || minPeptideLength < 1)
		{
			MessageBox::Show(L"The minimum peptide length must be a positive integer");
			return false;
		}

		int maxPeptideLength;
		if (!int::TryParse(max, maxPeptideLength) || maxPeptideLength < 1)
		{
			MessageBox::Show(L"The maximum peptide length must be a positive integer");
			return false;
		}

		if (std::stoi(min) > std::stoi(max))
		{
			MessageBox::Show(L"The maximum peptide length must be greater than or equal to the minimum peptide length.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckMaxModsPerPeptide(const std::wstring &text)
	{
		int maxModsPerPeptide;
		if (!int::TryParse(text, maxModsPerPeptide) || maxModsPerPeptide < 0)
		{
			MessageBox::Show(L"The max mods per peptide allowed is invalid. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"') + L"\n Please enter a number greater than or equal to zero.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckMaxFragementMass(const std::wstring &text)
	{
		int maxFragmentMass;
		if (!int::TryParse(text, maxFragmentMass) || maxFragmentMass < 0)
		{
			MessageBox::Show(L"The max fragment mass is invalid. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"') + L"\n Please enter a positive, non-zero number.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckMaxThreads(const std::wstring &text)
	{
		int maxThreads;
		if (!int::TryParse(text, maxThreads) || maxThreads > Environment::ProcessorCount || maxThreads < 1)
		{
			MessageBox::Show(L"Your current device has " + std::to_wstring(Environment::ProcessorCount) + L" processors. \n Please select a positive value less than or equal to this number.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckMinScoreAllowed(const std::wstring &text)
	{
		double minScore;
		if (!double::TryParse(text, minScore) || minScore < 1)
		{
			MessageBox::Show(L"The minimum score allowed is invalid. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"') + L"\n Please enter a positive, non-zero number.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckPeakFindingTolerance(const std::wstring &text)
	{
		double peakFindingTolerance;
		if (!double::TryParse(text, peakFindingTolerance) || peakFindingTolerance <= 0)
		{
			MessageBox::Show(L"The peak finding tolerance is invalid. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"') + L"\n Please enter a positive number.");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckHistogramBinWidth(const std::wstring &text)
	{
		float binWidth;
		if (!float::TryParse(text, binWidth) || binWidth < 0 || binWidth > 1)
		{
			MessageBox::Show(L"The histogram bin width must be between zero and one Daltons. \n You entered " + StringHelper::toString(L'"') + text + StringHelper::toString(L'"'));
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::CheckQValueFilter(const std::wstring &text)
	{
		double qValue;
		if (!double::TryParse(text, qValue) || qValue < 0 || qValue > 1)
		{
			MessageBox::Show(L"The q-value cutoff must be a number between 0 and 1");
			return false;
		}
		return true;
	}

	bool GlobalGuiSettings::VariableModCheck(std::vector<(std::wstring, std::wstring)*> &listOfModsVariable)
	{
		if (listOfModsVariable.size() > 1)
		{
			auto dialogResult = MessageBox::Show(L"More than 1 modification has been selected as variable. Using the GPTMD task to discover modifications is recommended instead. \n\nContinue anyway?", L"Multiple Variable Mods Detected", MessageBoxButton::OKCancel);
			if (dialogResult == MessageBoxResult::Cancel)
			{
				return false;
			}
		}
		return true;
	}
}
