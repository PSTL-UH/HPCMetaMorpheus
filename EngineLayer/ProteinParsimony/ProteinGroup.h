#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include "stringhelper.h"
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }

using namespace FlashLFQ;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	class ProteinGroup
	{
	private:
		bool privateIsDecoy = false;
		bool privateIsContaminant = false;
		std::vector<SpectraFileInfo*> privateFilesForQuantification;
		std::unordered_set<Protein*> privateProteins;
		std::wstring privateProteinGroupName;
		double privateProteinGroupScore = 0;
		std::unordered_set<PeptideWithSetModifications*> privateAllPeptides;
		std::unordered_set<PeptideWithSetModifications*> privateUniquePeptides;
		std::unordered_set<PeptideSpectralMatch*> privateAllPsmsBelowOnePercentFDR;
		std::vector<double> privateSequenceCoveragePercent;
		std::vector<std::wstring> privateSequenceCoverageDisplayList;
		std::vector<std::wstring> privateSequenceCoverageDisplayListWithMods;
		double privateQValue = 0;
		double privateBestPeptideQValue = 0;
		double privateBestPeptideScore = 0;
		int privateCumulativeTarget = 0;
		int privateCumulativeDecoy = 0;
		bool privateDisplayModsOnPeptides = false;
		std::vector<std::wstring> privateModsInfo;
		std::unordered_map<SpectraFileInfo*, double> privateIntensitiesByFile;

	public:
		ProteinGroup(std::unordered_set<Protein*> &proteins, std::unordered_set<PeptideWithSetModifications*> &peptides, std::unordered_set<PeptideWithSetModifications*> &uniquePeptides);

		bool getIsDecoy() const;

		bool getIsContaminant() const;

		std::vector<SpectraFileInfo*> getFilesForQuantification() const;
		void setFilesForQuantification(const std::vector<SpectraFileInfo*> &value);

		std::unordered_set<Protein*> getProteins() const;
		void setProteins(const std::unordered_set<Protein*> &value);

		std::wstring getProteinGroupName() const;
		void setProteinGroupName(const std::wstring &value);

		double getProteinGroupScore() const;
		void setProteinGroupScore(double value);

		std::unordered_set<PeptideWithSetModifications*> getAllPeptides() const;
		void setAllPeptides(const std::unordered_set<PeptideWithSetModifications*> &value);

		std::unordered_set<PeptideWithSetModifications*> getUniquePeptides() const;
		void setUniquePeptides(const std::unordered_set<PeptideWithSetModifications*> &value);

		std::unordered_set<PeptideSpectralMatch*> getAllPsmsBelowOnePercentFDR() const;
		void setAllPsmsBelowOnePercentFDR(const std::unordered_set<PeptideSpectralMatch*> &value);

		std::vector<double> getSequenceCoveragePercent() const;
		void setSequenceCoveragePercent(const std::vector<double> &value);

		std::vector<std::wstring> getSequenceCoverageDisplayList() const;
		void setSequenceCoverageDisplayList(const std::vector<std::wstring> &value);

		std::vector<std::wstring> getSequenceCoverageDisplayListWithMods() const;
		void setSequenceCoverageDisplayListWithMods(const std::vector<std::wstring> &value);

		double getQValue() const;
		void setQValue(double value);

		double getBestPeptideQValue() const;
		void setBestPeptideQValue(double value);

		double getBestPeptideScore() const;
		void setBestPeptideScore(double value);

		int getCumulativeTarget() const;
		void setCumulativeTarget(int value);

		int getCumulativeDecoy() const;
		void setCumulativeDecoy(int value);

		bool getDisplayModsOnPeptides() const;
		void setDisplayModsOnPeptides(bool value);

		std::vector<std::wstring> getModsInfo() const;
		void setModsInfo(const std::vector<std::wstring> &value);

		std::unordered_map<SpectraFileInfo*, double> getIntensitiesByFile() const;
		void setIntensitiesByFile(const std::unordered_map<SpectraFileInfo*, double> &value);

	private:
		std::vector<Protein*> ListOfProteinsOrderedByAccession;

	public:
		std::wstring GetTabSeparatedHeader();

		std::wstring ToString() override;

		// this method is only used internally, to make protein grouping faster
		// this is NOT an output and is NOT used for protein FDR calculations
		void Score();

		void CalculateSequenceCoverage();

		void MergeProteinGroupWith(ProteinGroup *other);

		ProteinGroup *ConstructSubsetProteinGroup(const std::wstring &fullFilePath);
	};
}
