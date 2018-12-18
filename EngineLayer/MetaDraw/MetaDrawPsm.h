#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <any>
#include <optional>
#include "stringhelper.h"

using namespace Proteomics::Fragmentation;
using namespace Chemistry;

namespace EngineLayer
{
	class MetaDrawPsm
	{
	private:
		std::wstring privateFullSequence;
		int privateMs2ScanNumber = 0;
		std::wstring privateFilename;
		int privatePrecursorScanNum = 0;
		int privatePrecursorCharge = 0;
		double privatePrecursorMz = 0;
		double privatePrecursorMass = 0;
		double privateScore = 0;
		std::wstring privateProteinAccession;
		std::vector<MatchedFragmentIon*> privateMatchedIons;
		double privateQValue = 0;
		std::optional<double> privateTotalIonCurrent;
		std::optional<double> privateDeltaScore;
		std::wstring privateNotch;
		std::wstring privateBaseSeq;
		std::wstring privateEssentialSeq;
		std::wstring privateMissedCleavage;
		std::wstring privatePeptideMonoMass;
		std::wstring privateMassDiffDa;
		std::wstring privateMassDiffPpm;
		std::wstring privateProteinName;
		std::wstring privateGeneName;
		std::wstring privateOrganismName;
		std::wstring privatePeptideDesicription;
		std::wstring privateStartAndEndResiduesInProtein;
		std::wstring privatePreviousAminoAcid;
		std::wstring privateNextAminoAcid;
		std::wstring privateDecoyContamTarget;
		std::optional<double> privateQValueNotch;
		std::wstring privateCrossType;
		std::wstring privateLinkResidues;
		std::optional<int> privateProteinLinkSite;
		std::optional<int> privateRank;
		std::wstring privateBetaPeptideProteinAccession;
		std::optional<int> privateBetaPeptideProteinLinkSite;
		std::wstring privateBetaPeptideBaseSequence;
		std::wstring privateBetaPeptideFullSequence;
		std::wstring privateBetaPeptideTheoreticalMass;
		std::optional<double> privateBetaPeptideScore;
		std::optional<int> privateBetaPeptideRank;
		std::vector<MatchedFragmentIon*> privateBetaPeptideMatchedIons;
		std::optional<double> privateXLTotalScore;
		std::wstring privateParentIons;

		static Regex *const IonParser;
		static std::vector<wchar_t> const MzSplit;

	public:
		std::wstring getFullSequence() const;
		int getMs2ScanNumber() const;
		std::wstring getFilename() const;
		int getPrecursorScanNum() const;
		int getPrecursorCharge() const;
		double getPrecursorMz() const;
		double getPrecursorMass() const;
		double getScore() const;
		std::wstring getProteinAccession() const;
		std::vector<MatchedFragmentIon*> getMatchedIons() const;
		double getQValue() const;

		std::optional<double> getTotalIonCurrent() const;
		std::optional<double> getDeltaScore() const;
		std::wstring getNotch() const;
		std::wstring getBaseSeq() const;
		std::wstring getEssentialSeq() const;
		std::wstring getMissedCleavage() const;
		std::wstring getPeptideMonoMass() const;
		std::wstring getMassDiffDa() const;
		std::wstring getMassDiffPpm() const;
		std::wstring getProteinName() const;
		std::wstring getGeneName() const;
		std::wstring getOrganismName() const;
		std::wstring getPeptideDesicription() const;
		std::wstring getStartAndEndResiduesInProtein() const;
		std::wstring getPreviousAminoAcid() const;
		std::wstring getNextAminoAcid() const;
		std::wstring getDecoyContamTarget() const;
		std::optional<double> getQValueNotch() const;

		//For crosslink
		std::wstring getCrossType() const;
		std::wstring getLinkResidues() const;
		std::optional<int> getProteinLinkSite() const;
		std::optional<int> getRank() const;
		std::wstring getBetaPeptideProteinAccession() const;
		std::optional<int> getBetaPeptideProteinLinkSite() const;
		std::wstring getBetaPeptideBaseSequence() const;
		std::wstring getBetaPeptideFullSequence() const;
		std::wstring getBetaPeptideTheoreticalMass() const;
		std::optional<double> getBetaPeptideScore() const;
		std::optional<int> getBetaPeptideRank() const;
		std::vector<MatchedFragmentIon*> getBetaPeptideMatchedIons() const;
		std::optional<double> getXLTotalScore() const;
		std::wstring getParentIons() const;

		MetaDrawPsm(const std::wstring &line, std::vector<wchar_t> &split, std::unordered_map<std::wstring, int> &parsedHeader);

	private:
		static std::vector<MatchedFragmentIon*> ReadFragmentIonsFromString(const std::wstring &matchedMzString, const std::wstring &peptideBaseSequence);
	};
}
