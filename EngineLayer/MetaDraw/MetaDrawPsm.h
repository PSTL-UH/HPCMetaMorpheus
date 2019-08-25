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
		std::string privateFullSequence;
		int privateMs2ScanNumber = 0;
		std::string privateFilename;
		int privatePrecursorScanNum = 0;
		int privatePrecursorCharge = 0;
		double privatePrecursorMz = 0;
		double privatePrecursorMass = 0;
		double privateScore = 0;
		std::string privateProteinAccession;
		std::vector<MatchedFragmentIon*> privateMatchedIons;
		double privateQValue = 0;
		std::optional<double> privateTotalIonCurrent;
		std::optional<double> privateDeltaScore;
		std::string privateNotch;
		std::string privateBaseSeq;
		std::string privateEssentialSeq;
		std::string privateMissedCleavage;
		std::string privatePeptideMonoMass;
		std::string privateMassDiffDa;
		std::string privateMassDiffPpm;
		std::string privateProteinName;
		std::string privateGeneName;
		std::string privateOrganismName;
		std::string privatePeptideDesicription;
		std::string privateStartAndEndResiduesInProtein;
		std::string privatePreviousAminoAcid;
		std::string privateNextAminoAcid;
		std::string privateDecoyContamTarget;
		std::optional<double> privateQValueNotch;
		std::string privateCrossType;
		std::string privateLinkResidues;
		std::optional<int> privateProteinLinkSite;
		std::optional<int> privateRank;
		std::string privateBetaPeptideProteinAccession;
		std::optional<int> privateBetaPeptideProteinLinkSite;
		std::string privateBetaPeptideBaseSequence;
		std::string privateBetaPeptideFullSequence;
		std::string privateBetaPeptideTheoreticalMass;
		std::optional<double> privateBetaPeptideScore;
		std::optional<int> privateBetaPeptideRank;
		std::vector<MatchedFragmentIon*> privateBetaPeptideMatchedIons;
		std::optional<double> privateXLTotalScore;
		std::string privateParentIons;

		static Regex *const IonParser;
		static std::vector<char> const MzSplit;

	public:
		std::string getFullSequence() const;
		int getMs2ScanNumber() const;
		std::string getFilename() const;
		int getPrecursorScanNum() const;
		int getPrecursorCharge() const;
		double getPrecursorMz() const;
		double getPrecursorMass() const;
		double getScore() const;
		std::string getProteinAccession() const;
		std::vector<MatchedFragmentIon*> getMatchedIons() const;
		double getQValue() const;

		std::optional<double> getTotalIonCurrent() const;
		std::optional<double> getDeltaScore() const;
		std::string getNotch() const;
		std::string getBaseSeq() const;
		std::string getEssentialSeq() const;
		std::string getMissedCleavage() const;
		std::string getPeptideMonoMass() const;
		std::string getMassDiffDa() const;
		std::string getMassDiffPpm() const;
		std::string getProteinName() const;
		std::string getGeneName() const;
		std::string getOrganismName() const;
		std::string getPeptideDesicription() const;
		std::string getStartAndEndResiduesInProtein() const;
		std::string getPreviousAminoAcid() const;
		std::string getNextAminoAcid() const;
		std::string getDecoyContamTarget() const;
		std::optional<double> getQValueNotch() const;

		//For crosslink
		std::string getCrossType() const;
		std::string getLinkResidues() const;
		std::optional<int> getProteinLinkSite() const;
		std::optional<int> getRank() const;
		std::string getBetaPeptideProteinAccession() const;
		std::optional<int> getBetaPeptideProteinLinkSite() const;
		std::string getBetaPeptideBaseSequence() const;
		std::string getBetaPeptideFullSequence() const;
		std::string getBetaPeptideTheoreticalMass() const;
		std::optional<double> getBetaPeptideScore() const;
		std::optional<int> getBetaPeptideRank() const;
		std::vector<MatchedFragmentIon*> getBetaPeptideMatchedIons() const;
		std::optional<double> getXLTotalScore() const;
		std::string getParentIons() const;

		MetaDrawPsm(const std::string &line, std::vector<char> &split, std::unordered_map<std::string, int> &parsedHeader);

	private:
		static std::vector<MatchedFragmentIon*> ReadFragmentIonsFromString(const std::string &matchedMzString, const std::string &peptideBaseSequence);
	};
}
