#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MetaDrawPsm; }


namespace EngineLayer
{
	class TsvResultReader
	{
	public:
		static const std::string FullSequenceLabel;
		static const std::string Ms2ScanNumberLabel;
		static const std::string FilenameLabel;
		static const std::string TotalIonCurrentLabel;
		static const std::string PrecursorScanNumLabel;
		static const std::string PrecursorChargeLabel;
		static const std::string PrecursorMzLabel;
		static const std::string PrecursorMassLabel;
		static const std::string ScoreLabel;
		static const std::string DeltaScoreLabel;
		static const std::string NotchLabel;
		static const std::string BaseSeqLabel;
		static const std::string EssentialSeqLabel;
		static const std::string MissedCleavageLabel;
		static const std::string PeptideMonoMassLabel;
		static const std::string MassDiffDaLabel;
		static const std::string MassDiffPpmLabel;
		static const std::string ProteinAccessionLabel;
		static const std::string ProteinNameLabel;
		static const std::string GeneNameLabel;
		static const std::string OrganismNameLabel;
		static const std::string PeptideDesicriptionLabel;
		static const std::string StartAndEndResiduesInProteinLabel;
		static const std::string PreviousAminoAcidLabel;
		static const std::string NextAminoAcidLabel;
		static const std::string DecoyContamTargetLabel;
		static const std::string MatchedIonsLabel;
		static const std::string QValueLabel;
		static const std::string QValueNotchLabel;

		//Crosslinks
		static const std::string CrossTypeLabel;
		static const std::string LinkResiduesLabel;
		static const std::string ProteinLinkSiteLabel;
		static const std::string RankLabel;
		static const std::string BetaPeptideProteinAccessionLabel;
		static const std::string BetaPeptideProteinLinkSiteLabel;
		static const std::string BetaPeptideBaseSequenceLabel;
		static const std::string BetaPeptideFullSequenceLabel;
		static const std::string BetaPeptideTheoreticalMassLabel;
		static const std::string BetaPeptideScoreLabel;
		static const std::string BetaPeptideRankLabel;
		static const std::string BetaPeptideMatchedIonsLabel;
		static const std::string XLTotalScoreLabel;
		static const std::string ParentIonsLabel;

	private:
		static std::vector<char> const Split;

	public:
		static std::vector<MetaDrawPsm*> ReadTsv(const std::string &filePath, std::vector<std::string> &warnings);

	private:
		static std::unordered_map<std::string, int> ParseHeader(const std::string &header);
	};
}
