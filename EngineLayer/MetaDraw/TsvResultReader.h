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
		static const std::wstring FullSequenceLabel;
		static const std::wstring Ms2ScanNumberLabel;
		static const std::wstring FilenameLabel;
		static const std::wstring TotalIonCurrentLabel;
		static const std::wstring PrecursorScanNumLabel;
		static const std::wstring PrecursorChargeLabel;
		static const std::wstring PrecursorMzLabel;
		static const std::wstring PrecursorMassLabel;
		static const std::wstring ScoreLabel;
		static const std::wstring DeltaScoreLabel;
		static const std::wstring NotchLabel;
		static const std::wstring BaseSeqLabel;
		static const std::wstring EssentialSeqLabel;
		static const std::wstring MissedCleavageLabel;
		static const std::wstring PeptideMonoMassLabel;
		static const std::wstring MassDiffDaLabel;
		static const std::wstring MassDiffPpmLabel;
		static const std::wstring ProteinAccessionLabel;
		static const std::wstring ProteinNameLabel;
		static const std::wstring GeneNameLabel;
		static const std::wstring OrganismNameLabel;
		static const std::wstring PeptideDesicriptionLabel;
		static const std::wstring StartAndEndResiduesInProteinLabel;
		static const std::wstring PreviousAminoAcidLabel;
		static const std::wstring NextAminoAcidLabel;
		static const std::wstring DecoyContamTargetLabel;
		static const std::wstring MatchedIonsLabel;
		static const std::wstring QValueLabel;
		static const std::wstring QValueNotchLabel;

		//Crosslinks
		static const std::wstring CrossTypeLabel;
		static const std::wstring LinkResiduesLabel;
		static const std::wstring ProteinLinkSiteLabel;
		static const std::wstring RankLabel;
		static const std::wstring BetaPeptideProteinAccessionLabel;
		static const std::wstring BetaPeptideProteinLinkSiteLabel;
		static const std::wstring BetaPeptideBaseSequenceLabel;
		static const std::wstring BetaPeptideFullSequenceLabel;
		static const std::wstring BetaPeptideTheoreticalMassLabel;
		static const std::wstring BetaPeptideScoreLabel;
		static const std::wstring BetaPeptideRankLabel;
		static const std::wstring BetaPeptideMatchedIonsLabel;
		static const std::wstring XLTotalScoreLabel;
		static const std::wstring ParentIonsLabel;

	private:
		static std::vector<wchar_t> const Split;

	public:
		static std::vector<MetaDrawPsm*> ReadTsv(const std::wstring &filePath, std::vector<std::wstring> &warnings);

	private:
		static std::unordered_map<std::wstring, int> ParseHeader(const std::wstring &header);
	};
}
