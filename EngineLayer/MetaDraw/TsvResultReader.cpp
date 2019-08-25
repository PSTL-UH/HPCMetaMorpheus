#include "TsvResultReader.h"
#include "MetaDrawPsm.h"
#include "../MetaMorpheusException.h"


namespace EngineLayer
{

const std::string TsvResultReader::FullSequenceLabel = "Full Sequence";
const std::string TsvResultReader::Ms2ScanNumberLabel = "Scan Number";
const std::string TsvResultReader::FilenameLabel = "File Name";
const std::string TsvResultReader::TotalIonCurrentLabel = "Total Ion Current";
const std::string TsvResultReader::PrecursorScanNumLabel = "Precursor Scan Number";
const std::string TsvResultReader::PrecursorChargeLabel = "Precursor Charge";
const std::string TsvResultReader::PrecursorMzLabel = "Precursor MZ";
const std::string TsvResultReader::PrecursorMassLabel = "Precursor Mass";
const std::string TsvResultReader::ScoreLabel = "Score";
const std::string TsvResultReader::DeltaScoreLabel = "Delta Score";
const std::string TsvResultReader::NotchLabel = "Notch";
const std::string TsvResultReader::BaseSeqLabel = "Base Sequence";
const std::string TsvResultReader::EssentialSeqLabel = "Essential Sequence";
const std::string TsvResultReader::MissedCleavageLabel = "Missed Cleavages";
const std::string TsvResultReader::PeptideMonoMassLabel = "Peptide Monoisotopic Mass";
const std::string TsvResultReader::MassDiffDaLabel = "Mass Diff (Da)";
const std::string TsvResultReader::MassDiffPpmLabel = "Mass Diff (ppm)";
const std::string TsvResultReader::ProteinAccessionLabel = "Protein Accession";
const std::string TsvResultReader::ProteinNameLabel = "Protein Name";
const std::string TsvResultReader::GeneNameLabel = "Gene Name";
const std::string TsvResultReader::OrganismNameLabel = "Organism Name";
const std::string TsvResultReader::PeptideDesicriptionLabel = "Peptide Description";
const std::string TsvResultReader::StartAndEndResiduesInProteinLabel = "Start and End Residues In Protein";
const std::string TsvResultReader::PreviousAminoAcidLabel = "Previous Amino Acid";
const std::string TsvResultReader::NextAminoAcidLabel = "Next Amino Acid";
const std::string TsvResultReader::DecoyContamTargetLabel = "Decoy/Contaminant/Target";
const std::string TsvResultReader::MatchedIonsLabel = "Matched Ion Mass-To-Charge Ratios";
const std::string TsvResultReader::QValueLabel = "QValue";
const std::string TsvResultReader::QValueNotchLabel = "QValue Notch";
const std::string TsvResultReader::CrossTypeLabel = "Cross Type";
const std::string TsvResultReader::LinkResiduesLabel = "Link Residues";
const std::string TsvResultReader::ProteinLinkSiteLabel = "Protein Link Site";
const std::string TsvResultReader::RankLabel = "Rank";
const std::string TsvResultReader::BetaPeptideProteinAccessionLabel = "Beta Peptide Protein Accession";
const std::string TsvResultReader::BetaPeptideProteinLinkSiteLabel = "Beta Peptide Protein LinkSite";
const std::string TsvResultReader::BetaPeptideBaseSequenceLabel = "Beta Peptide Base Sequence";
const std::string TsvResultReader::BetaPeptideFullSequenceLabel = "Beta Peptide Full Sequence";
const std::string TsvResultReader::BetaPeptideTheoreticalMassLabel = "Beta Peptide Theoretical Mass";
const std::string TsvResultReader::BetaPeptideScoreLabel = "Beta Peptide Score";
const std::string TsvResultReader::BetaPeptideRankLabel = "Beta Peptide Rank";
const std::string TsvResultReader::BetaPeptideMatchedIonsLabel = "Beta Peptide Matched Ion Mass-To-Charge Ratios";
const std::string TsvResultReader::XLTotalScoreLabel = "XL Total Score";
const std::string TsvResultReader::ParentIonsLabel = "Parent Ions";
std::vector<char> const TsvResultReader::Split = {'\t'};

	std::vector<MetaDrawPsm*> TsvResultReader::ReadTsv(const std::string &filePath, std::vector<std::string> &warnings)
	{
		std::vector<MetaDrawPsm*> psms;
		warnings = std::vector<std::string>();

		StreamReader *reader = nullptr;
		try
		{
			reader = new StreamReader(filePath);
		}
		catch (const std::runtime_error &e)
		{
		   delete reader;
		   throw MetaMorpheusException("Could not read file: " + e.what());
		}

		int lineCount = 0;

		std::string line;
		std::unordered_map<std::string, int> parsedHeader;

		while (reader->Peek() > 0)
		{
			lineCount++;

			line = reader->ReadLine();

			if (lineCount == 1)
			{
				parsedHeader = ParseHeader(line);
				continue;
			}

			try
			{
				MetaDrawPsm tempVar(line, Split, parsedHeader);
				psms.push_back(&tempVar);
			}
			catch (const std::runtime_error &e1)
			{
				warnings.push_back("Could not read line: " + std::to_string(lineCount));
			}
		}

		reader->Close();

		if ((lineCount - 1) != psms.size())
		{
			warnings.push_back("Warning: " + std::to_string((lineCount - 1) - psms.size()) + " PSMs were not read.");
		}

		delete reader;
		return psms;
	}

	std::unordered_map<std::string, int> TsvResultReader::ParseHeader(const std::string &header)
	{
		auto parsedHeader = std::unordered_map<std::string, int>();
		auto spl = header.Split(Split);

		parsedHeader.emplace(FullSequenceLabel, Array::IndexOf(spl, FullSequenceLabel));
		parsedHeader.emplace(Ms2ScanNumberLabel, Array::IndexOf(spl, Ms2ScanNumberLabel));
		parsedHeader.emplace(FilenameLabel, Array::IndexOf(spl, FilenameLabel));
		parsedHeader.emplace(TotalIonCurrentLabel, Array::IndexOf(spl, TotalIonCurrentLabel));
		parsedHeader.emplace(PrecursorScanNumLabel, Array::IndexOf(spl, PrecursorScanNumLabel));
		parsedHeader.emplace(PrecursorChargeLabel, Array::IndexOf(spl, PrecursorChargeLabel));
		parsedHeader.emplace(PrecursorMzLabel, Array::IndexOf(spl, PrecursorMzLabel));
		parsedHeader.emplace(PrecursorMassLabel, Array::IndexOf(spl, PrecursorMassLabel));
		parsedHeader.emplace(ScoreLabel, Array::IndexOf(spl, ScoreLabel));
		parsedHeader.emplace(DeltaScoreLabel, Array::IndexOf(spl, DeltaScoreLabel));
		parsedHeader.emplace(NotchLabel, Array::IndexOf(spl, NotchLabel));
		parsedHeader.emplace(BaseSeqLabel, Array::IndexOf(spl, BaseSeqLabel));
		parsedHeader.emplace(EssentialSeqLabel, Array::IndexOf(spl, EssentialSeqLabel));
		parsedHeader.emplace(MissedCleavageLabel, Array::IndexOf(spl, MissedCleavageLabel));
		parsedHeader.emplace(PeptideMonoMassLabel, Array::IndexOf(spl, PeptideMonoMassLabel));
		parsedHeader.emplace(MassDiffDaLabel, Array::IndexOf(spl, MassDiffDaLabel));
		parsedHeader.emplace(MassDiffPpmLabel, Array::IndexOf(spl, MassDiffPpmLabel));
		parsedHeader.emplace(ProteinAccessionLabel, Array::IndexOf(spl, ProteinAccessionLabel));
		parsedHeader.emplace(ProteinNameLabel, Array::IndexOf(spl, ProteinNameLabel));
		parsedHeader.emplace(GeneNameLabel, Array::IndexOf(spl, GeneNameLabel));
		parsedHeader.emplace(OrganismNameLabel, Array::IndexOf(spl, OrganismNameLabel));
		parsedHeader.emplace(PeptideDesicriptionLabel, Array::IndexOf(spl, PeptideDesicriptionLabel));
		parsedHeader.emplace(StartAndEndResiduesInProteinLabel, Array::IndexOf(spl, StartAndEndResiduesInProteinLabel));
		parsedHeader.emplace(PreviousAminoAcidLabel, Array::IndexOf(spl, PreviousAminoAcidLabel));
		parsedHeader.emplace(NextAminoAcidLabel, Array::IndexOf(spl, NextAminoAcidLabel));
		parsedHeader.emplace(DecoyContamTargetLabel, Array::IndexOf(spl, DecoyContamTargetLabel));
		parsedHeader.emplace(MatchedIonsLabel, Array::IndexOf(spl, MatchedIonsLabel));
		parsedHeader.emplace(QValueLabel, Array::IndexOf(spl, QValueLabel));
		parsedHeader.emplace(QValueNotchLabel, Array::IndexOf(spl, QValueNotchLabel));

		parsedHeader.emplace(CrossTypeLabel, Array::IndexOf(spl, CrossTypeLabel));
		parsedHeader.emplace(LinkResiduesLabel, Array::IndexOf(spl, LinkResiduesLabel));
		parsedHeader.emplace(ProteinLinkSiteLabel, Array::IndexOf(spl, ProteinLinkSiteLabel));
		parsedHeader.emplace(RankLabel, Array::IndexOf(spl, RankLabel));
		parsedHeader.emplace(BetaPeptideProteinAccessionLabel, Array::IndexOf(spl, BetaPeptideProteinAccessionLabel));
		parsedHeader.emplace(BetaPeptideProteinLinkSiteLabel, Array::IndexOf(spl, BetaPeptideProteinLinkSiteLabel));
		parsedHeader.emplace(BetaPeptideBaseSequenceLabel, Array::IndexOf(spl, BetaPeptideBaseSequenceLabel));
		parsedHeader.emplace(BetaPeptideFullSequenceLabel, Array::IndexOf(spl, BetaPeptideFullSequenceLabel));
		parsedHeader.emplace(BetaPeptideTheoreticalMassLabel, Array::IndexOf(spl, BetaPeptideTheoreticalMassLabel));
		parsedHeader.emplace(BetaPeptideScoreLabel, Array::IndexOf(spl, BetaPeptideScoreLabel));
		parsedHeader.emplace(BetaPeptideRankLabel, Array::IndexOf(spl, BetaPeptideRankLabel));
		parsedHeader.emplace(BetaPeptideMatchedIonsLabel, Array::IndexOf(spl, BetaPeptideMatchedIonsLabel));
		parsedHeader.emplace(XLTotalScoreLabel, Array::IndexOf(spl, XLTotalScoreLabel));
		parsedHeader.emplace(ParentIonsLabel, Array::IndexOf(spl, ParentIonsLabel));

		return parsedHeader;
	}
}
