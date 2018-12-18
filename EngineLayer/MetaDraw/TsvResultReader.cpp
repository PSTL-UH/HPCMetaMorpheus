#include "TsvResultReader.h"
#include "MetaDrawPsm.h"
#include "../MetaMorpheusException.h"


namespace EngineLayer
{

const std::wstring TsvResultReader::FullSequenceLabel = L"Full Sequence";
const std::wstring TsvResultReader::Ms2ScanNumberLabel = L"Scan Number";
const std::wstring TsvResultReader::FilenameLabel = L"File Name";
const std::wstring TsvResultReader::TotalIonCurrentLabel = L"Total Ion Current";
const std::wstring TsvResultReader::PrecursorScanNumLabel = L"Precursor Scan Number";
const std::wstring TsvResultReader::PrecursorChargeLabel = L"Precursor Charge";
const std::wstring TsvResultReader::PrecursorMzLabel = L"Precursor MZ";
const std::wstring TsvResultReader::PrecursorMassLabel = L"Precursor Mass";
const std::wstring TsvResultReader::ScoreLabel = L"Score";
const std::wstring TsvResultReader::DeltaScoreLabel = L"Delta Score";
const std::wstring TsvResultReader::NotchLabel = L"Notch";
const std::wstring TsvResultReader::BaseSeqLabel = L"Base Sequence";
const std::wstring TsvResultReader::EssentialSeqLabel = L"Essential Sequence";
const std::wstring TsvResultReader::MissedCleavageLabel = L"Missed Cleavages";
const std::wstring TsvResultReader::PeptideMonoMassLabel = L"Peptide Monoisotopic Mass";
const std::wstring TsvResultReader::MassDiffDaLabel = L"Mass Diff (Da)";
const std::wstring TsvResultReader::MassDiffPpmLabel = L"Mass Diff (ppm)";
const std::wstring TsvResultReader::ProteinAccessionLabel = L"Protein Accession";
const std::wstring TsvResultReader::ProteinNameLabel = L"Protein Name";
const std::wstring TsvResultReader::GeneNameLabel = L"Gene Name";
const std::wstring TsvResultReader::OrganismNameLabel = L"Organism Name";
const std::wstring TsvResultReader::PeptideDesicriptionLabel = L"Peptide Description";
const std::wstring TsvResultReader::StartAndEndResiduesInProteinLabel = L"Start and End Residues In Protein";
const std::wstring TsvResultReader::PreviousAminoAcidLabel = L"Previous Amino Acid";
const std::wstring TsvResultReader::NextAminoAcidLabel = L"Next Amino Acid";
const std::wstring TsvResultReader::DecoyContamTargetLabel = L"Decoy/Contaminant/Target";
const std::wstring TsvResultReader::MatchedIonsLabel = L"Matched Ion Mass-To-Charge Ratios";
const std::wstring TsvResultReader::QValueLabel = L"QValue";
const std::wstring TsvResultReader::QValueNotchLabel = L"QValue Notch";
const std::wstring TsvResultReader::CrossTypeLabel = L"Cross Type";
const std::wstring TsvResultReader::LinkResiduesLabel = L"Link Residues";
const std::wstring TsvResultReader::ProteinLinkSiteLabel = L"Protein Link Site";
const std::wstring TsvResultReader::RankLabel = L"Rank";
const std::wstring TsvResultReader::BetaPeptideProteinAccessionLabel = L"Beta Peptide Protein Accession";
const std::wstring TsvResultReader::BetaPeptideProteinLinkSiteLabel = L"Beta Peptide Protein LinkSite";
const std::wstring TsvResultReader::BetaPeptideBaseSequenceLabel = L"Beta Peptide Base Sequence";
const std::wstring TsvResultReader::BetaPeptideFullSequenceLabel = L"Beta Peptide Full Sequence";
const std::wstring TsvResultReader::BetaPeptideTheoreticalMassLabel = L"Beta Peptide Theoretical Mass";
const std::wstring TsvResultReader::BetaPeptideScoreLabel = L"Beta Peptide Score";
const std::wstring TsvResultReader::BetaPeptideRankLabel = L"Beta Peptide Rank";
const std::wstring TsvResultReader::BetaPeptideMatchedIonsLabel = L"Beta Peptide Matched Ion Mass-To-Charge Ratios";
const std::wstring TsvResultReader::XLTotalScoreLabel = L"XL Total Score";
const std::wstring TsvResultReader::ParentIonsLabel = L"Parent Ions";
std::vector<wchar_t> const TsvResultReader::Split = {L'\t'};

	std::vector<MetaDrawPsm*> TsvResultReader::ReadTsv(const std::wstring &filePath, std::vector<std::wstring> &warnings)
	{
		std::vector<MetaDrawPsm*> psms;
		warnings = std::vector<std::wstring>();

		StreamReader *reader = nullptr;
		try
		{
			reader = new StreamReader(filePath);
		}
		catch (const std::runtime_error &e)
		{
		   delete reader;
		   throw MetaMorpheusException(L"Could not read file: " + e.what());
		}

		int lineCount = 0;

		std::wstring line;
		std::unordered_map<std::wstring, int> parsedHeader;

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
				warnings.push_back(L"Could not read line: " + std::to_wstring(lineCount));
			}
		}

		reader->Close();

		if ((lineCount - 1) != psms.size())
		{
			warnings.push_back(L"Warning: " + std::to_wstring((lineCount - 1) - psms.size()) + L" PSMs were not read.");
		}

		delete reader;
		return psms;
	}

	std::unordered_map<std::wstring, int> TsvResultReader::ParseHeader(const std::wstring &header)
	{
		auto parsedHeader = std::unordered_map<std::wstring, int>();
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
