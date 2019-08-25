#include "MetaDrawPsm.h"
#include "TsvResultReader.h"

using namespace Proteomics::Fragmentation;
using namespace Chemistry;

namespace EngineLayer
{

Regex *const MetaDrawPsm::IonParser = new Regex(LR"(([a-zA-Z]+)(\d+))");
std::vector<char> const MetaDrawPsm::MzSplit = {'[', ',', ']', ';'};

	std::string MetaDrawPsm::getFullSequence() const
	{
		return privateFullSequence;
	}

	int MetaDrawPsm::getMs2ScanNumber() const
	{
		return privateMs2ScanNumber;
	}

	std::string MetaDrawPsm::getFilename() const
	{
		return privateFilename;
	}

	int MetaDrawPsm::getPrecursorScanNum() const
	{
		return privatePrecursorScanNum;
	}

	int MetaDrawPsm::getPrecursorCharge() const
	{
		return privatePrecursorCharge;
	}

	double MetaDrawPsm::getPrecursorMz() const
	{
		return privatePrecursorMz;
	}

	double MetaDrawPsm::getPrecursorMass() const
	{
		return privatePrecursorMass;
	}

	double MetaDrawPsm::getScore() const
	{
		return privateScore;
	}

	std::string MetaDrawPsm::getProteinAccession() const
	{
		return privateProteinAccession;
	}

	std::vector<MatchedFragmentIon*> MetaDrawPsm::getMatchedIons() const
	{
		return privateMatchedIons;
	}

	double MetaDrawPsm::getQValue() const
	{
		return privateQValue;
	}

	std::optional<double> MetaDrawPsm::getTotalIonCurrent() const
	{
		return privateTotalIonCurrent;
	}

	std::optional<double> MetaDrawPsm::getDeltaScore() const
	{
		return privateDeltaScore;
	}

	std::string MetaDrawPsm::getNotch() const
	{
		return privateNotch;
	}

	std::string MetaDrawPsm::getBaseSeq() const
	{
		return privateBaseSeq;
	}

	std::string MetaDrawPsm::getEssentialSeq() const
	{
		return privateEssentialSeq;
	}

	std::string MetaDrawPsm::getMissedCleavage() const
	{
		return privateMissedCleavage;
	}

	std::string MetaDrawPsm::getPeptideMonoMass() const
	{
		return privatePeptideMonoMass;
	}

	std::string MetaDrawPsm::getMassDiffDa() const
	{
		return privateMassDiffDa;
	}

	std::string MetaDrawPsm::getMassDiffPpm() const
	{
		return privateMassDiffPpm;
	}

	std::string MetaDrawPsm::getProteinName() const
	{
		return privateProteinName;
	}

	std::string MetaDrawPsm::getGeneName() const
	{
		return privateGeneName;
	}

	std::string MetaDrawPsm::getOrganismName() const
	{
		return privateOrganismName;
	}

	std::string MetaDrawPsm::getPeptideDesicription() const
	{
		return privatePeptideDesicription;
	}

	std::string MetaDrawPsm::getStartAndEndResiduesInProtein() const
	{
		return privateStartAndEndResiduesInProtein;
	}

	std::string MetaDrawPsm::getPreviousAminoAcid() const
	{
		return privatePreviousAminoAcid;
	}

	std::string MetaDrawPsm::getNextAminoAcid() const
	{
		return privateNextAminoAcid;
	}

	std::string MetaDrawPsm::getDecoyContamTarget() const
	{
		return privateDecoyContamTarget;
	}

	std::optional<double> MetaDrawPsm::getQValueNotch() const
	{
		return privateQValueNotch;
	}

	std::string MetaDrawPsm::getCrossType() const
	{
		return privateCrossType;
	}

	std::string MetaDrawPsm::getLinkResidues() const
	{
		return privateLinkResidues;
	}

	std::optional<int> MetaDrawPsm::getProteinLinkSite() const
	{
		return privateProteinLinkSite;
	}

	std::optional<int> MetaDrawPsm::getRank() const
	{
		return privateRank;
	}

	std::string MetaDrawPsm::getBetaPeptideProteinAccession() const
	{
		return privateBetaPeptideProteinAccession;
	}

	std::optional<int> MetaDrawPsm::getBetaPeptideProteinLinkSite() const
	{
		return privateBetaPeptideProteinLinkSite;
	}

	std::string MetaDrawPsm::getBetaPeptideBaseSequence() const
	{
		return privateBetaPeptideBaseSequence;
	}

	std::string MetaDrawPsm::getBetaPeptideFullSequence() const
	{
		return privateBetaPeptideFullSequence;
	}

	std::string MetaDrawPsm::getBetaPeptideTheoreticalMass() const
	{
		return privateBetaPeptideTheoreticalMass;
	}

	std::optional<double> MetaDrawPsm::getBetaPeptideScore() const
	{
		return privateBetaPeptideScore;
	}

	std::optional<int> MetaDrawPsm::getBetaPeptideRank() const
	{
		return privateBetaPeptideRank;
	}

	std::vector<MatchedFragmentIon*> MetaDrawPsm::getBetaPeptideMatchedIons() const
	{
		return privateBetaPeptideMatchedIons;
	}

	std::optional<double> MetaDrawPsm::getXLTotalScore() const
	{
		return privateXLTotalScore;
	}

	std::string MetaDrawPsm::getParentIons() const
	{
		return privateParentIons;
	}

	MetaDrawPsm::MetaDrawPsm(const std::string &line, std::vector<char> &split, std::unordered_map<std::string, int> &parsedHeader)
	{
		auto spl = line.Split(split);

		//Required properties
		Filename = StringHelper::trim(spl[parsedHeader[TsvResultReader::FilenameLabel]]);
		Ms2ScanNumber = std::stoi(spl[parsedHeader[TsvResultReader::Ms2ScanNumberLabel]]);
		PrecursorScanNum = std::stoi(StringHelper::trim(spl[parsedHeader[TsvResultReader::PrecursorScanNumLabel]]));
		PrecursorCharge = static_cast<int>(std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::PrecursorChargeLabel]])));
		PrecursorMz = std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::PrecursorMzLabel]]));
		PrecursorMass = std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::PrecursorMassLabel]]));
		BaseSeq = StringHelper::trim(spl[parsedHeader[TsvResultReader::BaseSeqLabel]]);
		FullSequence = spl[parsedHeader[TsvResultReader::FullSequenceLabel]];
		PeptideMonoMass = StringHelper::trim(spl[parsedHeader[TsvResultReader::PeptideMonoMassLabel]]);
		Score = std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::ScoreLabel]]));
		DecoyContamTarget = StringHelper::trim(spl[parsedHeader[TsvResultReader::DecoyContamTargetLabel]]);
		QValue = std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::QValueLabel]]));
		MatchedIons = ReadFragmentIonsFromString(StringHelper::trim(spl[parsedHeader[TsvResultReader::MatchedIonsLabel]]), getBaseSeq());

		//For general psms
		TotalIonCurrent = (parsedHeader[TsvResultReader::TotalIonCurrentLabel] < 0) ? std::nullopt : static_cast<std::optional<double>>(std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::TotalIonCurrentLabel]])));
		DeltaScore = (parsedHeader[TsvResultReader::DeltaScoreLabel] < 0) ? std::nullopt : static_cast<std::optional<double>>(std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::DeltaScoreLabel]])));
		Notch = (parsedHeader[TsvResultReader::NotchLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::NotchLabel]]);
		EssentialSeq = (parsedHeader[TsvResultReader::EssentialSeqLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::EssentialSeqLabel]]);
		MissedCleavage = (parsedHeader[TsvResultReader::MissedCleavageLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::MissedCleavageLabel]]);
		MassDiffDa = (parsedHeader[TsvResultReader::MassDiffDaLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::MassDiffDaLabel]]);
		MassDiffPpm = (parsedHeader[TsvResultReader::MassDiffPpmLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::MassDiffPpmLabel]]);
		ProteinAccession = (parsedHeader[TsvResultReader::ProteinAccessionLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::ProteinAccessionLabel]]);
		ProteinName = (parsedHeader[TsvResultReader::ProteinNameLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::ProteinNameLabel]]);
		GeneName = (parsedHeader[TsvResultReader::GeneNameLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::GeneNameLabel]]);
		OrganismName = (parsedHeader[TsvResultReader::OrganismNameLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::OrganismNameLabel]]);
		PeptideDesicription = (parsedHeader[TsvResultReader::PeptideDesicriptionLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::PeptideDesicriptionLabel]]);
		StartAndEndResiduesInProtein = (parsedHeader[TsvResultReader::StartAndEndResiduesInProteinLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::StartAndEndResiduesInProteinLabel]]);
		PreviousAminoAcid = (parsedHeader[TsvResultReader::PreviousAminoAcidLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::PreviousAminoAcidLabel]]);
		NextAminoAcid = (parsedHeader[TsvResultReader::NextAminoAcidLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::NextAminoAcidLabel]]);
		QValueNotch = (parsedHeader[TsvResultReader::QValueNotchLabel] < 0) ? std::nullopt : static_cast<std::optional<double>>(std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::QValueNotchLabel]])));

		//For crosslinks
		CrossType = (parsedHeader[TsvResultReader::CrossTypeLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::CrossTypeLabel]]);
		LinkResidues = (parsedHeader[TsvResultReader::LinkResiduesLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::LinkResiduesLabel]]);
		ProteinLinkSite = (parsedHeader[TsvResultReader::ProteinLinkSiteLabel] < 0) ? std::nullopt : static_cast<std::optional<int>>(std::stoi(StringHelper::trim(spl[parsedHeader[TsvResultReader::ProteinLinkSiteLabel]])));
		Rank = (parsedHeader[TsvResultReader::RankLabel] < 0) ? std::nullopt : static_cast<std::optional<int>>(std::stoi(StringHelper::trim(spl[parsedHeader[TsvResultReader::RankLabel]])));
		BetaPeptideProteinAccession = (parsedHeader[TsvResultReader::BetaPeptideProteinAccessionLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::BetaPeptideProteinAccessionLabel]]);
		BetaPeptideProteinLinkSite = (parsedHeader[TsvResultReader::BetaPeptideProteinLinkSiteLabel] < 0) ? std::nullopt : static_cast<std::optional<int>>(std::stoi(StringHelper::trim(spl[parsedHeader[TsvResultReader::BetaPeptideProteinLinkSiteLabel]])));
		BetaPeptideBaseSequence = (parsedHeader[TsvResultReader::BetaPeptideBaseSequenceLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::BetaPeptideBaseSequenceLabel]]);
		BetaPeptideFullSequence = (parsedHeader[TsvResultReader::BetaPeptideFullSequenceLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::BetaPeptideFullSequenceLabel]]);
		BetaPeptideTheoreticalMass = (parsedHeader[TsvResultReader::BetaPeptideTheoreticalMassLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::BetaPeptideTheoreticalMassLabel]]);
		BetaPeptideScore = (parsedHeader[TsvResultReader::BetaPeptideScoreLabel] < 0) ? std::nullopt : static_cast<std::optional<double>>(std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::BetaPeptideScoreLabel]])));
		BetaPeptideRank = (parsedHeader[TsvResultReader::BetaPeptideRankLabel] < 0) ? std::nullopt : static_cast<std::optional<int>>(std::stoi(StringHelper::trim(spl[parsedHeader[TsvResultReader::BetaPeptideRankLabel]])));
		BetaPeptideMatchedIons = (parsedHeader[TsvResultReader::BetaPeptideMatchedIonsLabel] < 0) ? nullptr : ReadFragmentIonsFromString(StringHelper::trim(spl[parsedHeader[TsvResultReader::BetaPeptideMatchedIonsLabel]]), getBetaPeptideBaseSequence());
		XLTotalScore = (parsedHeader[TsvResultReader::XLTotalScoreLabel] < 0) ? std::nullopt : static_cast<std::optional<double>>(std::stod(StringHelper::trim(spl[parsedHeader[TsvResultReader::XLTotalScoreLabel]])));
		ParentIons = (parsedHeader[TsvResultReader::ParentIonsLabel] < 0) ? "" : StringHelper::trim(spl[parsedHeader[TsvResultReader::ParentIonsLabel]]);
	}

	std::vector<MatchedFragmentIon*> MetaDrawPsm::ReadFragmentIonsFromString(const std::string &matchedMzString, const std::string &peptideBaseSequence)
	{
		auto peaks = matchedMzString.Split(MzSplit, StringSplitOptions::RemoveEmptyEntries).Select([&] (std::any v)
		{
			v->Trim();
		}).ToList();
		peaks.RemoveAll([&] (std::any p)
		{
			p->Contains("\"");
		});

		std::vector<MatchedFragmentIon*> matchedIons;

		for (auto peak : peaks)
		{
			auto split = peak.Split(std::vector<char> {'+', ':'});

			std::string ionTypeAndNumber = split[0];
			Match *result = IonParser->Match(ionTypeAndNumber);

			ProductType *productType = std::any_cast<ProductType*>(Enum::Parse(ProductType::typeid, result->Groups[1]->Value));

			int fragmentNumber = std::stoi(result->Groups[2]->Value);
			int z = std::stoi(split[1]);
			double mz = std::stod(split[2]);
			double neutralLoss = 0;

			// check for neutral loss
			if (ionTypeAndNumber.find("-") != std::string::npos)
			{
				std::string temp = StringHelper::replace(ionTypeAndNumber, "(", "");
				temp = StringHelper::replace(temp, ")", "");
				auto split2 = StringHelper::split(temp, '-');
				neutralLoss = std::stod(split2[1]);
			}

			FragmentationTerminus *terminus = FragmentationTerminus::None;
			if (TerminusSpecificProductTypes::ProductTypeToFragmentationTerminus->ContainsKey(productType))
			{
				terminus = TerminusSpecificProductTypes::ProductTypeToFragmentationTerminus[productType];
			}

			int aminoAcidPosition = fragmentNumber;
			if (terminus == FragmentationTerminus::C)
			{
				aminoAcidPosition = peptideBaseSequence.length() - fragmentNumber;
			}

			auto t = new NeutralTerminusFragment(terminus, mz.ToMass(z) - DissociationTypeCollection::GetMassShiftFromProductType(productType), fragmentNumber, aminoAcidPosition);
			Product *p = new Product(productType, t, neutralLoss);
			MatchedFragmentIon tempVar(p, mz, 1.0, z);
			matchedIons.push_back(&tempVar);

//C# TO C++ CONVERTER TODO TASK: A 'delete p' statement was not added since p was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete t' statement was not added since t was passed to a method or constructor. Handle memory management manually.
		}

		return matchedIons;
	}
}
