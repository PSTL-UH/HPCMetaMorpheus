#include "SearchParameters.h"

using namespace UsefulProteomicsDatabases;
using namespace EngineLayer;

namespace TaskLayer
{

	SearchParameters::SearchParameters()
	{
		// default search task parameters
		setDisposeOfFileWhenDone(true);
		setDoParsimony(true);
		setNoOneHitWonders(false);
		setModPeptidesAreDifferent(false);
		setDoQuantification(true);
		setQuantifyPpmTol(5);
		setSearchTarget(true);
		setDecoyType(getDecoyType()->Reverse);
		setDoHistogramAnalysis(false);
		setHistogramBinTolInDaltons(0.003);
		setDoLocalizationAnalysis(true);
		setWritePrunedDatabase(false);
		setKeepAllUniprotMods(true);
		setMassDiffAcceptorType(getMassDiffAcceptorType()::OneMM);
		setMaxFragmentSize(30000.0);
		setWriteMzId(true);
		setWritePepXml(false);

		setModsToWriteSelection(std::unordered_map<std::string, int>
		{
			{"N-linked glycosylation", 3},
			{"O-linked glycosylation", 3},
			{"Other glycosylation", 3},
			{"Common Biological", 3},
			{"Less Common", 3},
			{"Metal", 3},
			{"2+ nucleotide substitution", 3},
			{"1 nucleotide substitution", 3},
			{"UniProt", 2}
		});

		setWriteDecoys(true);
		setWriteContaminants(true);
		setLocalFdrCategories(std::vector<FdrCategory> {FdrCategory::FullySpecific});
	}

	bool SearchParameters::getDisposeOfFileWhenDone() const
	{
		return privateDisposeOfFileWhenDone;
	}

	void SearchParameters::setDisposeOfFileWhenDone(bool value)
	{
		privateDisposeOfFileWhenDone = value;
	}

	bool SearchParameters::getDoParsimony() const
	{
		return privateDoParsimony;
	}

	void SearchParameters::setDoParsimony(bool value)
	{
		privateDoParsimony = value;
	}

	bool SearchParameters::getModPeptidesAreDifferent() const
	{
		return privateModPeptidesAreDifferent;
	}

	void SearchParameters::setModPeptidesAreDifferent(bool value)
	{
		privateModPeptidesAreDifferent = value;
	}

	bool SearchParameters::getNoOneHitWonders() const
	{
		return privateNoOneHitWonders;
	}

	void SearchParameters::setNoOneHitWonders(bool value)
	{
		privateNoOneHitWonders = value;
	}

	bool SearchParameters::getMatchBetweenRuns() const
	{
		return privateMatchBetweenRuns;
	}

	void SearchParameters::setMatchBetweenRuns(bool value)
	{
		privateMatchBetweenRuns = value;
	}

	bool SearchParameters::getNormalize() const
	{
		return privateNormalize;
	}

	void SearchParameters::setNormalize(bool value)
	{
		privateNormalize = value;
	}

	double SearchParameters::getQuantifyPpmTol() const
	{
		return privateQuantifyPpmTol;
	}

	void SearchParameters::setQuantifyPpmTol(double value)
	{
		privateQuantifyPpmTol = value;
	}

	bool SearchParameters::getDoHistogramAnalysis() const
	{
		return privateDoHistogramAnalysis;
	}

	void SearchParameters::setDoHistogramAnalysis(bool value)
	{
		privateDoHistogramAnalysis = value;
	}

	bool SearchParameters::getSearchTarget() const
	{
		return privateSearchTarget;
	}

	void SearchParameters::setSearchTarget(bool value)
	{
		privateSearchTarget = value;
	}

	DecoyType *SearchParameters::getDecoyType() const
	{
		return privateDecoyType;
	}

	void SearchParameters::setDecoyType(DecoyType *value)
	{
		privateDecoyType = value;
	}

	TaskLayer::MassDiffAcceptorType SearchParameters::getMassDiffAcceptorType() const
	{
		return privateMassDiffAcceptorType;
	}

	void SearchParameters::setMassDiffAcceptorType(TaskLayer::MassDiffAcceptorType value)
	{
		privateMassDiffAcceptorType = value;
	}

	bool SearchParameters::getWritePrunedDatabase() const
	{
		return privateWritePrunedDatabase;
	}

	void SearchParameters::setWritePrunedDatabase(bool value)
	{
		privateWritePrunedDatabase = value;
	}

	bool SearchParameters::getKeepAllUniprotMods() const
	{
		return privateKeepAllUniprotMods;
	}

	void SearchParameters::setKeepAllUniprotMods(bool value)
	{
		privateKeepAllUniprotMods = value;
	}

	bool SearchParameters::getDoLocalizationAnalysis() const
	{
		return privateDoLocalizationAnalysis;
	}

	void SearchParameters::setDoLocalizationAnalysis(bool value)
	{
		privateDoLocalizationAnalysis = value;
	}

	bool SearchParameters::getDoQuantification() const
	{
		return privateDoQuantification;
	}

	void SearchParameters::setDoQuantification(bool value)
	{
		privateDoQuantification = value;
	}

	TaskLayer::SearchType SearchParameters::getSearchType() const
	{
		return privateSearchType;
	}

	void SearchParameters::setSearchType(TaskLayer::SearchType value)
	{
		privateSearchType = value;
	}

	std::vector<FdrCategory> SearchParameters::getLocalFdrCategories() const
	{
		return privateLocalFdrCategories;
	}

	void SearchParameters::setLocalFdrCategories(const std::vector<FdrCategory> &value)
	{
		privateLocalFdrCategories = value;
	}

	std::string SearchParameters::getCustomMdac() const
	{
		return privateCustomMdac;
	}

	void SearchParameters::setCustomMdac(const std::string &value)
	{
		privateCustomMdac = value;
	}

	double SearchParameters::getMaxFragmentSize() const
	{
		return privateMaxFragmentSize;
	}

	void SearchParameters::setMaxFragmentSize(double value)
	{
		privateMaxFragmentSize = value;
	}

	double SearchParameters::getHistogramBinTolInDaltons() const
	{
		return privateHistogramBinTolInDaltons;
	}

	void SearchParameters::setHistogramBinTolInDaltons(double value)
	{
		privateHistogramBinTolInDaltons = value;
	}

	std::unordered_map<std::string, int> SearchParameters::getModsToWriteSelection() const
	{
		return privateModsToWriteSelection;
	}

	void SearchParameters::setModsToWriteSelection(const std::unordered_map<std::string, int> &value)
	{
		privateModsToWriteSelection = value;
	}

	double SearchParameters::getMaximumMassThatFragmentIonScoreIsDoubled() const
	{
		return privateMaximumMassThatFragmentIonScoreIsDoubled;
	}

	void SearchParameters::setMaximumMassThatFragmentIonScoreIsDoubled(double value)
	{
		privateMaximumMassThatFragmentIonScoreIsDoubled = value;
	}

	bool SearchParameters::getWriteMzId() const
	{
		return privateWriteMzId;
	}

	void SearchParameters::setWriteMzId(bool value)
	{
		privateWriteMzId = value;
	}

	bool SearchParameters::getWritePepXml() const
	{
		return privateWritePepXml;
	}

	void SearchParameters::setWritePepXml(bool value)
	{
		privateWritePepXml = value;
	}

	bool SearchParameters::getWriteDecoys() const
	{
		return privateWriteDecoys;
	}

	void SearchParameters::setWriteDecoys(bool value)
	{
		privateWriteDecoys = value;
	}

	bool SearchParameters::getWriteContaminants() const
	{
		return privateWriteContaminants;
	}

	void SearchParameters::setWriteContaminants(bool value)
	{
		privateWriteContaminants = value;
	}
}
