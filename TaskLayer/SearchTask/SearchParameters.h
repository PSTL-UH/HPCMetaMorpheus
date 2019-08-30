#pragma once

#include "MassDiffAcceptorType.h"
#include "SearchType.h"
#include "../../EngineLayer/ProteinScoringAndFdr/FdrCategory.h"
#include <string>
#include <unordered_map>
#include <vector>

using namespace UsefulProteomicsDatabases;
using namespace EngineLayer;

namespace TaskLayer
{
	class SearchParameters
	{
	private:
		bool privateDisposeOfFileWhenDone = false;
		bool privateDoParsimony = false;
		bool privateModPeptidesAreDifferent = false;
		bool privateNoOneHitWonders = false;
		bool privateMatchBetweenRuns = false;
		bool privateNormalize = false;
		double privateQuantifyPpmTol = 0;
		bool privateDoHistogramAnalysis = false;
		bool privateSearchTarget = false;
		DecoyType *privateDecoyType;
		TaskLayer::MassDiffAcceptorType privateMassDiffAcceptorType = static_cast<TaskLayer::MassDiffAcceptorType>(0);
		bool privateWritePrunedDatabase = false;
		bool privateKeepAllUniprotMods = false;
		bool privateDoLocalizationAnalysis = false;
		bool privateDoQuantification = false;
		TaskLayer::SearchType privateSearchType = static_cast<TaskLayer::SearchType>(0);
		std::vector<FdrCategory> privateLocalFdrCategories;
		std::string privateCustomMdac;
		double privateMaxFragmentSize = 0;
		double privateHistogramBinTolInDaltons = 0;
		std::unordered_map<std::string, int> privateModsToWriteSelection;
		double privateMaximumMassThatFragmentIonScoreIsDoubled = 0;
		bool privateWriteMzId = false;
		bool privateWritePepXml = false;
		bool privateWriteDecoys = false;
		bool privateWriteContaminants = false;

	public:
		SearchParameters();

		bool getDisposeOfFileWhenDone() const;
		void setDisposeOfFileWhenDone(bool value);
		bool getDoParsimony() const;
		void setDoParsimony(bool value);
		bool getModPeptidesAreDifferent() const;
		void setModPeptidesAreDifferent(bool value);
		bool getNoOneHitWonders() const;
		void setNoOneHitWonders(bool value);
		bool getMatchBetweenRuns() const;
		void setMatchBetweenRuns(bool value);
		bool getNormalize() const;
		void setNormalize(bool value);
		double getQuantifyPpmTol() const;
		void setQuantifyPpmTol(double value);
		bool getDoHistogramAnalysis() const;
		void setDoHistogramAnalysis(bool value);
		bool getSearchTarget() const;
		void setSearchTarget(bool value);
		DecoyType *getDecoyType() const;
		void setDecoyType(DecoyType *value);
		TaskLayer::MassDiffAcceptorType getMassDiffAcceptorType() const;
		void setMassDiffAcceptorType(TaskLayer::MassDiffAcceptorType value);
		bool getWritePrunedDatabase() const;
		void setWritePrunedDatabase(bool value);
		bool getKeepAllUniprotMods() const;
		void setKeepAllUniprotMods(bool value);
		bool getDoLocalizationAnalysis() const;
		void setDoLocalizationAnalysis(bool value);
		bool getDoQuantification() const;
		void setDoQuantification(bool value);
		TaskLayer::SearchType getSearchType() const;
		void setSearchType(TaskLayer::SearchType value);
		std::vector<FdrCategory> getLocalFdrCategories() const;
		void setLocalFdrCategories(const std::vector<FdrCategory> &value);
		std::string getCustomMdac() const;
		void setCustomMdac(const std::string &value);
		double getMaxFragmentSize() const;
		void setMaxFragmentSize(double value);
		double getHistogramBinTolInDaltons() const;
		void setHistogramBinTolInDaltons(double value);
		std::unordered_map<std::string, int> getModsToWriteSelection() const;
		void setModsToWriteSelection(const std::unordered_map<std::string, int> &value);
		double getMaximumMassThatFragmentIonScoreIsDoubled() const;
		void setMaximumMassThatFragmentIonScoreIsDoubled(double value);
		bool getWriteMzId() const;
		void setWriteMzId(bool value);
		bool getWritePepXml() const;
		void setWritePepXml(bool value);
		bool getWriteDecoys() const;
		void setWriteDecoys(bool value);
		bool getWriteContaminants() const;
		void setWriteContaminants(bool value);
	};
}
