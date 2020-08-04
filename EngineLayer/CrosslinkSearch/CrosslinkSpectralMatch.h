#pragma once

#include "../PeptideSpectralMatch.h"
#include "PsmCrossType.h"
#include <string>
#include <unordered_map>
#include <vector>
#include "stringhelper.h"
#include "stringbuilder.h"

#include "../Ms2ScanWithSpecificMass.h"
#include "../PeptideSpectralMatch.h"

#include "Proteomics/Proteomics.h"
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    namespace CrosslinkSearch
    {
        class CrosslinkSpectralMatch : public PeptideSpectralMatch
        {
        private:
            CrosslinkSpectralMatch *privateBetaPeptide=nullptr;
            std::vector<int> privateLinkPositions;
            double privateDeltaScore = 0;
            double privateXLTotalScore = 0;
            int privateXlProteinPos = 0;
            std::vector<int> privateXlRank;
            std::string privateParentIonExist;
            int privateParentIonExistNum = 0;
            std::vector<int> privateParentIonMaxIntensityRanks;
            PsmCrossType privateCrossType = static_cast<PsmCrossType>(0);
            
        public:
            CrosslinkSpectralMatch(PeptideWithSetModifications *theBestPeptide, int notch, double score, int scanIndex,
                                   Ms2ScanWithSpecificMass *scan, DigestionParams *digestionParams,
                                   std::vector<MatchedFragmentIon*> &matchedFragmentIons);
            
            CrosslinkSpectralMatch *getBetaPeptide() const;
            void setBetaPeptide(CrosslinkSpectralMatch *value);
            std::vector<int> getLinkPositions() const;
            void setLinkPositions(const std::vector<int> &value);
            double getDeltaScore() const;
            void setDeltaScore(double value);
            double getXLTotalScore() const;
            void setXLTotalScore(double value);
            int getXlProteinPos() const;
            void setXlProteinPos(int value);
            std::vector<int> getXlRank() const;
            void setXlRank(const std::vector<int> &value);
            std::string getParentIonExist() const;
            void setParentIonExist(const std::string &value);
            int getParentIonExistNum() const;
            void setParentIonExistNum(int value);
            std::vector<int> getParentIonMaxIntensityRanks() const;
            void setParentIonMaxIntensityRanks(const std::vector<int> &value);
            PsmCrossType getCrossType() const;
            void setCrossType(PsmCrossType value);
            
            static std::vector<int> GetPossibleCrosslinkerModSites(std::vector<char> &crosslinkerModSites,
                                                                   PeptideWithSetModifications *peptide);
            
            /// <summary>
			/// Rank experimental mass spectral peaks by intensity
			/// </summary>
			static std::vector<int> GenerateIntensityRanks(std::vector<double> &experimental_intensities);

			static std::string GetTabSepHeaderCross();

			static std::string GetTabSepHeaderSingle();

			static std::string GetTabSepHeaderGlyco();

			std::string ToString();

			static std::vector<std::tuple<std::string, std::string>> MatchedIonDataDictionary(PeptideSpectralMatch *psm);
		};
	}
}
