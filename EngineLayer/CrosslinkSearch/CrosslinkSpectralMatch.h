#pragma once

#include "../PeptideSpectralMatch.h"
#include "PsmCrossType.h"
#include <string>
#include <unordered_map>
#include <vector>
#include "stringhelper.h"
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class Ms2ScanWithSpecificMass; }
namespace EngineLayer { class PeptideSpectralMatch; }

using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace CrosslinkSearch
	{
		class CrosslinkSpectralMatch : public PeptideSpectralMatch
		{
		private:
			CrosslinkSpectralMatch *privateBetaPeptide;
			std::vector<int> privateLinkPositions;
			double privateDeltaScore = 0;
			double privateXLTotalScore = 0;
			int privateXlProteinPos = 0;
			std::vector<int> privateXlRank;
			std::wstring privateParentIonExist;
			int privateParentIonExistNum = 0;
			std::vector<int> privateParentIonMaxIntensityRanks;
			PsmCrossType privateCrossType = static_cast<PsmCrossType>(0);

		public:
			CrosslinkSpectralMatch(PeptideWithSetModifications *theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass *scan, DigestionParams *digestionParams, std::vector<MatchedFragmentIon*> &matchedFragmentIons);

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
				std::wstring getParentIonExist() const;
				void setParentIonExist(const std::wstring &value);
				int getParentIonExistNum() const;
				void setParentIonExistNum(int value);
				std::vector<int> getParentIonMaxIntensityRanks() const;
				void setParentIonMaxIntensityRanks(const std::vector<int> &value);
				PsmCrossType getCrossType() const;
				void setCrossType(PsmCrossType value);

			static std::vector<int> GetPossibleCrosslinkerModSites(std::vector<wchar_t> &crosslinkerModSites, PeptideWithSetModifications *peptide);

			/// <summary>
			/// Rank experimental mass spectral peaks by intensity
			/// </summary>
			static std::vector<int> GenerateIntensityRanks(std::vector<double> &experimental_intensities);

			static std::wstring GetTabSepHeaderCross();

			static std::wstring GetTabSepHeaderSingle();

			static std::wstring GetTabSepHeaderGlyco();

			std::wstring ToString() override;

			static std::unordered_map<std::wstring, std::wstring> MatchedIonDataDictionary(PeptideSpectralMatch *psm);
		};
	}
}
