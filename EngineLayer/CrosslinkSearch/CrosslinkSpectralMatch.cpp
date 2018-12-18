#include "CrosslinkSpectralMatch.h"
#include "../Ms2ScanWithSpecificMass.h"

using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
	namespace CrosslinkSearch
	{

		CrosslinkSpectralMatch::CrosslinkSpectralMatch(PeptideWithSetModifications *theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass *scan, DigestionParams *digestionParams, std::vector<MatchedFragmentIon*> &matchedFragmentIons) : PeptideSpectralMatch(theBestPeptide, notch, score, scanIndex, scan, digestionParams, matchedFragmentIons)
		{
			this->setXLTotalScore(score);
		}

		CrosslinkSpectralMatch *CrosslinkSpectralMatch::getBetaPeptide() const
		{
			return privateBetaPeptide;
		}

		void CrosslinkSpectralMatch::setBetaPeptide(CrosslinkSpectralMatch *value)
		{
			privateBetaPeptide = value;
		}

		std::vector<int> CrosslinkSpectralMatch::getLinkPositions() const
		{
			return privateLinkPositions;
		}

		void CrosslinkSpectralMatch::setLinkPositions(const std::vector<int> &value)
		{
			privateLinkPositions = value;
		}

		double CrosslinkSpectralMatch::getDeltaScore() const
		{
			return privateDeltaScore;
		}

		void CrosslinkSpectralMatch::setDeltaScore(double value)
		{
			privateDeltaScore = value;
		}

		double CrosslinkSpectralMatch::getXLTotalScore() const
		{
			return privateXLTotalScore;
		}

		void CrosslinkSpectralMatch::setXLTotalScore(double value)
		{
			privateXLTotalScore = value;
		}

		int CrosslinkSpectralMatch::getXlProteinPos() const
		{
			return privateXlProteinPos;
		}

		void CrosslinkSpectralMatch::setXlProteinPos(int value)
		{
			privateXlProteinPos = value;
		}

		std::vector<int> CrosslinkSpectralMatch::getXlRank() const
		{
			return privateXlRank;
		}

		void CrosslinkSpectralMatch::setXlRank(const std::vector<int> &value)
		{
			privateXlRank = value;
		}

		std::wstring CrosslinkSpectralMatch::getParentIonExist() const
		{
			return privateParentIonExist;
		}

		void CrosslinkSpectralMatch::setParentIonExist(const std::wstring &value)
		{
			privateParentIonExist = value;
		}

		int CrosslinkSpectralMatch::getParentIonExistNum() const
		{
			return privateParentIonExistNum;
		}

		void CrosslinkSpectralMatch::setParentIonExistNum(int value)
		{
			privateParentIonExistNum = value;
		}

		std::vector<int> CrosslinkSpectralMatch::getParentIonMaxIntensityRanks() const
		{
			return privateParentIonMaxIntensityRanks;
		}

		void CrosslinkSpectralMatch::setParentIonMaxIntensityRanks(const std::vector<int> &value)
		{
			privateParentIonMaxIntensityRanks = value;
		}

		PsmCrossType CrosslinkSpectralMatch::getCrossType() const
		{
			return privateCrossType;
		}

		void CrosslinkSpectralMatch::setCrossType(PsmCrossType value)
		{
			privateCrossType = value;
		}

		std::vector<int> CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(std::vector<wchar_t> &crosslinkerModSites, PeptideWithSetModifications *peptide)
		{
			std::vector<int> possibleXlPositions;
			bool wildcard = crosslinkerModSites.Any([&] (std::any p)
			{
				return p == L'X';
			});

			for (int r = 0; r < peptide->BaseSequence->Length; r++)
			{
				if (crosslinkerModSites.Contains(peptide->BaseSequence[r]) || wildcard)
				{
					possibleXlPositions.push_back(r + 1);
				}
			}

			return possibleXlPositions;
		}

		std::vector<int> CrosslinkSpectralMatch::GenerateIntensityRanks(std::vector<double> &experimental_intensities)
		{
			auto y = experimental_intensities.ToArray();
			auto x = Enumerable::Range(1, y.size()).OrderBy([&] (std::any p)
			{
				return p;
			})->ToArray();
			Array::Sort(y, x);
			auto experimental_intensities_rank = Enumerable::Range(1, y.size()).OrderByDescending([&] (std::any p)
			{
				return p;
			})->ToArray();
			Array::Sort(x, experimental_intensities_rank);
			return experimental_intensities_rank;
		}

		std::wstring CrosslinkSpectralMatch::GetTabSepHeaderCross()
		{
			auto sb = new StringBuilder();
			sb->append(L"File Name" + StringHelper::toString(L'\t'));
			sb->append(L"Scan Number" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor Scan Number" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor MZ" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor Charge" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor Mass" + StringHelper::toString(L'\t'));
			sb->append(L"Cross Type" + StringHelper::toString(L'\t'));
			sb->append(std::wstring(L"Link Residues") + L"\t");

			sb->append(L"Peptide" + StringHelper::toString(L'\t'));
			sb->append(L"Protein Accession" + StringHelper::toString(L'\t'));
			sb->append(L"Protein Link Site" + StringHelper::toString(L'\t'));
			sb->append(L"Base Sequence" + StringHelper::toString(L'\t'));
			sb->append(L"Full Sequence" + StringHelper::toString(L'\t'));
			sb->append(L"Peptide Monoisotopic Mass" + StringHelper::toString(L'\t'));
			sb->append(L"Score" + StringHelper::toString(L'\t'));
			sb->append(L"Rank" + StringHelper::toString(L'\t'));

			sb->append(L"Matched Ion Series" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Mass-To-Charge Ratios" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Mass Diff (Ppm)" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Intensities" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Counts" + StringHelper::toString(L'\t'));

			sb->append(L"Beta Peptide" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Protein Accession" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Protein LinkSite" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Base Sequence" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Full Sequence" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Theoretical Mass" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Score" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Rank" + StringHelper::toString(L'\t'));

			sb->append(L"Beta Peptide Matched Ion Series" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Matched Ion Mass-To-Charge Ratios" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Matched Ion Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Matched Ion Mass Diff (Ppm)" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Matched Ion Intensities" + StringHelper::toString(L'\t'));
			sb->append(L"Beta Peptide Matched Ion Counts" + StringHelper::toString(L'\t'));

			sb->append(L"Summary" + StringHelper::toString(L'\t'));
			sb->append(L"XL Total Score" + StringHelper::toString(L'\t'));
			sb->append(L"Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append(L"Parent Ions" + StringHelper::toString(L'\t'));
			sb->append(L"ParentIonsNum" + StringHelper::toString(L'\t'));
			sb->append(L"ParentIonMaxIntensityRank" + StringHelper::toString(L'\t'));
			sb->append(L"Decoy/Contaminant/Target" + StringHelper::toString(L'\t'));
			sb->append(L"QValue" + StringHelper::toString(L'\t'));


			delete sb;
			return sb->toString();
		}

		std::wstring CrosslinkSpectralMatch::GetTabSepHeaderSingle()
		{
			auto sb = new StringBuilder();
			sb->append(L"File Name" + StringHelper::toString(L'\t'));
			sb->append(L"Scan Number" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor Scan Number" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor MZ" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor Charge" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor Mass" + StringHelper::toString(L'\t'));
			sb->append(L"Cross Type" + StringHelper::toString(L'\t'));
			sb->append(std::wstring(L"Link Residues") + L"\t");

			sb->append(L"Peptide" + StringHelper::toString(L'\t'));
			sb->append(L"Protein Accession" + StringHelper::toString(L'\t'));
			sb->append(L"Protein Link Site" + StringHelper::toString(L'\t'));
			sb->append(L"Base Sequence" + StringHelper::toString(L'\t'));
			sb->append(L"Full Sequence" + StringHelper::toString(L'\t'));
			sb->append(L"Peptide Monoisotopic Mass" + StringHelper::toString(L'\t'));
			sb->append(L"Score" + StringHelper::toString(L'\t'));
			sb->append(L"Rank" + StringHelper::toString(L'\t'));

			sb->append(L"Matched Ion Series" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Mass-To-Charge Ratios" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Mass Diff (Ppm)" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Intensities" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Counts" + StringHelper::toString(L'\t'));
			sb->append(L"Decoy/Contaminant/Target" + StringHelper::toString(L'\t'));
			sb->append(L"QValue" + StringHelper::toString(L'\t'));

			delete sb;
			return sb->toString();
		}

		std::wstring CrosslinkSpectralMatch::GetTabSepHeaderGlyco()
		{
			auto sb = new StringBuilder();
			sb->append(L"File Name" + StringHelper::toString(L'\t'));
			sb->append(L"Scan Number" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor Scan Number" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor MZ" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor Charge" + StringHelper::toString(L'\t'));
			sb->append(L"Precursor Mass" + StringHelper::toString(L'\t'));
			sb->append(L"Cross Type" + StringHelper::toString(L'\t'));
			sb->append(std::wstring(L"Link Residues") + L"\t");

			sb->append(L"Peptide" + StringHelper::toString(L'\t'));
			sb->append(L"Protein Accession" + StringHelper::toString(L'\t'));
			sb->append(L"Protein Link Site" + StringHelper::toString(L'\t'));
			sb->append(L"Base Sequence" + StringHelper::toString(L'\t'));
			sb->append(L"Full Sequence" + StringHelper::toString(L'\t'));
			sb->append(L"Peptide Monoisotopic Mass" + StringHelper::toString(L'\t'));
			sb->append(L"Score" + StringHelper::toString(L'\t'));
			sb->append(L"Rank" + StringHelper::toString(L'\t'));

			sb->append(L"Matched Ion Series" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Mass-To-Charge Ratios" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Mass Diff (Ppm)" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Intensities" + StringHelper::toString(L'\t'));
			sb->append(L"Matched Ion Counts" + StringHelper::toString(L'\t'));

			sb->append(L"Decoy/Contaminant/Target" + StringHelper::toString(L'\t'));
			sb->append(L"QValue" + StringHelper::toString(L'\t'));

			sb->append(L"GlyID" + StringHelper::toString(L'\t'));
			sb->append(L"GlyMass" + StringHelper::toString(L'\t'));
			sb->append(L"GlyStruct(H,N,A,G,F)" + StringHelper::toString(L'\t'));

			delete sb;
			return sb->toString();
		}

		std::wstring CrosslinkSpectralMatch::ToString()
		{
			std::wstring position = L"";
			switch (getCrossType())
			{
				case PsmCrossType::Single:
					break;

				case PsmCrossType::Loop:
					position = L"(" + std::to_wstring(getLinkPositions()[0]) + L"-" + std::to_wstring(getLinkPositions()[1]) + L")";
					break;

				default:
					position = L"(" + std::to_wstring(getLinkPositions()[0]) + L")";
					break;
			}

			auto sb = new StringBuilder();
			sb->append(getFullFilePath() + L"\t");
			sb->append(std::to_wstring(getScanNumber()) + L"\t");
			sb->append(getPrecursorScanNumber() + L"\t");
			sb->append(std::to_wstring(getScanPrecursorMonoisotopicPeakMz()) + L"\t");
			sb->append(std::to_wstring(getScanPrecursorCharge()) + L"\t");
			sb->append(std::to_wstring(getScanPrecursorMass()) + L"\t");
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->append(getCrossType().ToString() + L"\t");



			if (getLinkPositions().size() > 0)
			{
				if (getCrossType() == PsmCrossType::Loop)
				{
					sb->append(getBaseSequence()[getLinkPositions()[0] - 1] + L";" + getBaseSequence()[getLinkPositions()[1] - 1] + L"\t");
				}
				else if (getCrossType() == PsmCrossType::Inter || getCrossType() == PsmCrossType::Intra)
				{
					sb->append(getBaseSequence()[getLinkPositions()[0] - 1] + L";" + getBetaPeptide()->getBaseSequence()[getBetaPeptide()->getLinkPositions()[0] - 1] + L"\t");
				}
				else
				{
					// deadend
					sb->append(getBaseSequence()[getLinkPositions()[0] - 1] + L"\t");
				}
			}
			else
			{
				sb->append(L"\t");
			}

			sb->append(L"\t");
			sb->append(getProteinAccession() + L"\t");
			sb->append(std::to_wstring(getXlProteinPos()) + L"\t");
			sb->append(getBaseSequence() + L"\t");
			sb->append(getFullSequence() + position + L"\t");
			sb->append((getPeptideMonisotopicMass().HasValue ? std::to_wstring(getPeptideMonisotopicMass().Value) : L"---"));
			sb->append(L"\t");
			sb->append(std::to_wstring(getScore()) + L"\t");
			sb->append(std::to_wstring(getXlRank()[0]) + L"\t");

			for (auto mid : MatchedIonDataDictionary(this))
			{
				sb->append(mid.Value);
				sb->append(L"\t");
			}

			if (getBetaPeptide() != nullptr)
			{
				sb->append(L"\t");
				sb->append(getBetaPeptide()->getProteinAccession() + L"\t");
				sb->append(std::to_wstring(getBetaPeptide()->getXlProteinPos()) + L"\t");
				sb->append(getBetaPeptide()->getBaseSequence() + L"\t");
				sb->append(getBetaPeptide()->getFullSequence() + L"(" + std::to_wstring(getBetaPeptide()->getLinkPositions()[0]) + L")" + L"\t");
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				sb->append(getBetaPeptide()->getPeptideMonisotopicMass().ToString() + L"\t");
				sb->append(std::to_wstring(getBetaPeptide()->getScore()) + L"\t");
				sb->append(std::to_wstring(getXlRank()[1]) + L"\t");

				for (auto betamid : MatchedIonDataDictionary(this->getBetaPeptide()))
				{
					sb->append(betamid.Value);
					sb->append(L"\t");
				}

				sb->append(L"\t");
				sb->append(std::to_wstring(getXLTotalScore()) + L"\t");

				// mass of crosslinker
				sb->append(((getPeptideMonisotopicMass().HasValue) ? std::to_wstring(getScanPrecursorMass() - getBetaPeptide()->getPeptideMonisotopicMass() - getPeptideMonisotopicMass().Value) : L"---"));
				sb->append(L"\t");

				int alphaNumParentIons = getMatchedFragmentIons().size()([&] (std::any p)
				{
				delete sb;
					return p::NeutralTheoreticalProduct->ProductType == ProductType::M;
				});
				int betaNumParentIons = getBetaPeptide()->getMatchedFragmentIons().size()([&] (std::any p)
				{
				delete sb;
					return p::NeutralTheoreticalProduct->ProductType == ProductType::M;
				});

				sb->append(std::to_wstring(alphaNumParentIons) + L";" + std::to_wstring(betaNumParentIons) + L"\t");
				sb->append(std::to_wstring(alphaNumParentIons) + std::to_wstring(betaNumParentIons) + L"\t");
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				sb->append(((getParentIonMaxIntensityRanks().size() > 0) && (getParentIonMaxIntensityRanks().Any()) ? getParentIonMaxIntensityRanks().Min().ToString() : L"-"));
				sb->append(L"\t");

			}

			if (getBetaPeptide() == nullptr)
			{
				sb->append((getIsDecoy()) ? L"D" : (getIsContaminant()) ? L"C" : L"T");
				sb->append(L"\t");
			}
			else
			{
				sb->append((getIsDecoy() || getBetaPeptide()->getIsDecoy()) ? L"D" : (getIsContaminant() || getBetaPeptide()->getIsContaminant()) ? L"C" : L"T");
				sb->append(L"\t");
			}

			sb->append(std::to_wstring(getFdrInfo()->getQValue()));
			sb->append(L"\t");


			delete sb;
			return sb->toString();
		}

		std::unordered_map<std::wstring, std::wstring> CrosslinkSpectralMatch::MatchedIonDataDictionary(PeptideSpectralMatch *psm)
		{
			std::unordered_map<std::wstring, std::wstring> s;
			AddMatchedIonsData(s, psm);
			return s;
		}
	}
}
