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

		std::string CrosslinkSpectralMatch::getParentIonExist() const
		{
			return privateParentIonExist;
		}

		void CrosslinkSpectralMatch::setParentIonExist(const std::string &value)
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

		std::vector<int> CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(std::vector<char> &crosslinkerModSites, PeptideWithSetModifications *peptide)
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

		std::string CrosslinkSpectralMatch::GetTabSepHeaderCross()
		{
			auto sb = new StringBuilder();
			sb->append("File Name" + StringHelper::toString(L'\t'));
			sb->append("Scan Number" + StringHelper::toString(L'\t'));
			sb->append("Precursor Scan Number" + StringHelper::toString(L'\t'));
			sb->append("Precursor MZ" + StringHelper::toString(L'\t'));
			sb->append("Precursor Charge" + StringHelper::toString(L'\t'));
			sb->append("Precursor Mass" + StringHelper::toString(L'\t'));
			sb->append("Cross Type" + StringHelper::toString(L'\t'));
			sb->append(std::string("Link Residues") + "\t");

			sb->append("Peptide" + StringHelper::toString(L'\t'));
			sb->append("Protein Accession" + StringHelper::toString(L'\t'));
			sb->append("Protein Link Site" + StringHelper::toString(L'\t'));
			sb->append("Base Sequence" + StringHelper::toString(L'\t'));
			sb->append("Full Sequence" + StringHelper::toString(L'\t'));
			sb->append("Peptide Monoisotopic Mass" + StringHelper::toString(L'\t'));
			sb->append("Score" + StringHelper::toString(L'\t'));
			sb->append("Rank" + StringHelper::toString(L'\t'));

			sb->append("Matched Ion Series" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Mass-To-Charge Ratios" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Mass Diff (Ppm)" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Intensities" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Counts" + StringHelper::toString(L'\t'));

			sb->append("Beta Peptide" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Protein Accession" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Protein LinkSite" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Base Sequence" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Full Sequence" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Theoretical Mass" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Score" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Rank" + StringHelper::toString(L'\t'));

			sb->append("Beta Peptide Matched Ion Series" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Matched Ion Mass-To-Charge Ratios" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Matched Ion Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Matched Ion Mass Diff (Ppm)" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Matched Ion Intensities" + StringHelper::toString(L'\t'));
			sb->append("Beta Peptide Matched Ion Counts" + StringHelper::toString(L'\t'));

			sb->append("Summary" + StringHelper::toString(L'\t'));
			sb->append("XL Total Score" + StringHelper::toString(L'\t'));
			sb->append("Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append("Parent Ions" + StringHelper::toString(L'\t'));
			sb->append("ParentIonsNum" + StringHelper::toString(L'\t'));
			sb->append("ParentIonMaxIntensityRank" + StringHelper::toString(L'\t'));
			sb->append("Decoy/Contaminant/Target" + StringHelper::toString(L'\t'));
			sb->append("QValue" + StringHelper::toString(L'\t'));


			delete sb;
			return sb->toString();
		}

		std::string CrosslinkSpectralMatch::GetTabSepHeaderSingle()
		{
			auto sb = new StringBuilder();
			sb->append("File Name" + StringHelper::toString(L'\t'));
			sb->append("Scan Number" + StringHelper::toString(L'\t'));
			sb->append("Precursor Scan Number" + StringHelper::toString(L'\t'));
			sb->append("Precursor MZ" + StringHelper::toString(L'\t'));
			sb->append("Precursor Charge" + StringHelper::toString(L'\t'));
			sb->append("Precursor Mass" + StringHelper::toString(L'\t'));
			sb->append("Cross Type" + StringHelper::toString(L'\t'));
			sb->append(std::string("Link Residues") + "\t");

			sb->append("Peptide" + StringHelper::toString(L'\t'));
			sb->append("Protein Accession" + StringHelper::toString(L'\t'));
			sb->append("Protein Link Site" + StringHelper::toString(L'\t'));
			sb->append("Base Sequence" + StringHelper::toString(L'\t'));
			sb->append("Full Sequence" + StringHelper::toString(L'\t'));
			sb->append("Peptide Monoisotopic Mass" + StringHelper::toString(L'\t'));
			sb->append("Score" + StringHelper::toString(L'\t'));
			sb->append("Rank" + StringHelper::toString(L'\t'));

			sb->append("Matched Ion Series" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Mass-To-Charge Ratios" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Mass Diff (Ppm)" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Intensities" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Counts" + StringHelper::toString(L'\t'));
			sb->append("Decoy/Contaminant/Target" + StringHelper::toString(L'\t'));
			sb->append("QValue" + StringHelper::toString(L'\t'));

			delete sb;
			return sb->toString();
		}

		std::string CrosslinkSpectralMatch::GetTabSepHeaderGlyco()
		{
			auto sb = new StringBuilder();
			sb->append("File Name" + StringHelper::toString(L'\t'));
			sb->append("Scan Number" + StringHelper::toString(L'\t'));
			sb->append("Precursor Scan Number" + StringHelper::toString(L'\t'));
			sb->append("Precursor MZ" + StringHelper::toString(L'\t'));
			sb->append("Precursor Charge" + StringHelper::toString(L'\t'));
			sb->append("Precursor Mass" + StringHelper::toString(L'\t'));
			sb->append("Cross Type" + StringHelper::toString(L'\t'));
			sb->append(std::string("Link Residues") + "\t");

			sb->append("Peptide" + StringHelper::toString(L'\t'));
			sb->append("Protein Accession" + StringHelper::toString(L'\t'));
			sb->append("Protein Link Site" + StringHelper::toString(L'\t'));
			sb->append("Base Sequence" + StringHelper::toString(L'\t'));
			sb->append("Full Sequence" + StringHelper::toString(L'\t'));
			sb->append("Peptide Monoisotopic Mass" + StringHelper::toString(L'\t'));
			sb->append("Score" + StringHelper::toString(L'\t'));
			sb->append("Rank" + StringHelper::toString(L'\t'));

			sb->append("Matched Ion Series" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Mass-To-Charge Ratios" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Mass Diff (Da)" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Mass Diff (Ppm)" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Intensities" + StringHelper::toString(L'\t'));
			sb->append("Matched Ion Counts" + StringHelper::toString(L'\t'));

			sb->append("Decoy/Contaminant/Target" + StringHelper::toString(L'\t'));
			sb->append("QValue" + StringHelper::toString(L'\t'));

			sb->append("GlyID" + StringHelper::toString(L'\t'));
			sb->append("GlyMass" + StringHelper::toString(L'\t'));
			sb->append("GlyStruct(H,N,A,G,F)" + StringHelper::toString(L'\t'));

			delete sb;
			return sb->toString();
		}

		std::string CrosslinkSpectralMatch::ToString()
		{
			std::string position = "";
			switch (getCrossType())
			{
				case PsmCrossType::Single:
					break;

				case PsmCrossType::Loop:
					position = "(" + std::to_string(getLinkPositions()[0]) + "-" + std::to_string(getLinkPositions()[1]) + ")";
					break;

				default:
					position = "(" + std::to_string(getLinkPositions()[0]) + ")";
					break;
			}

			auto sb = new StringBuilder();
			sb->append(getFullFilePath() + "\t");
			sb->append(std::to_string(getScanNumber()) + "\t");
			sb->append(getPrecursorScanNumber() + "\t");
			sb->append(std::to_string(getScanPrecursorMonoisotopicPeakMz()) + "\t");
			sb->append(std::to_string(getScanPrecursorCharge()) + "\t");
			sb->append(std::to_string(getScanPrecursorMass()) + "\t");
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->append(getCrossType().ToString() + "\t");



			if (getLinkPositions().size() > 0)
			{
				if (getCrossType() == PsmCrossType::Loop)
				{
					sb->append(getBaseSequence()[getLinkPositions()[0] - 1] + ";" + getBaseSequence()[getLinkPositions()[1] - 1] + "\t");
				}
				else if (getCrossType() == PsmCrossType::Inter || getCrossType() == PsmCrossType::Intra)
				{
					sb->append(getBaseSequence()[getLinkPositions()[0] - 1] + ";" + getBetaPeptide()->getBaseSequence()[getBetaPeptide()->getLinkPositions()[0] - 1] + "\t");
				}
				else
				{
					// deadend
					sb->append(getBaseSequence()[getLinkPositions()[0] - 1] + "\t");
				}
			}
			else
			{
				sb->append("\t");
			}

			sb->append("\t");
			sb->append(getProteinAccession() + "\t");
			sb->append(std::to_string(getXlProteinPos()) + "\t");
			sb->append(getBaseSequence() + "\t");
			sb->append(getFullSequence() + position + "\t");
			sb->append((getPeptideMonisotopicMass().HasValue ? std::to_string(getPeptideMonisotopicMass().Value) : "---"));
			sb->append("\t");
			sb->append(std::to_string(getScore()) + "\t");
			sb->append(std::to_string(getXlRank()[0]) + "\t");

			for (auto mid : MatchedIonDataDictionary(this))
			{
				sb->append(mid.Value);
				sb->append("\t");
			}

			if (getBetaPeptide() != nullptr)
			{
				sb->append("\t");
				sb->append(getBetaPeptide()->getProteinAccession() + "\t");
				sb->append(std::to_string(getBetaPeptide()->getXlProteinPos()) + "\t");
				sb->append(getBetaPeptide()->getBaseSequence() + "\t");
				sb->append(getBetaPeptide()->getFullSequence() + "(" + std::to_string(getBetaPeptide()->getLinkPositions()[0]) + ")" + "\t");
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				sb->append(getBetaPeptide()->getPeptideMonisotopicMass().ToString() + "\t");
				sb->append(std::to_string(getBetaPeptide()->getScore()) + "\t");
				sb->append(std::to_string(getXlRank()[1]) + "\t");

				for (auto betamid : MatchedIonDataDictionary(this->getBetaPeptide()))
				{
					sb->append(betamid.Value);
					sb->append("\t");
				}

				sb->append("\t");
				sb->append(std::to_string(getXLTotalScore()) + "\t");

				// mass of crosslinker
				sb->append(((getPeptideMonisotopicMass().HasValue) ? std::to_string(getScanPrecursorMass() - getBetaPeptide()->getPeptideMonisotopicMass() - getPeptideMonisotopicMass().Value) : "---"));
				sb->append("\t");

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

				sb->append(std::to_string(alphaNumParentIons) + ";" + std::to_string(betaNumParentIons) + "\t");
				sb->append(std::to_string(alphaNumParentIons) + std::to_string(betaNumParentIons) + "\t");
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				sb->append(((getParentIonMaxIntensityRanks().size() > 0) && (getParentIonMaxIntensityRanks().Any()) ? getParentIonMaxIntensityRanks().Min().ToString() : "-"));
				sb->append("\t");

			}

			if (getBetaPeptide() == nullptr)
			{
				sb->append((getIsDecoy()) ? "D" : (getIsContaminant()) ? "C" : "T");
				sb->append("\t");
			}
			else
			{
				sb->append((getIsDecoy() || getBetaPeptide()->getIsDecoy()) ? "D" : (getIsContaminant() || getBetaPeptide()->getIsContaminant()) ? "C" : "T");
				sb->append("\t");
			}

			sb->append(std::to_string(getFdrInfo()->getQValue()));
			sb->append("\t");


			delete sb;
			return sb->toString();
		}

		std::unordered_map<std::string, std::string> CrosslinkSpectralMatch::MatchedIonDataDictionary(PeptideSpectralMatch *psm)
		{
			std::unordered_map<std::string, std::string> s;
			AddMatchedIonsData(s, psm);
			return s;
		}
	}
}
