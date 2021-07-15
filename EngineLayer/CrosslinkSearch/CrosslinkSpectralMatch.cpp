#include "CrosslinkSpectralMatch.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "Crosslinker.h"

#include "Sort.h"

#include <numeric>
#include <algorithm>
#include <sstream>
#include <iomanip>

using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
    namespace CrosslinkSearch
    {
        
        CrosslinkSpectralMatch::CrosslinkSpectralMatch(PeptideWithSetModifications *theBestPeptide, int notch,
                                                       double score, int scanIndex,
                                                       Ms2ScanWithSpecificMass *scan,
                                                       DigestionParams *digestionParams,
                                                       std::vector<MatchedFragmentIon*> &matchedFragmentIons) :
            PeptideSpectralMatch(theBestPeptide, notch, score, scanIndex, scan, digestionParams, matchedFragmentIons)
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
        
        std::vector<int> CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(std::vector<char> &crosslinkerModSites,
                                                                                PeptideWithSetModifications *peptide)
        {
            std::vector<int> possibleXlPositions;

            bool wildcard =false;
            for ( char p: crosslinkerModSites ){
                if (  p == 'X' ) {
                    wildcard = true;
                    break;
                }                      
            }
            
            for (int r = 0; r < (int) peptide->getBaseSequence().size(); r++)
            {
                //if (crosslinkerModSites.Contains(peptide->getBaseSequence()[r]) || wildcard)
                if ( std::find(crosslinkerModSites.begin(), crosslinkerModSites.end(), peptide->getBaseSequence()[r]) !=
                     crosslinkerModSites.end() || wildcard ) 
                {
                    possibleXlPositions.push_back(r + 1);
                }
            }
            
            return possibleXlPositions;
        }
        
        std::vector<int> CrosslinkSpectralMatch::GenerateIntensityRanks(std::vector<double> &experimental_intensities)
        {
            auto y = experimental_intensities;
#ifdef ORIG
            auto x = Enumerable::Range(1, y.size()).OrderBy([&] (std::any p) {
                    return p;
                })->ToArray();
            Array::Sort(y, x);
#endif
            std::vector<int> x(y.size());
            std::iota(x.begin(), x.end(), 1);
            Sort::SortPairs(y, x, y.size() );
            
#ifdef ORIG
            auto experimental_intensities_rank = Enumerable::Range(1, y.size()).OrderByDescending([&] (std::any p)  {
                    return p;
                })->ToArray();
            Array::Sort(x, experimental_intensities_rank);
#endif
            std::vector<int> experimental_intensities_rank(y.size() );
            int n = y.size();
            std::generate(experimental_intensities_rank.begin(), experimental_intensities_rank.end(), [&] () {return n--;});
            Sort::SortPairs(x, experimental_intensities_rank, x.size() );
            
            return experimental_intensities_rank;
        }
        
        std::string CrosslinkSpectralMatch::GetTabSepHeaderCross()
        {
            auto sb = new StringBuilder();
            sb->append("File Name\t");
            sb->append("Scan Number\t");
            sb->append("Precursor Scan Number\t");
            sb->append("Precursor MZ\t");
            sb->append("Precursor Charge\t");
            sb->append("Precursor Mass\t");
            sb->append("Cross Type\t");
            sb->append("Link Residues\t");
            
            sb->append("Peptide\t");
            sb->append("Protein Accession\t");
            sb->append("Protein Link Site\t");
            sb->append("Base Sequence\t");
            sb->append("Full Sequence\t");
            sb->append("Peptide Monoisotopic Mass\t");
            sb->append("Score\t");
            sb->append("Rank\t");
            
            sb->append("Matched Ion Series\t");
            sb->append("Matched Ion Mass-To-Charge Ratios\t");
            sb->append("Matched Ion Mass Diff (Da)\t");
            sb->append("Matched Ion Mass Diff (Ppm)\t");
            sb->append("Matched Ion Intensities\t");
            sb->append("Matched Ion Counts\t");
            
            sb->append("Beta Peptide\t");
            sb->append("Beta Peptide Protein Accession\t");
            sb->append("Beta Peptide Protein LinkSite\t");
            sb->append("Beta Peptide Base Sequence\t");
            sb->append("Beta Peptide Full Sequence\t");
            sb->append("Beta Peptide Theoretical Mass\t");
            sb->append("Beta Peptide Score\t");
            sb->append("Beta Peptide Rank\t");
            
            sb->append("Beta Peptide Matched Ion Series\t");
            sb->append("Beta Peptide Matched Ion Mass-To-Charge Ratios\t");
            sb->append("Beta Peptide Matched Ion Mass Diff (Da)\t");
            sb->append("Beta Peptide Matched Ion Mass Diff (Ppm)\t");
            sb->append("Beta Peptide Matched Ion Intensities\t");
            sb->append("Beta Peptide Matched Ion Counts\t");
            
            sb->append("Summary\t");
            sb->append("XL Total Score\t");
            sb->append("Mass Diff (Da)\t");
            sb->append("Parent Ions\t");
            sb->append("ParentIonsNum\t");
            sb->append("ParentIonMaxIntensityRank\t");
            sb->append("Decoy/Contaminant/Target\t");
            sb->append("QValue\t");
            
            
            std::string s = sb->toString();
            delete sb;
            return s;
        }
        
        std::string CrosslinkSpectralMatch::GetTabSepHeaderSingle()
        {
            auto sb = new StringBuilder();
            sb->append("File Name\t");
            sb->append("Scan Number\t");
            sb->append("Precursor Scan Number\t");
            sb->append("Precursor MZ\t");
            sb->append("Precursor Charge\t");
            sb->append("Precursor Mass\t");
            sb->append("Cross Type\t");
            sb->append("Link Residues\t");

            sb->append("Protein Accession\t");
            sb->append("Protein Link Site\t");
            sb->append("Base Sequence\t");
            sb->append("Full Sequence\t");
            sb->append("Peptide Monoisotopic Mass\t");
            sb->append("Score\t");
            sb->append("Rank\t");
            
            sb->append("Matched Ion Series\t");
            sb->append("Matched Ion Mass-To-Charge Ratios\t");
            sb->append("Matched Ion Mass Diff (Da)\t");
            sb->append("Matched Ion Mass Diff (Ppm)\t");
            sb->append("Matched Ion Intensities\t");
            sb->append("Matched Ion Counts\t");
            sb->append("Decoy/Contaminant/Target\t");
            sb->append("QValue\t");
            
            std::string s = sb->toString();
            delete sb;
            return s;
        }
        
        std::string CrosslinkSpectralMatch::GetTabSepHeaderGlyco()
        {
            auto sb = new StringBuilder();
            sb->append("File Name\t");
            sb->append("Scan Number\t");
            sb->append("Precursor Scan Number\t");
            sb->append("Precursor MZ\t");
            sb->append("Precursor Charge\t");
            sb->append("Precursor Mass\t");
            sb->append("Cross Type\t");
            sb->append("Link Residues\t");

            sb->append("Protein Accession\t");
            sb->append("Protein Link Site\t");
            sb->append("Base Sequence\t");
            sb->append("Full Sequence\t");
            sb->append("Peptide Monoisotopic Mass\t");
            sb->append("Score\t");
            sb->append("Rank\t");
            
            sb->append("Matched Ion Series\t");
            sb->append("Matched Ion Mass-To-Charge Ratios\t");
            sb->append("Matched Ion Mass Diff (Da)\t");
            sb->append("Matched Ion Mass Diff (Ppm)\t");
            sb->append("Matched Ion Intensities\t");
            sb->append("Matched Ion Counts\t");
            
            sb->append("Decoy/Contaminant/Target\t");
            sb->append("QValue\t");
            
            sb->append("GlyID\t");
            sb->append("GlyMass\t");
            sb->append("GlyStruct(H,N,A,G,F)\t");
            
            std::string s = sb->toString();
            delete sb;
            return s;
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
            if ( getPrecursorScanNumber().has_value() ) {
                sb->append(std::to_string(getPrecursorScanNumber().value()) + "\t");
            }
            else {
                std::string s = "-\t";
                sb->append(s);
            }
            sb->append(std::to_string(getScanPrecursorMonoisotopicPeakMz()) + "\t");
            sb->append(std::to_string(getScanPrecursorCharge()) + "\t");
            sb->append(std::to_string(getScanPrecursorMass()) + "\t");
            auto crosslinktype = getCrossType();
            sb->append(PsmCrossTypeToString(crosslinktype) + "\t");
            
            if (getLinkPositions().size() > 0)
            {
                if (getCrossType() == PsmCrossType::Loop)
                {
                    std::stringstream ss;
                    ss << getBaseSequence()[getLinkPositions()[0] - 1] << ';' <<
                        getBaseSequence()[getLinkPositions()[1] - 1] << '\t';
                    sb->append(ss.str() );
                }
                else if (getCrossType() == PsmCrossType::Inter || getCrossType() == PsmCrossType::Intra)
                {
                    std::stringstream ss;
                    ss << getBaseSequence()[getLinkPositions()[0] - 1] << ';' <<
                        getBetaPeptide()->getBaseSequence()[getBetaPeptide()->getLinkPositions()[0] - 1] << '\t';
                    sb->append(ss.str());
                }
                else
                {
                    // deadend
                    std::stringstream ss;
                    ss << getBaseSequence()[getLinkPositions()[0] - 1] << "\t";
                    sb->append(ss.str() );
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
            sb->append((getPeptideMonisotopicMass().has_value() ? std::to_string(getPeptideMonisotopicMass().value()) : "---"));
            sb->append("\t");
            sb->append(std::to_string(getScore()) + "\t");
            sb->append(std::to_string(getXlRank()[0]) + "\t");
            
            for (auto mid : MatchedIonDataDictionary(this))
            {
                sb->append(std::get<1>(mid));
                sb->append("\t");
            }
            
            if (getBetaPeptide() != nullptr)
            {
                auto betaPeptide = getBetaPeptide();
                
                sb->append("\t");
                sb->append(betaPeptide->getProteinAccession() + "\t");
                sb->append(std::to_string(betaPeptide->getXlProteinPos()) + "\t");
                sb->append(betaPeptide->getBaseSequence() + "\t");
                sb->append(betaPeptide->getFullSequence() + "(" + std::to_string(betaPeptide->getLinkPositions()[0]) +
                           ")" + "\t");
                sb->append(std::to_string(betaPeptide->getPeptideMonisotopicMass().value()) + "\t");
                sb->append(std::to_string(betaPeptide->getScore()) + "\t");
                sb->append(std::to_string(getXlRank()[1]) + "\t");
                
                for (auto betamid : MatchedIonDataDictionary(this->getBetaPeptide()))
                {
                    sb->append(std::get<1>(betamid));
                    sb->append("\t");
                }
                
                sb->append("\t");
                sb->append(std::to_string(getXLTotalScore()) + "\t");
                
                // mass of crosslinker
                sb->append(((getPeptideMonisotopicMass().has_value()) ? std::to_string(getScanPrecursorMass() -
                            betaPeptide->getPeptideMonisotopicMass().value() - getPeptideMonisotopicMass().value()) : "---"));
                sb->append("\t");
                
                int alphaNumParentIons = 0;
                for ( auto p : getMatchedFragmentIons() ) {
                    if ( p->NeutralTheoreticalProduct->productType == ProductType::M ) {
                        alphaNumParentIons++;
                    }
                }

                int betaNumParentIons = 0;
                for ( auto p :  betaPeptide->getMatchedFragmentIons() ) {
                    if ( p->NeutralTheoreticalProduct->productType == ProductType::M ) {
                        alphaNumParentIons++;
                    }
                }

                
                sb->append(std::to_string(alphaNumParentIons) + ";" + std::to_string(betaNumParentIons) + "\t");
                sb->append(std::to_string(alphaNumParentIons) + std::to_string(betaNumParentIons) + "\t");
                sb->append(((getParentIonMaxIntensityRanks().size() > 0) && (!getParentIonMaxIntensityRanks().empty()) ?
                   std::to_string(*std::min_element(getParentIonMaxIntensityRanks().begin(), getParentIonMaxIntensityRanks().end()))
                            : "-"));
                sb->append("\t");                            
            }
            
            if (getBetaPeptide() == nullptr)
            {
                sb->append((getIsDecoy()) ? "D" : (getIsContaminant()) ? "C" : "T");
                sb->append("\t");
            }
            else
            {
                sb->append((getIsDecoy() || getBetaPeptide()->getIsDecoy()) ? "D" :
                           (getIsContaminant() || getBetaPeptide()->getIsContaminant()) ? "C" : "T");
                sb->append("\t");
            }
            
            sb->append(std::to_string(getFdrInfo()->getQValue()));
            sb->append("\t");
            
            
            std::string s= sb->toString();
            delete sb;
            return s;
        }
        
        std::vector<std::tuple<std::string, std::string>> CrosslinkSpectralMatch::MatchedIonDataDictionary(PeptideSpectralMatch *psm)
        {
            std::vector<std::tuple<std::string, std::string>> s;
            AddMatchedIonsData(s, psm);
            return s;
        }


        int CrosslinkSpectralMatch::Pack(char *buf, size_t &buf_len,
                                         const std::vector<CrosslinkSpectralMatch *> &csmVec)
        {
            size_t pos = 0;
            int ret;

            for ( auto csm: csmVec ) {
                size_t len = buf_len - pos;
                ret = CrosslinkSpectralMatch::Pack_internal(buf+pos, len, csm);
                if ( ret == -1 ) {
                    buf_len = pos + len;
                    return ret;
                }
                pos += ret;
                auto betaPeptide = csm->getBetaPeptide();
                if ( betaPeptide != nullptr ) {
                    len = buf_len - pos;
                    ret = CrosslinkSpectralMatch::Pack_internal(buf+pos, len, betaPeptide);
                    if ( ret == -1 ) {
                        buf_len = pos + len;
                        return ret;
                    }
                    pos += ret;                    
                }
            }
            buf_len = pos;
            return pos;
        }

        int CrosslinkSpectralMatch::Pack(char *buf, size_t &buf_len, CrosslinkSpectralMatch *csm)
        {
            size_t pos = 0;
            int ret;

            size_t len = buf_len - pos;
            ret = CrosslinkSpectralMatch::Pack_internal(buf+pos, len, csm);
            if ( ret == -1 ) {
                buf_len = pos + len;
                return ret;
            }
            pos += ret;
            auto betaPeptide = csm->getBetaPeptide();
            if ( betaPeptide != nullptr ) {
                len = buf_len - pos;
                ret = CrosslinkSpectralMatch::Pack_internal(buf+pos, len, betaPeptide);
                if ( ret == -1 ) {
                    buf_len = pos + len;
                    return ret;
                }
                pos += ret;                    
            }
            buf_len = pos;
            return pos;            
        }
        
        int CrosslinkSpectralMatch::Pack_internal(char *buf, size_t &buf_len, CrosslinkSpectralMatch *csm)
        {
            size_t bufpos = 0;
            size_t pos = BinaryPack::LineStartOffset;

            auto mFrIons = csm->getMatchedFragmentIons ();
            auto dp = csm->digestionParams;
            auto uMapPep = csm->getPeptidesToMatchingFragments();
            std::vector<int> lPositions  = csm->getLinkPositions();
            std::vector<int> xlRanks = csm->getXlRank();
            bool has_beta_peptide = csm->getBetaPeptide() != nullptr;          
            
            // line 1
            retlen = BinaryPack::PackBool(tmpbuf+pos, csm->getNotch().has_value());
            pos += retlen;
            if ( csm->getNotch().has_value() ) {
                retlen = BinaryPack::PackInt(tmpbuf+pos, csm->getNotch().value());
                pos += retlen;
            }
            
            //output << std::setprecision(12) << csm->getXLTotalScore() << "\t" <<
            //    csm->getDeltaScore() << "\t" <<
            //    csm->getScanIndex() << "\t" <<
            //    csm->getScanNumber() << "\t" <<
            //    csm->getXlProteinPos() << "\t" <<
            //    csm->getScore() << "\t" <<
            //    csm->getRunnerUpScore() << "\t";
            retlen = BinaryPack::PackDouble(tmpbuf+pos, csm->getXLTotalScore() );
            pos += retlen;
            retlen = BinaryPack::PackDouble(tmpbuf+pos, csm->getDeltaScore() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, csm->getScanIndex() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, csm->getScanNumber() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, csm->getXlProteinPos() );
            pos += retlen;
            retlen = BinaryPack::PackDouble(tmpbuf+pos, csm->getScore() );
            pos += retlen;
            retlen = BinaryPack::PackDouble(tmpbuf+pos, csm->getRunnerUpScore() );
            pos += retlen;
            
            PsmCrossType ctype = csm->getCrossType();
            //output << PsmCrossTypeToString(ctype) << "\t" <<
            //    csm->getMatchedFragmentIons().size() << "\t" << 
            //    lPositions.size() << "\t" <<
            //    xlRanks.size() << "\t";
            retlen = BinaryPack::PackString(tmpbuf+pos, PsmCrossTypeToString(ctype) );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, csm->getMatchedFragmentIons().size() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, lPositions.size() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, xlRanks.size() );
            pos += retlen;
                        
            //if ( has_beta_peptide ) {
            //    output << "true" << std::endl;
            //}
            //else {
            //    output << "false" << std::endl;
            //}
            retlen = BinaryPack::PackBool(tmpbuf+pos, has_beta_peptide );
            pos += retlen;
            
            // set lenght of line right at the beginning.
            BinaryPack::SetLineLength(tmpbuf, pos);

            if ( bufpos+pos > buf_len ) {
                buf_len = bufpos+pos;
                return -1;
            }
            memcpy ( buf+bufpos, tmpbuf, pos );
            bufpos += pos;

            //Line 2: FdrInfo related Data
            pos = BinaryPack::LineStartOffset;
            memset(tmpbuf, 0, 256);
            
            FdrInfo *fdr = csm->getFdrInfo();
            retlen = BinaryPack::PackBool(tmpbuf+pos, (fdr != nullptr) );
            pos += retlen;
            
            if ( fdr != nullptr ) {
                //output <<  fdr->getCumulativeTarget() <<  "\t" <<
                //    fdr->getCumulativeDecoy() <<  "\t" <<
                //    fdr->getQValue() <<  "\t" <<
                //    fdr->getCumulativeTargetNotch() <<  "\t" <<
                //    fdr->getCumulativeDecoyNotch() <<  "\t" <<
                //    fdr->getQValueNotch() <<  "\t" <<
                //    fdr->getMaximumLikelihood() <<  "\t" <<
                //    fdr->getEValue() <<  "\t" <<
                //    fdr->getEScore() <<  "\t";
                //if ( fdr->getCalculateEValue() ) {
                //    output << "true\n";
                // }
                //else {
                //    output << "false\n";
                //}
                retlen = BinaryPack::PackDouble(tmpbuf+pos, fdr->getCumulativeTarget() );
                pos += retlen;
                retlen = BinaryPack::PackDouble(tmpbuf+pos, fdr->getCumulativeDecoy() );
                pos += retlen;
                retlen = BinaryPack::PackDouble(tmpbuf+pos, fdr->getQValue() );
                pos += retlen;
                retlen = BinaryPack::PackDouble(tmpbuf+pos, fdr->getCumulativeTargetNotch() );
                pos += retlen;
                retlen = BinaryPack::PackDouble(tmpbuf+pos, fdr->getCumulativeDecoyNotch() );
                pos += retlen;
                retlen = BinaryPack::PackDouble(tmpbuf+pos, fdr->getQValueNotch() );
                pos += retlen;
                retlen = BinaryPack::PackDouble(tmpbuf+pos, fdr->getMaximumLikelihood() );
                pos += retlen;
                retlen = BinaryPack::PackDouble(tmpbuf+pos, fdr->getEValue() );
                pos += retlen;
                retlen = BinaryPack::PackDouble(tmpbuf+pos, fdr->getEScore() );
                pos += retlen;
                retlen = BinaryPack::PackBool(tmpbuf+pos, fdr->getCalculateEValue() );
                pos += retlen;                
            }
            // set lenght of line right at the beginning.
            BinaryPack::SetLineLength(tmpbuf, pos);

            if ( bufpos+pos > buf_len ) {
                buf_len = bufpos+pos;
                return -1;
            }
            memcpy ( buf+bufpos, tmpbuf, pos );
            bufpos += pos;
                        
            //line 3: LinkPositions
            pos = BinaryPack::LineStartOffset;
            memset(tmpbuf, 0, 256);

            for ( auto lp: lPositions ) {
                //output << lp << "\t";
                retlen = BinaryPack::PackInt(tmpbuf+pos, lp );
                pos += retlen;
            }
            // set lenght of line right at the beginning.
            BinaryPack::SetLineLength(tmpbuf, pos);

            if ( bufpos+pos > buf_len ) {
                buf_len = bufpos+pos;
                return -1;
            }
            memcpy ( buf+bufpos, tmpbuf, pos );
            bufpos += pos;
            
            //line 4: xlRank
            pos = BinaryPack::LineStartOffset;
            memset(tmpbuf, 0, 256);

            for ( auto xl: xlRanks) {
                //output << xl << "\t";
                retlen = BinaryPack::PackInt(tmpbuf+pos, xl );
                pos += retlen;
            }

            // set lenght of line right at the beginning.
            BinaryPack::SetLineLength(tmpbuf, pos);

            if ( bufpos+pos > buf_len ) {
                buf_len = bufpos+pos;
                return -1;
            }
            memcpy ( buf+bufpos, tmpbuf, pos );
            bufpos += pos;
            
            //line 5: DigestionParams();
            pos = BinaryPack::LineStartOffset;
            memset(tmpbuf, 0, 256);
            
            std::string s = dp->ToString();
            retlen = BinaryPack::PackString(tmpbuf+pos, s);
            pos += retlen;

            // set lenght of line right at the beginning.
            BinaryPack::SetLineLength(tmpbuf, pos);

            if ( bufpos+pos > buf_len ) {
                buf_len = bufpos+pos;
                return -1;
            }
            memcpy ( buf+bufpos, tmpbuf, pos );
            bufpos += pos;
                        
            //line 6-10: PeptideWithSetModifications;
            //Assuming right now only a single PeptideWithSetModifications
            if ( uMapPep.size() != 1 ) {
                std::cout << "CrosslinkSpectralMatch::Pack: Error - unordered_map has more than one entry!\n";
            }
            auto pep = std::get<0>(*uMapPep.begin());
            size_t tmp_len = buf_len - bufpos;
            // need to check whether the routine below adds a \n at the end!
            int ret = PeptideWithSetModifications::Pack(buf+bufpos, tmp_len, pep );
            if ( ret == -1 ) {
                buf_len += tmp_len - (buf_len - bufpos);
                return -1;
            }
            bufpos += tmp_len;
            
            //line 11-x: one line for each MatchedFragmentIon
            for ( auto i=0; i< mFrIons.size(); i++ ) {
                tmp_len = buf_len - bufpos;
                // need to check whether the routine below adds a \n at the end!
                ret = MatchedFragmentIon::Pack ( buf+bufpos, tmp_len, mFrIons[i]);
                if ( ret == -1 ) {
                    buf_len += tmp_len - (buf_len - bufpos);
                    return -1;
                }
                bufpos += tmp_len;
            }
            
            return (int)bufpos;
        }

        void CrosslinkSpectralMatch::Unpack (char *buf, size_t buf_len, int count, size_t &len,
                                             std::vector<CrosslinkSpectralMatch*> &pepVec,
                                             const std::vector<Ms2ScanWithSpecificMass*> &ms2Scans,
                                             const std::vector<Protein *> &proteinList )
        {
            std::string input_buf (buf);
            std::vector<std::string> lines = StringHelper::split(input_buf, '\n');

            size_t total_len=0;
            int counter=0;
            for (auto  i=0; i < lines.size();  ) {
                size_t tmp_len=0;
                CrosslinkSpectralMatch *pep;
                bool has_beta_peptide=false;
                CrosslinkSpectralMatch::Unpack_internal(lines, i, tmp_len, &pep, ms2Scans, proteinList,
                                                        has_beta_peptide );
                total_len += tmp_len;
                pepVec.push_back(pep);
                if ( has_beta_peptide ) {
                    CrosslinkSpectralMatch *beta_pep;
                    CrosslinkSpectralMatch::Unpack_internal(lines, i, tmp_len, &beta_pep, ms2Scans, proteinList,
                                                            has_beta_peptide );
                    pep->setBetaPeptide(beta_pep);
                    total_len += tmp_len;
                }
                counter ++;
                if ( counter == count ) break;
            }
            len = total_len;
        }

        void CrosslinkSpectralMatch::Unpack (char *buf, size_t buf_len, size_t &len,
                                             CrosslinkSpectralMatch** newCsm,
                                             const std::vector<Ms2ScanWithSpecificMass*> &ms2Scans,
                                             const std::vector<Protein *> &proteinList )
        {
            std::string input(buf);
            std::vector<std::string> lines = StringHelper::split(input, '\n');
            int index=0;
            if ( lines.size() < 10 ) {
                std::cout << "CrosslinkSpectralMatch::Unpack : input does not contains enough information to " <<
                    "reconstruct the CrosslinkSpectralMatch. " << std::endl;
                return;
            }
            bool has_beta_peptide=false;            
            CrosslinkSpectralMatch::Unpack_internal ( lines, index, len, newCsm, ms2Scans, proteinList,
                                                      has_beta_peptide );
            if ( has_beta_peptide) {
                CrosslinkSpectralMatch* beta_pep;
                size_t tmp_len=0;
                CrosslinkSpectralMatch::Unpack_internal ( lines, index, tmp_len, &beta_pep, ms2Scans, proteinList,
                                                          has_beta_peptide );
                (*newCsm)->setBetaPeptide(beta_pep);
                len += tmp_len;
            }            
        }

        void CrosslinkSpectralMatch::Unpack_internal (std::vector<std::string> &input,
                                                      int &index, size_t &len,
                                                      CrosslinkSpectralMatch** newCsm,
                                                      const std::vector<Ms2ScanWithSpecificMass*> &ms2Scans,
                                                      const std::vector<Protein *> &proteinList,
                                                      bool &has_beta_peptide )
        {
            size_t total_len = 0;

            //Dissect line 1: generic information
            std::vector<std::string> splits = StringHelper::split(input[index], '\t');
            //+1 for the endl symbol removed when converting the char* into a std::vector<std::string>
            total_len += input[index].length() + 1; 
            index++;
            int notch=-1, scanindex, scannumber, proteinPos, matchedFragmentIonsVecsize, lpositionsize, xlranksize;
            double  deltaScore, XLTotalScore, score, runnerUpScore;

            if ( splits[0] == "true" ) {
                notch = std::stoi(splits[1]);
            }
            XLTotalScore = std::stod (splits[2]);
            deltaScore   = std::stod (splits[3]);
            scanindex    = std::stoi (splits[4]);
            scannumber   = std::stoi (splits[5]);
            proteinPos   = std::stoi (splits[6]);
            score        = std::stod (splits[7]);
            runnerUpScore = std::stod(splits[8]);
            PsmCrossType ctype = PsmCrossTypeFromString(splits[9]);
            matchedFragmentIonsVecsize = std::stoi(splits[10]);
            lpositionsize = std::stoi(splits[11]);
            xlranksize     = std::stoi(splits[12]);
            if ( splits[13] == "true" ) {
                has_beta_peptide = true;
            }

            //line 2: FdrInfo related data
            double cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetNotch;
            double cumulativeDecoyNotch, qValueNotch, maximumLikelihood, eValue, eScore;
            bool calculateEValue=false;
            bool have_fdr = false;
            
            splits.clear();
            splits = StringHelper::split(input[index], '\t');
            total_len += input[index].length() + 1;            
            index++;
            if ( splits.size() > 1 ) {
                have_fdr = true;
                
                cumulativeTarget = std::stod(splits[0]);
                cumulativeDecoy  = std::stod(splits[1]);
                qValue           = std::stod(splits[2]);
                cumulativeTargetNotch = std::stod(splits[3]);
                cumulativeDecoyNotch  = std::stod(splits[4]);
                qValueNotch           = std::stod(splits[5]);
                maximumLikelihood     = std::stod(splits[6]);
                eValue                = std::stod(splits[7]);
                eScore                = std::stod(splits[8]);
                if ( splits[9] == "true" ) {
                    calculateEValue = true;
                }
            }

            //line 3: linkPositions
            std::vector<int> linkPosvec;
            splits.clear();
            splits = StringHelper::split(input[index], '\t');
            total_len += input[index].length() + 1; 
            index++;
            for ( auto i=0; i<lpositionsize; i++ ) {
                linkPosvec.push_back(std::stoi(splits[i]));
            }
            
            //line 4: xlRank
            std::vector<int> xlRankVec;
            splits.clear();
            splits = StringHelper::split(input[index], '\t');
            total_len += input[index].length() + 1; 
            index++;
            for ( auto i=0; i<xlranksize; i++ ) {
                xlRankVec.push_back(std::stoi(splits[i]));
            }
            
            //line 5: DigestionParams
            DigestionParams *dp = DigestionParams::FromString(input[index]);
            total_len += input[index].length() + 1; 
            index++;

            //line 6-10: PeptideWithSetModifications
            PeptideWithSetModifications* pep;
            size_t tmp_len=0;
            PeptideWithSetModifications::Unpack(input, index, tmp_len, &pep);
            pep->SetNonSerializedPeptideInfo ( proteinList );
            total_len += tmp_len;
            index += 4;

            
            // line 11-x: Vector of MatchedFragmentIons
            std::vector<MatchedFragmentIon*> matchedFragmentIonsVec;
            for ( auto i=0; i< matchedFragmentIonsVecsize; i++ ) {
                MatchedFragmentIon *ion;
                tmp_len=0;
                MatchedFragmentIon::Unpack(input[index], tmp_len, &ion);
                matchedFragmentIonsVec.push_back(ion);
                index++;
                total_len += tmp_len;                    
            }

            Ms2ScanWithSpecificMass *scan = nullptr;
            for ( auto i = 0; i < ms2Scans.size() ; i ++ ) {
                if ( ms2Scans[i]->getOneBasedScanNumber() == scannumber ) {
                    scan = ms2Scans[i];
                    break;
                }
            }
            CrosslinkSpectralMatch *csm = new CrosslinkSpectralMatch ( pep, notch, XLTotalScore, scanindex, scan, dp,
                                                                       matchedFragmentIonsVec );
            csm->setXLTotalScore(XLTotalScore);
            csm->setDeltaScore(deltaScore);
            csm->setXlProteinPos(proteinPos);
            csm->setScore(score);
            csm->setRunnerUpScore(runnerUpScore);
            csm->setCrossType (ctype);
            csm->setXlRank(xlRankVec);
            csm->setLinkPositions(linkPosvec);

            if ( have_fdr) {
                csm->SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetNotch,
                                  cumulativeDecoyNotch, qValueNotch, maximumLikelihood, 
                                  eValue, eScore, calculateEValue);
            }
            csm->ResolveAllAmbiguities();

            *newCsm = csm;
            len = total_len;
            return ;
        }
        
    }
}

