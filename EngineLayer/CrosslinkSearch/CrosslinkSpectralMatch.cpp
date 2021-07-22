#include "CrosslinkSpectralMatch.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "Crosslinker.h"

#include "Sort.h"
#include "BinaryPack.h"

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

            auto mFrIons = csm->getMatchedFragmentIons ();
            auto dp = csm->digestionParams;
            auto uMapPep = csm->getPeptidesToMatchingFragments();
            std::vector<int> lPositions  = csm->getLinkPositions();
            std::vector<int> xlRanks = csm->getXlRank();
            bool has_beta_peptide = csm->getBetaPeptide() != nullptr;          

            size_t pos = BinaryPack::LineStartOffset;
            char tmpbuf[256];
            int retlen;
            
            // line 1
            retlen = BinaryPack::PackBool(tmpbuf+pos, csm->getNotch().has_value());
            pos += retlen;
            if ( csm->getNotch().has_value() ) {
                retlen = BinaryPack::PackInt(tmpbuf+pos, csm->getNotch().value());
                pos += retlen;
            }
            
            retlen = BinaryPack::PackDouble(tmpbuf+pos, csm->getXLTotalScore() );
            pos += retlen;
            retlen = BinaryPack::PackDouble(tmpbuf+pos, csm->getDeltaScore() );
            pos += retlen;
            retlen = BinaryPack::PackDouble(tmpbuf+pos, csm->getScore() );
            pos += retlen;
            retlen = BinaryPack::PackDouble(tmpbuf+pos, csm->getRunnerUpScore() );
            pos += retlen;
            retlen = BinaryPack::PackDouble(tmpbuf+pos, csm->getPeptideMonisotopicMass().value() );
            pos += retlen;

            retlen = BinaryPack::PackInt(tmpbuf+pos, csm->getScanIndex() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, csm->getScanNumber() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, csm->getXlProteinPos() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, (int)csm->getMatchedFragmentIons().size() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, (int)lPositions.size() );
            pos += retlen;
            retlen = BinaryPack::PackInt(tmpbuf+pos, (int)xlRanks.size() );
            pos += retlen;
                        
            retlen = BinaryPack::PackBool(tmpbuf+pos, has_beta_peptide );
            pos += retlen;

            PsmCrossType ctype = csm->getCrossType();
            retlen = BinaryPack::PackString(tmpbuf+pos, PsmCrossTypeToString(ctype) );
            pos += retlen;

            
            // set lenght of line right at the beginning.
            BinaryPack::SetLineLength(tmpbuf, pos);

            if ( bufpos+pos > buf_len ) {
                buf_len = bufpos+pos;
                return -1;
            }
            memcpy ( buf+bufpos, tmpbuf, pos );
            bufpos += pos;

            //Line 2: pack FdrInfo 
            FdrInfo *fdr = csm->getFdrInfo();
            size_t tmp_len = buf_len - bufpos;

            // this routine sets all the required aspects of a packed line (e.g. header, length)
            int ret = FdrInfo::Pack(buf+bufpos, tmp_len, fdr );
            if ( ret == -1 ) {
                buf_len += tmp_len - (buf_len - bufpos);
                return -1;
            }
            bufpos += tmp_len;

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
            tmp_len = buf_len - bufpos;

            // this routine sets all the required aspects of a packed line (e.g. header, length)
            ret = PeptideWithSetModifications::Pack(buf+bufpos, tmp_len, pep );
            if ( ret == -1 ) {
                buf_len += tmp_len - (buf_len - bufpos);
                return -1;
            }
            bufpos += tmp_len;
            
            //line 11-x: one line for each MatchedFragmentIon
            for ( auto i=0; i< mFrIons.size(); i++ ) {
                tmp_len = buf_len - bufpos;

                // dito here, header and length should be set
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
                                             const std::vector<Modification*> &mods,                                             
                                             const std::vector<Ms2ScanWithSpecificMass*> &ms2Scans,
                                             const std::vector<Protein *> &proteinList )
        {
            std::vector<char *> lines = BinaryPack::SplitLines(buf, buf_len);

            size_t total_len=0;
            int counter=0;
            for (auto  i=0; i < lines.size();  ) {
                size_t tmp_len=0;
                CrosslinkSpectralMatch *pep;
                bool has_beta_peptide=false;
                CrosslinkSpectralMatch::Unpack_internal(lines, i, tmp_len, &pep, mods, ms2Scans, proteinList,
                                                        has_beta_peptide );
                total_len += tmp_len;
                pepVec.push_back(pep);
                if ( has_beta_peptide ) {
                    CrosslinkSpectralMatch *beta_pep;
                    CrosslinkSpectralMatch::Unpack_internal(lines, i, tmp_len, &beta_pep, mods, ms2Scans, proteinList,
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
                                             const std::vector<Modification*> &mods,
                                             const std::vector<Ms2ScanWithSpecificMass*> &ms2Scans,
                                             const std::vector<Protein *> &proteinList )
        {
            std::vector<char *> lines = BinaryPack::SplitLines(buf, buf_len);
            int index=0;
            if ( lines.size() < 10 ) {
                std::cout << "CrosslinkSpectralMatch::Unpack : input does not contain enough information to " <<
                    "reconstruct the CrosslinkSpectralMatch. " << std::endl;
                return;
            }
            bool has_beta_peptide=false;            
            CrosslinkSpectralMatch::Unpack_internal ( lines, index, len, newCsm, mods, ms2Scans, proteinList,
                                                      has_beta_peptide );
            if ( has_beta_peptide) {
                CrosslinkSpectralMatch* beta_pep;
                size_t tmp_len=0;
                CrosslinkSpectralMatch::Unpack_internal ( lines, index, tmp_len, &beta_pep, mods, ms2Scans, proteinList,
                                                          has_beta_peptide );
                (*newCsm)->setBetaPeptide(beta_pep);
                len += tmp_len;
            }            
        }

        void CrosslinkSpectralMatch::Unpack_internal (std::vector<char*> &input,
                                                      int &index, size_t &len,
                                                      CrosslinkSpectralMatch** newCsm,
                                                      const std::vector<Modification*> &mods,
                                                      const std::vector<Ms2ScanWithSpecificMass*> &ms2Scans,
                                                      const std::vector<Protein *> &proteinList,
                                                      bool &has_beta_peptide )
        {
            size_t total_len = 0;
            int linelen=0;
            int retlen, pos=0;
            char *buf=NULL;
            //Dissect line 1: generic information
            buf = input[index];
            retlen = BinaryPack::GetLineLength( buf, linelen);
            pos += retlen;
            total_len += linelen; 
            index++;
            
            int notch=-1, scanindex, scannumber, proteinPos, matchedFragmentIonsVecsize, lpositionsize, xlranksize;
            double  deltaScore, XLTotalScore, score, runnerUpScore, peptideMonisotopicMass;
            bool  tmpvar;
            
            retlen = BinaryPack::UnpackBool ( buf+pos, tmpvar );
            pos += retlen;
            if ( tmpvar ) {
                retlen = BinaryPack::UnpackInt(buf+pos, notch );
                pos += retlen;
            }
            retlen = BinaryPack::UnpackDouble ( buf+pos, XLTotalScore );
            pos += retlen;
            retlen = BinaryPack::UnpackDouble ( buf+pos, deltaScore );
            pos += retlen;
            retlen = BinaryPack::UnpackDouble ( buf+pos, score );
            pos += retlen;
            retlen = BinaryPack::UnpackDouble ( buf+pos, runnerUpScore );
            pos += retlen;
            retlen = BinaryPack::UnpackDouble ( buf+pos, peptideMonisotopicMass );
            pos += retlen;

            retlen = BinaryPack::UnpackInt ( buf+pos, scanindex );
            pos += retlen;
            retlen = BinaryPack::UnpackInt ( buf+pos, scannumber );
            pos += retlen;
            retlen = BinaryPack::UnpackInt ( buf+pos, proteinPos );
            pos += retlen;
            retlen = BinaryPack::UnpackInt ( buf+pos, matchedFragmentIonsVecsize );
            pos += retlen;
            retlen = BinaryPack::UnpackInt ( buf+pos, lpositionsize );
            pos += retlen;
            retlen = BinaryPack::UnpackInt ( buf+pos, xlranksize );
            pos += retlen;

            retlen = BinaryPack::UnpackBool ( buf+pos, has_beta_peptide );
            pos += retlen;

            std::string tmpstring;
            retlen = BinaryPack::UnpackString ( buf+pos, tmpstring );
            pos += retlen;
            PsmCrossType ctype = PsmCrossTypeFromString(tmpstring);
            
            //line 2: FdrInfo related data
            FdrInfo* fdr=nullptr;
            size_t tmp_len=0;
            FdrInfo::Unpack(input[index], tmp_len, &fdr);
            total_len += tmp_len;
            index ++;

            //line 3: linkPositions
            pos = 0;
            buf = input[index];
            index++;
            retlen = BinaryPack::GetLineLength(buf, linelen);
            pos += retlen;
            total_len += linelen;            

            std::vector<int> linkPosvec;
            for ( auto i=0; i<lpositionsize; i++ ) {
                int tmpint;
                retlen = BinaryPack::UnpackInt ( buf+pos, tmpint );
                pos += retlen;
                linkPosvec.push_back(tmpint);
            }
            
            //line 4: xlRank
            pos = 0;
            buf = input[index];
            index++;
            retlen = BinaryPack::GetLineLength(buf, linelen);
            pos += retlen;
            total_len += linelen;            

            std::vector<int> xlRankVec;
            for ( auto i=0; i<xlranksize; i++ ) {
                int tmpint;
                retlen = BinaryPack::UnpackInt ( buf+pos, tmpint );
                pos += retlen;
                xlRankVec.push_back(tmpint);
            }
            
            //line 5: DigestionParams
            pos = 0;
            buf = input[index];
            index++;
            retlen = BinaryPack::GetLineLength(buf, linelen);
            pos += retlen;
            total_len += linelen;            

            std::string dpstring;
            retlen = BinaryPack::UnpackString(buf+pos, dpstring);
            pos += retlen;
            DigestionParams *dp = DigestionParams::FromString(dpstring);

            //line 6-10: PeptideWithSetModifications
            PeptideWithSetModifications* pep;
            tmp_len=0;
            PeptideWithSetModifications::Unpack(input, index, tmp_len, &pep);
            pep->SetNonSerializedPeptideInfo ( mods, proteinList );
            total_len += tmp_len;
            index += 4;

#ifdef DEBUG
            // Safety check:
            double monIsotopicMass = pep->getMonoisotopicMass();
            if ( monIsotopicMass != peptideMonisotopicMass ) {
                std::cout << "Safety check failed when reconstructing PeptideWithSetMOdifications in CrosslinkSpectralMatch::Unpack(). " <<
                    "monIsotopicMass is " << monIsotopicMass << " should be " << peptideMonisotopicMass << std::endl;
            }
#endif
            
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
            if ( fdr != nullptr ) {
                csm->setFdrInfo(fdr);
            }
            csm->ResolveAllAmbiguities();

            // This is a hack. Need to find a proper solution.
            // In some situations involving Mods, the PeptideMonoisotopicMass
            // is not correct after deserialization.
            csm->setPeptideMonisotopicMass(std::make_optional(peptideMonisotopicMass));

            
            *newCsm = csm;
            len = total_len;
            return ;
        }
        
    }
}

