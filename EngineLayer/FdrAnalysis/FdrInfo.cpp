#include "FdrInfo.h"

#include <string.h>
#include "BinaryPack.h"

namespace EngineLayer
{
	namespace FdrAnalysis
	{

            double FdrInfo::getCumulativeTarget() const
            {
                return privateCumulativeTarget;
            }
            
            void FdrInfo::setCumulativeTarget(double value)
            {
                privateCumulativeTarget = value;
            }
            
            double FdrInfo::getCumulativeDecoy() const
            {
                return privateCumulativeDecoy;
            }
            
            void FdrInfo::setCumulativeDecoy(double value)
            {
                privateCumulativeDecoy = value;
            }
            
            double FdrInfo::getCumulativeTargetNotch() const
            {
                return privateCumulativeTargetNotch;
            }
            
            void FdrInfo::setCumulativeTargetNotch(double value)
            {
                privateCumulativeTargetNotch = value;
            }
            
            double FdrInfo::getCumulativeDecoyNotch() const
            {
                return privateCumulativeDecoyNotch;
            }
            
            void FdrInfo::setCumulativeDecoyNotch(double value)
            {
                privateCumulativeDecoyNotch = value;
            }
            
            double FdrInfo::getQValue() const
            {
                return privateQValue;
            }
            
            void FdrInfo::setQValue(double value)
            {
                privateQValue = value;
            }
            
            double FdrInfo::getQValueNotch() const
            {
                return privateQValueNotch;
            }
            
            void FdrInfo::setQValueNotch(double value)
            {
                privateQValueNotch = value;
            }
            
            bool FdrInfo::getCalculateEValue() const
            {
                return privateCalculateEValue;
            }
            
            void FdrInfo::setCalculateEValue(bool value)
            {
                privateCalculateEValue = value;
            }
            
            double FdrInfo::getMaximumLikelihood() const
            {
                return privateMaximumLikelihood;
            }
            
            void FdrInfo::setMaximumLikelihood(double value)
            {
                privateMaximumLikelihood = value;
            }
            
            double FdrInfo::getEValue() const
            {
                return privateEValue;
            }
            
            void FdrInfo::setEValue(double value)
            {
                privateEValue = value;
            }
            
            double FdrInfo::getEScore() const
            {
                return privateEScore;
            }
            
            void FdrInfo::setEScore(double value)
            {
                privateEScore = value;
            }
            
            int FdrInfo::Pack( char* buf, size_t &buf_len, FdrInfo *fdr)
            {
                char tmpbuf[128];
                size_t pos = BinaryPack::LineStartOffset;
                int retlen;
                memset(tmpbuf, 0, 128);
                
                retlen = BinaryPack::PackBool(tmpbuf+pos, (fdr != nullptr) );
                pos += retlen;
                
                if ( fdr != nullptr ) {
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
                
                buf_len = pos;

                if ( pos > buf_len ) {
                    return -1;
                }
                memcpy ( buf, tmpbuf, pos );
                return (int)pos;
            }

            void FdrInfo::Unpack( char* buf, size_t &len, FdrInfo **newfdr )
            {
                double cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetNotch;
                double cumulativeDecoyNotch, qValueNotch, maximumLikelihood, eValue, eScore;
                bool calculateEValue=false;
                bool has_fdr = false;
                
                size_t pos = 0;
                int linelen, retlen;
                retlen = BinaryPack::GetLineLength(buf, linelen);
                pos += retlen;
                
                retlen = BinaryPack::UnpackBool ( buf+pos, has_fdr );
                pos += retlen;
                
                if ( has_fdr ) {
                    retlen = BinaryPack::UnpackDouble ( buf+pos, cumulativeTarget );
                    pos += retlen;
                    retlen = BinaryPack::UnpackDouble ( buf+pos, cumulativeDecoy );
                    pos += retlen;
                    retlen = BinaryPack::UnpackDouble ( buf+pos, qValue );
                    pos += retlen;
                    retlen = BinaryPack::UnpackDouble ( buf+pos, cumulativeTargetNotch );
                    pos += retlen;
                    retlen = BinaryPack::UnpackDouble ( buf+pos, cumulativeDecoyNotch );
                    pos += retlen;
                    retlen = BinaryPack::UnpackDouble ( buf+pos, qValueNotch );
                    pos += retlen;
                    retlen = BinaryPack::UnpackDouble ( buf+pos, maximumLikelihood );
                    pos += retlen;
                    retlen = BinaryPack::UnpackDouble ( buf+pos, eValue );
                    pos += retlen;
                    retlen = BinaryPack::UnpackDouble ( buf+pos, eScore );
                    pos += retlen;
                    
                    retlen = BinaryPack::UnpackBool ( buf+pos, calculateEValue);
                    pos += retlen;
                }
                
                if ( has_fdr ) {
                    FdrInfo* tempVar= new FdrInfo();
                    tempVar->setCumulativeTarget(cumulativeTarget);
                    tempVar->setCumulativeDecoy(cumulativeDecoy);
                    tempVar->setQValue(qValue);
                    tempVar->setCumulativeTargetNotch(cumulativeTargetNotch);
                    tempVar->setCumulativeDecoyNotch(cumulativeDecoyNotch);
                    tempVar->setQValueNotch(qValueNotch);
                    tempVar->setMaximumLikelihood(maximumLikelihood);
                    tempVar->setEScore(eScore);
                    tempVar->setEValue(eValue);
                    tempVar->setCalculateEValue(calculateEValue);

                    *newfdr = tempVar;
                }
                else {
                    *newfdr = nullptr;
                }

                len = linelen;
            }
        }
}
