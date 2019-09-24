#include "FdrAnalysisResults.h"
#include "FdrAnalysisEngine.h"

namespace EngineLayer
{
    namespace FdrAnalysis
    {
        
        FdrAnalysisResults::FdrAnalysisResults(FdrAnalysisEngine *s) : MetaMorpheusEngineResults(s)
        {
            setDeltaScoreImprovement(false);
        }
        
        int FdrAnalysisResults::getPsmsWithin1PercentFdr() const
        {
            return privatePsmsWithin1PercentFdr;
        }
        
        void FdrAnalysisResults::setPsmsWithin1PercentFdr(int value)
        {
            privatePsmsWithin1PercentFdr = value;
        }
        
        bool FdrAnalysisResults::getDeltaScoreImprovement() const
        {
            return privateDeltaScoreImprovement;
        }
        
        void FdrAnalysisResults::setDeltaScoreImprovement(bool value)
        {
            privateDeltaScoreImprovement = value;
        }
        
        std::string FdrAnalysisResults::ToString()
        {
            auto sb = new StringBuilder();
            sb->appendLine(MetaMorpheusEngineResults::ToString());
            sb->appendLine("PSMs within 1% fdr: " + std::to_string(getPsmsWithin1PercentFdr()));
            sb->appendLine("Delta Score Used for FDR Analysis: " + StringHelper::toString(getDeltaScoreImprovement()));
            
            std::string s = sb->toString();
            delete sb;
            return s;
        }
    }
}
