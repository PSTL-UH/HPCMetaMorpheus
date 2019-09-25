#include "BestPeptideScoreNotch.h"

using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
    namespace CrosslinkSearch
    {
        
        BestPeptideScoreNotch::BestPeptideScoreNotch(PeptideWithSetModifications *bestPeptide, double bestScore, int bestNotch)
        {
            setBestPeptide(bestPeptide);
            setBestScore(bestScore);
            setBestNotch(bestNotch);
        }
        
        PeptideWithSetModifications *BestPeptideScoreNotch::getBestPeptide() const
        {
            return privateBestPeptide;
        }
        
        void BestPeptideScoreNotch::setBestPeptide(PeptideWithSetModifications *value)
        {
            privateBestPeptide = value;
        }
        
        double BestPeptideScoreNotch::getBestScore() const
        {
            return privateBestScore;
        }
        
        void BestPeptideScoreNotch::setBestScore(double value)
        {
            privateBestScore = value;
        }
        
        int BestPeptideScoreNotch::getBestNotch() const
        {
            return privateBestNotch;
        }
        
        void BestPeptideScoreNotch::setBestNotch(int value)
        {
            privateBestNotch = value;
        }
        
        std::vector<int> BestPeptideScoreNotch::getTopPosition() const
        {
            return privateTopPosition;
        }
        
        void BestPeptideScoreNotch::setTopPosition(const std::vector<int> &value)
        {
            privateTopPosition = value;
        }
    }
}
