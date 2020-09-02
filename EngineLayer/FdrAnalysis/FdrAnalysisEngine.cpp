#include "FdrAnalysisEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "FdrAnalysisResults.h"
#include "../GlobalVariables.h"

#include <boost/math/special_functions/gamma.hpp>

#include <numeric>
#include <math.h>
#include "Group.h"

namespace EngineLayer
{
    namespace FdrAnalysis
    {
        
        FdrAnalysisEngine::FdrAnalysisEngine(std::vector<PeptideSpectralMatch*> &psms,
                                             int massDiffAcceptorNumNotches,
                                             CommonParameters *commonParameters,
                                             std::vector<std::string> nestedIds,
                                             int verbosityLevel ) :
            MetaMorpheusEngine(commonParameters, nestedIds, verbosityLevel),
            MassDiffAcceptorNumNotches(massDiffAcceptorNumNotches),
            UseDeltaScore(commonParameters->getUseDeltaScore()),
            CalculateEValue(commonParameters->getCalculateEValue()),
            ScoreCutoff(commonParameters->getScoreCutoff())
        {
            AllPsms = psms;
        }
        
        MetaMorpheusEngineResults *FdrAnalysisEngine::RunSpecific()
        {
            FdrAnalysisResults *myAnalysisResults = new FdrAnalysisResults(this);
            
            Status("Running FDR analysis...");
            DoFalseDiscoveryRateAnalysis(myAnalysisResults);

            std::vector<PeptideSpectralMatch *> tmp;
            for ( auto b: AllPsms ) {
                if (b->getFdrInfo()->getQValue() < 0.01 ) {
                    tmp.push_back(b);
                }
            }
                                       
            myAnalysisResults->setPsmsWithin1PercentFdr(tmp.size() );
            
            return myAnalysisResults;
        }
        
        void FdrAnalysisEngine::DoFalseDiscoveryRateAnalysis(FdrAnalysisResults *myAnalysisResults)
        {
            // Stop if canceled
            if (GlobalVariables::getStopLoops())
            {
                return;
            }
            
            // calculate FDR on a per-protease basis (targets and decoys for a specific protease)
#ifdef ORIG
            auto psmsGroupedByProtease = AllPsms.GroupBy([&] (std::any p){
                    p::DigestionParams::Protease;
                });
#endif
            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f1 = [&]
                (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                return l->digestionParams->getProtease() < r->digestionParams->getProtease(); } ;
            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f2 = [&]
                (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                return l->digestionParams->getProtease() != r->digestionParams->getProtease(); } ;
            std::vector<std::vector<PeptideSpectralMatch *>> psmsGroupedByProtease = Group::GroupBy( AllPsms, f1, f2 );
            
            for (auto psms : psmsGroupedByProtease)
            {

                // generate the null distribution for e-value calculations
                double globalMeanScore = 0;
                int globalMeanCount = 0;
                
#ifdef ORIG
                //if (CalculateEValue && psms.Any())
#endif
                if (CalculateEValue && !psms.empty())
                {
                    std::vector<double> combinedScores;
                    for (auto psm : psms)
                    {
                        auto allscores = psm->getAllScores();
                        std::sort(allscores.begin(), allscores.end());
                        combinedScores.insert(combinedScores.end(), allscores.begin(), allscores.end());
                        
                        //remove top scoring peptide
#ifdef ORIG
                        //if (combinedScores.Any())
#endif
                        if ( !combinedScores.empty() ) 
                        {
                            combinedScores.pop_back();
                        }
                    }
                    
#ifdef ORIG
                    //if (combinedScores.Any())
#endif
                    if ( !combinedScores.empty() ) 
                    {
                        //globalMeanScore = combinedScores.Average();
                        globalMeanScore = std::accumulate(combinedScores.begin(),combinedScores.end(), 0.0)/combinedScores.size(); 
                        globalMeanCount = static_cast<int>(static_cast<double>(combinedScores.size()) / psms.size());
                    }
                    else
                    {
                        // should be a very rare case... if there are PSMs but each PSM only has one hit
                        globalMeanScore = 0;
                        globalMeanCount = 0;
                    }
                }
                
                //Calculate delta scores for the psms (regardless of if we are using them)
                for (auto psm : psms)
                {
                    if (psm != nullptr)
                    {
                        psm->CalculateDeltaScore(ScoreCutoff);
                    }
                }
                
                //determine if Score or DeltaScore performs better
                if (UseDeltaScore)
                {
                    constexpr double qValueCutoff = 0.01; //optimize to get the most PSMs at a 1% FDR
                    
#ifdef ORIG
                    std::vector<PeptideSpectralMatch*> scoreSorted = psms.OrderByDescending([&] (std::any b)   {
                            b::Score;
                        }).ThenBy([&] (std::any b)   {
                                b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass-b::PeptideMonisotopicMass->Value)
                                :  std::numeric_limits<double>::max();
                            }).GroupBy([&] (std::any b) {
                                    return std::tuple <std::string, int, std::optional<double>>(b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass)}
                                )->Select([&] (std::any b) {
                                        b::First();
                                    }).ToList();

#endif
                    //OrderByDescending and ThenBy
                    std::stable_sort (psms.begin(), psms.end(), [&] (PeptideSpectralMatch *r, PeptideSpectralMatch *l) {
                            if ( r->getScore() > l->getScore() ) return true;
                            if ( r->getScore() < l->getScore() ) return false;

                            double lval= std::numeric_limits<double>::max(), rval=std::numeric_limits<double>::max();
                            if ( l->getPeptideMonisotopicMass().has_value() ) {
                                lval = std::abs(l->getScanPrecursorMass() - l->getPeptideMonisotopicMass().value());
                            }
                            if ( r->getPeptideMonisotopicMass().has_value() ) {
                                rval = std::abs(r->getScanPrecursorMass() - r->getPeptideMonisotopicMass().value());
                            }
                            return rval < lval;
                        });
                    //GroupBy tuple
                    std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f3 = [&]
                        (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                        return l->getFullFilePath() != r->getFullFilePath() &&
                        l->getScanNumber() != r->getScanNumber() &&
                        l->getPeptideMonisotopicMass() != r->getPeptideMonisotopicMass(); } ;
                    std::vector<std::vector<PeptideSpectralMatch *>>tvec;
                    //Edgar: Note: cannot use the Group::GroupBy functionality here because of the custom sorting.
                    for ( auto p : psms ) {
                        bool found = false;
                        for ( auto q: tvec ) {
                            if ( f3 ( p, q[0] ) ) {
                                found = true;
                                q.push_back(p);
                                break;
                            }
                        }
                        if ( !found ) {
                            auto vec = new std::vector<PeptideSpectralMatch *>();
                            vec->push_back(p);
                            tvec.push_back(*vec);
                        }
                    }


                    // Select first
                    std::vector<PeptideSpectralMatch*> scoreSorted;
                    for ( auto t: tvec ) {
                        scoreSorted.push_back(t[0]);
                    }
            
                    int ScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);
#ifdef ORIG
                    scoreSorted = psms.OrderByDescending([&] (std::any b)	{
                            b::DeltaScore;
                        }).ThenBy([&] (std::any b){
                                b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) :
                                std::numeric_limits<double>::max();
                            }).GroupBy([&] (std::any b)		{
                                    return std::tuple < std::string>;
                                }, int, std::optional<double>>(b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass))->Select([&] (std::any b)  {
                                        b::First();
                                    }).ToList();
#endif
                    //OrderByDescending and ThenBy
                    std::stable_sort (psms.begin(), psms.end(), [&] (PeptideSpectralMatch *r, PeptideSpectralMatch *l) {
                            if ( r->getDeltaScore() > l->getDeltaScore() ) return true;
                            if ( r->getDeltaScore() < l->getDeltaScore() ) return false;

                            double lval= std::numeric_limits<double>::max(), rval=std::numeric_limits<double>::max();
                            if ( l->getPeptideMonisotopicMass().has_value() ) {
                                lval = std::abs(l->getScanPrecursorMass() - l->getPeptideMonisotopicMass().value());
                            }
                            if ( r->getPeptideMonisotopicMass().has_value() ) {
                                rval = std::abs(r->getScanPrecursorMass() - r->getPeptideMonisotopicMass().value());
                            }
                            return rval < lval;
                        });
                    //GroupBy tuple
                    tvec.clear();
                    for ( auto p : psms ) {
                        bool found = false;
                        for ( auto q: tvec ) {
                            if ( f3 ( p, q[0] ) ) {
                                found = true;
                                q.push_back(p);
                                break;
                            }
                        }
                        if ( !found ) {
                            auto vec = new std::vector<PeptideSpectralMatch *>();
                            vec->push_back(p);
                            tvec.push_back(*vec);
                        }
                    }
                    // Select first
                    scoreSorted.clear();
                    for ( auto t: tvec ) {
                        scoreSorted.push_back(t[0]);
                    }
                    int DeltaScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);
                    
                    //sort by best method
                    myAnalysisResults->setDeltaScoreImprovement(DeltaScorePSMs > ScorePSMs);
#ifdef ORIG
                    psms = myAnalysisResults->getDeltaScoreImprovement() ? psms.OrderByDescending([&] (std::any b)     {
                            b::DeltaScore;
                        }).ThenBy([&] (std::any b)  {
                                b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
                            }).ToList() :
                        psms.OrderByDescending([&] (std::any b) {
                                b::Score;
                            }).ThenBy([&] (std::any b){
                                    b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
                                }).ToList();
#endif
                    if ( myAnalysisResults->getDeltaScoreImprovement() ) {
                        //OrderByDescending and ThenBy
                        std::stable_sort (psms.begin(), psms.end(), [&] (PeptideSpectralMatch *r, PeptideSpectralMatch *l) {
                                if ( r->getDeltaScore() > l->getDeltaScore() ) return true;
                                if ( r->getDeltaScore() < l->getDeltaScore() ) return false;
                                
                                double lval= std::numeric_limits<double>::max(), rval=std::numeric_limits<double>::max();
                                if ( l->getPeptideMonisotopicMass().has_value() ) {
                                    lval = std::abs(l->getScanPrecursorMass() - l->getPeptideMonisotopicMass().value());
                                }
                                if ( r->getPeptideMonisotopicMass().has_value() ) {
                                    rval = std::abs(r->getScanPrecursorMass() - r->getPeptideMonisotopicMass().value());
                                }
                                return rval < lval;
                            });
                    }
                    else {
                    //OrderByDescending and ThenBy
                    std::stable_sort (psms.begin(), psms.end(), [&] (PeptideSpectralMatch *r, PeptideSpectralMatch *l) {
                            if ( r->getScore() > l->getScore() ) return true;
                            if ( r->getScore() < l->getScore() ) return false;

                            double lval= std::numeric_limits<double>::max(), rval=std::numeric_limits<double>::max();
                            if ( l->getPeptideMonisotopicMass().has_value() ) {
                                lval = std::abs(l->getScanPrecursorMass() - l->getPeptideMonisotopicMass().value());
                            }
                            if ( r->getPeptideMonisotopicMass().has_value() ) {
                                rval = std::abs(r->getScanPrecursorMass() - r->getPeptideMonisotopicMass().value());
                            }
                            return rval < lval;
                        });
                    }
                    
                }
                else //sort by score
                {
#ifdef ORIG
                    psms = psms.OrderByDescending([&] (std::any b) {
                            b::Score;
                        }).ThenBy([&] (std::any b)   {
                                b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) :
                                std::numeric_limits<double>::max();
                            }).ToList();
#endif
                    //OrderByDescending and ThenBy
                    std::stable_sort (psms.begin(), psms.end(), [&] (PeptideSpectralMatch *r, PeptideSpectralMatch *l) {
                            if ( r->getScore() > l->getScore() ) return true;
                            if ( r->getScore() < l->getScore() ) return false;

                            double lval= std::numeric_limits<double>::max(), rval=std::numeric_limits<double>::max();
                            if ( l->getPeptideMonisotopicMass().has_value() ) {
                                lval = std::abs(l->getScanPrecursorMass() - l->getPeptideMonisotopicMass().value());
                            }
                            if ( r->getPeptideMonisotopicMass().has_value() ) {
                                rval = std::abs(r->getScanPrecursorMass() - r->getPeptideMonisotopicMass().value());
                            }
                            return rval < lval;
                        });
                }
                
                double cumulativeTarget = 0;
                double cumulativeDecoy = 0;
                
                //set up arrays for local FDRs
                std::vector<double> cumulativeTargetPerNotch(MassDiffAcceptorNumNotches + 1);
                std::vector<double> cumulativeDecoyPerNotch(MassDiffAcceptorNumNotches + 1);
                
                //Assign FDR values to PSMs
                for (int i = 0; i < (int)psms.size(); i++)
                {
                    // Stop if canceled
                    if (GlobalVariables::getStopLoops())
                    {
                        break;
                    }
                    
                    PeptideSpectralMatch *psm = psms[i];
                    std::optional<int> tempVar = psm->getNotch();
                    int notch = tempVar.has_value() ? tempVar.value() : MassDiffAcceptorNumNotches;
                    if (psm->getIsDecoy())
                    {
                        // the PSM can be ambiguous between a target and a decoy sequence
                        // in that case, count it as the fraction of decoy hits
                        // e.g. if the PSM matched to 1 target and 2 decoys, it counts as 2/3 decoy
                        double decoyHits = 0;
                        double totalHits = 0;
#ifdef ORIG
                        auto hits = psm->BestMatchingPeptides.GroupBy([&] (std::any p)       {
                                p::Peptide::FullSequence;
                            });
#endif
                        std::vector<std::vector<std::tuple<int, PeptideWithSetModifications *>>>hits;
                        int current =0;
                        std::string currFullSequence, prevFullSequence;
                        auto tmpsms = psm->getBestMatchingPeptides();
                        for ( auto p = tmpsms.begin(); p != tmpsms.end(); p++ ) {
                            if ( p == tmpsms.begin() ) {
                                std::vector<std::tuple<int, PeptideWithSetModifications *>> *v = new std::vector<std::tuple<int, PeptideWithSetModifications *>>;
                                hits.push_back(*v);
                                prevFullSequence = std::get<1>(*p)->getFullSequence();
                            }
                            else {
                                auto q = p - 1;
                                currFullSequence = std::get<1>(*p)->getFullSequence();
                                if ( currFullSequence != prevFullSequence ) {
                                    std::vector<std::tuple<int, PeptideWithSetModifications *>> *v = new std::vector<std::tuple<int, PeptideWithSetModifications *>>;
                                    hits.push_back(*v);
                                    current++;
                                    prevFullSequence = currFullSequence;
                                }
                            }
                            hits[current].push_back(*p);
                        }
                        for (auto hit : hits)
                        {
#ifdef ORIG
                            if (hit->First().Peptide.Protein.IsDecoy)
#endif
                                if (std::get<1>(hit.front())->getProtein()->getIsDecoy())
                            {
                                decoyHits++;
                            }
                            totalHits++;
                        }
                        
                        cumulativeDecoy += decoyHits / totalHits;
                        cumulativeDecoyPerNotch[notch] += decoyHits / totalHits;
                    }
                    else
                    {
                        cumulativeTarget++;
                        cumulativeTargetPerNotch[notch]++;
                    }
                    
                    double qValue = std::min(1.0, cumulativeDecoy / cumulativeTarget);
                    double qValueNotch = std::min(1.0, cumulativeDecoyPerNotch[notch] / cumulativeTargetPerNotch[notch]);
                    
                    double maximumLikelihood = 0;
                    double eValue = 0;
                    double eScore = 0;
                    if (CalculateEValue)
                    {
                        eValue = GetEValue(psm, globalMeanCount, globalMeanScore, maximumLikelihood);
                        //eScore = -Math::Log(eValue, 10);
                        eScore = -log10(eValue);
                    }
                    
                    psm->SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetPerNotch[notch],
                                      cumulativeDecoyPerNotch[notch], qValueNotch, maximumLikelihood, eValue, eScore, CalculateEValue);
                }
                
                // set q-value thresholds such that a lower scoring PSM can't have 
                // a higher confidence than a higher scoring PSM
                //Populate min qValues
                double qValueThreshold = 1.0;
                std::vector<double> qValueNotchThreshold(MassDiffAcceptorNumNotches + 1);
                for (int i = 0; i < (int)qValueNotchThreshold.size(); i++)
                {
                    qValueNotchThreshold[i] = 1.0;
                }
                
                for (int i = (int)psms.size() - 1; i >= 0; i--)
                {
                    PeptideSpectralMatch *psm = psms[i];
                    
                    // threshold q-values
                    if (psm->getFdrInfo()->getQValue() > qValueThreshold)
                    {
                        psm->getFdrInfo()->setQValue(qValueThreshold);
                    }
                    else if (psm->getFdrInfo()->getQValue() < qValueThreshold)
                    {
                        qValueThreshold = psm->getFdrInfo()->getQValue();
                    }
                    
                    // threshold notch q-values
                    std::optional<int> tempVar2 = psm->getNotch();
                    int notch = tempVar2.has_value() ? tempVar2.value() : MassDiffAcceptorNumNotches;
                    if (psm->getFdrInfo()->getQValueNotch() > qValueNotchThreshold[notch])
                    {
                        psm->getFdrInfo()->setQValueNotch(qValueNotchThreshold[notch]);
                    }
                    else if (psm->getFdrInfo()->getQValueNotch() < qValueNotchThreshold[notch])
                    {
                        qValueNotchThreshold[notch] = psm->getFdrInfo()->getQValueNotch();
                    }
                }
            }
        }
        
        double FdrAnalysisEngine::GetEValue(PeptideSpectralMatch *psm, int globalMeanCount, double globalMeanScore, double &maximumLikelihood)
        {
            // get all of the PSM's scores for all hits, sort them, then remove the last value (the best score)
            std::vector<double> scoresWithoutBestHit;
            scoresWithoutBestHit.insert(scoresWithoutBestHit.end(), psm->getAllScores().begin(), psm->getAllScores().end());
            std::sort(scoresWithoutBestHit.begin(), scoresWithoutBestHit.end());
            
#ifdef ORIG
            if (scoresWithoutBestHit.Any())
#endif
            if (!scoresWithoutBestHit.empty())

            {
                scoresWithoutBestHit.pop_back();
            }
            
            // this is the "default" case for when there are no scores except the best hit
            // it uses a global mean score (all scores across all PSMs) to generate the null Poisson distribution
            // this will be overriden by the next few lines if there are enough scores in this PSM to estimate a null distribution
            //double preValue = SpecialFunctions::GammaLowerRegularized(globalMeanScore, psm->getScore());
            //Edgar Note: the boost gamma lower regularized function gamma_p does not accept the first argument to be 0.
            //            I verified with the C# version that the result returned by SpecialFunctions::GammaLowerRegularized
            //            is in that case 1, independent of the second argument. Mimicking the same here.
            double preValue = 1;
            if ( globalMeanScore > 0 ) {
                preValue = boost::math::gamma_p (globalMeanScore, psm->getScore());
            }
            maximumLikelihood = globalMeanScore;
            
            // calculate single-spectrum evalue if there are enough hits besides the best scoring peptide
            if (psm->getScore() == 0)
            {
                preValue = 1;
                maximumLikelihood = 0;
            }
#ifdef ORIG
            else if (scoresWithoutBestHit.Any())
#endif
            else if (!scoresWithoutBestHit.empty())
            {
                //maximumLikelihood = scoresWithoutBestHit.Average();
                maximumLikelihood = std::accumulate(scoresWithoutBestHit.begin(), scoresWithoutBestHit.end(), 0.0)/scoresWithoutBestHit.size();
                
                // this is the cumulative distribution for the poisson at each score up to but not including the score of the winner.
                // This is the probability that the winner has of getting that score at random by matching against a SINGLE spectrum
                if (maximumLikelihood > 0)
                {
                    //preValue = SpecialFunctions::GammaLowerRegularized(maximumLikelihood, psm->getScore());
                    //preValue = Math::GammaLowerRegularized(maximumLikelihood, psm->getScore());
                    preValue = boost::math::gamma_p ( maximumLikelihood, psm->getScore());
                }
            }
            
            // Now the probability of getting the winner's score goes up for each spectrum match.
            // We multiply the preValue by the number of theoretical spectrum within the tolerance to get this new probability.
            int count = scoresWithoutBestHit.size();
            if (count == 0)
            {
                count = globalMeanCount;
            }
            
            double probabilityOfScore = 1 - std::pow(preValue, count);
            return count * probabilityOfScore;
        }
        
        int FdrAnalysisEngine::GetNumPSMsAtqValueCutoff(std::vector<PeptideSpectralMatch*> &psms, double qValueCutoff)
        {
            int cumulative_target = 0;
            int cumulative_decoy = 0;
            for (auto psm : psms)
            {
                if (psm->getIsDecoy())
                {
                    cumulative_decoy++;
                    if (static_cast<double>(cumulative_decoy) / cumulative_target >= qValueCutoff)
                    {
                        return cumulative_target;
                    }
                }
                else
                {
                    cumulative_target++;
                }
            }
            return cumulative_target;
        }
    }
}
