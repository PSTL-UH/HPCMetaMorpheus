#include "CrosslinkSearchEngine.h"
#include "CrosslinkSpectralMatch.h"
#include "Crosslinker.h"
#include "../PrecursorSearchModes/MassDiffAcceptor.h"
#include "../PrecursorSearchModes/OpenMassDiffAcceptor.h"
#include "../CommonParameters.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../PrecursorSearchModes/SingleAbsoluteAroundZeroMassDiffAcceptor.h"
#include "../MetaMorpheusEngineResults.h"
#include "../EventArgs/ProgressEventArgs.h"
#include "BestPeptideScoreNotch.h"
#include "PsmCrossType.h"
#include "CrosslinkedPeptides.h"

#include "../GlobalVariables.h"


#ifdef TIMING_INFO
#include <chrono>
#endif

#include "omp.h"

using namespace EngineLayer::ModernSearch;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
    namespace CrosslinkSearch
    {
        
        CrosslinkSearchEngine::CrosslinkSearchEngine(std::vector<CrosslinkSpectralMatch*> &globalCsms,
                                                     std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans,
                                                     std::vector<PeptideWithSetModifications*> &peptideIndex,
                                                     std::vector<std::vector<int>> &fragmentIndex,
                                                     int currentPartition,
                                                     CommonParameters *commonParameters,
                                                     Crosslinker *crosslinker,
                                                     bool CrosslinkSearchTop,
                                                     int CrosslinkSearchTopNum,
                                                     bool quench_H2O,
                                                     bool quench_NH2,
                                                     bool quench_Tris,
                                                     std::vector<std::string> nestedIds,
                                                     int verbosityLevel) :
            ModernSearchEngine(unusedPsms,
                               listOfSortedms2Scans,
                               peptideIndex,
                               fragmentIndex,
                               currentPartition,
                               commonParameters,
                               new OpenSearchMode(), 0.0,
                               nestedIds,
                               verbosityLevel ),
            GlobalCsms(globalCsms),
            CrosslinkSearchTopN(CrosslinkSearchTop),
            TopN(CrosslinkSearchTopNum),
            QuenchH2O(quench_H2O),
            QuenchNH2(quench_NH2),
            QuenchTris(quench_Tris)
        {
            privateCrosslinker = new Crosslinker(crosslinker);
            GenerateCrosslinkModifications(crosslinker);
            //AllCrosslinkerSites = privateCrosslinker->getCrosslinkerModSites().ToCharArray().Concat(
            //    privateCrosslinker->getCrosslinkerModSites2().ToCharArray())->Distinct()->ToArray();
            std::string s1 = privateCrosslinker->getCrosslinkerModSites();
            std::string s2 = privateCrosslinker->getCrosslinkerModSites2();
            std::vector<char> vs1 (s1.begin(), s1.end() );
            std::vector<char> vs2 (s2.begin(), s2.end() );
            
            AllCrosslinkerSites = vs1;
            for ( auto c : vs2 ) {
                bool found = false;
                for ( auto c2 : AllCrosslinkerSites ) {
                    if ( c == c2 ) {
                        found = true;
                        break;
                    }
                }
                if ( !found ) {
                    AllCrosslinkerSites.push_back(c);
                }
            }
            
            if (dynamic_cast<PpmTolerance*>(commonParameters->getPrecursorMassTolerance()) != nullptr)
            {
                XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(commonParameters->getPrecursorMassTolerance()->getValue());
            }
            else
            {
                XLPrecusorSearchMode = new SingleAbsoluteAroundZeroSearchMode(commonParameters->getPrecursorMassTolerance()->getValue());
            }
        }
        
        MetaMorpheusEngineResults *CrosslinkSearchEngine::RunSpecific()
        {
            int oldPercentProgress = 0;
#ifdef TIMING_INFO
            long long int indexedScoringTime=0, bestPeptideScoreNotchTime=0, findCrosslinkedPeptideTime=0, binsSearchTime=0;
            int iteration_count= (int)ListOfSortedMs2Scans.size();
            auto tx = std::chrono::high_resolution_clock::now();
#endif

            ReportEngineProgress("Performing crosslink search... " + std::to_string(CurrentPartition) + "/" +
                                 std::to_string(commonParameters->getTotalPartitions()), oldPercentProgress);
            
            unsigned char byteScoreCutoff = static_cast<unsigned char>(commonParameters->getScoreCutoff());
            
            int ListOfSortedMs2Scanssize = (int) ListOfSortedMs2Scans.size();
            int PeptideIndexsize = PeptideIndex.size();

#pragma omp parallel
         {
             double progress = 0;
             std::map <int, CrosslinkSpectralMatch *> localinfo; // per thread
             int tid = omp_get_thread_num();
             int num_threads = omp_get_num_threads();
             
#pragma omp for schedule(guided)
            for (int scanIndex = 0; scanIndex < ListOfSortedMs2Scanssize; scanIndex++)
            {            
                std::vector<unsigned char> scoringTable(PeptideIndexsize );
                std::vector<int> idsOfPeptidesPossiblyObserved;
                auto scan = ListOfSortedMs2Scans[scanIndex];

#ifdef TIMING_INFO
                auto t0 = std::chrono::high_resolution_clock::now();
#endif
                // get fragment bins for this scan
                std::vector<int> allBinsToSearch = GetBinsToSearch(scan);
                std::vector<BestPeptideScoreNotch*> bestPeptideScoreNotchList;
#ifdef TIMING_INFO
                auto t0e = std::chrono::high_resolution_clock::now();
                binsSearchTime += std::chrono::duration_cast<std::chrono::microseconds>( t0e - t0 ).count();
#endif
                
                // first-pass scoring
#ifdef TIMING_INFO
                auto t1 = std::chrono::high_resolution_clock::now();
#endif
                IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan->getPrecursorMass(),
                               -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                               PeptideIndex, massDiffAcceptor, 0);
#ifdef TIMING_INFO
                auto t1e = std::chrono::high_resolution_clock::now();
                indexedScoringTime += std::chrono::duration_cast<std::chrono::microseconds>( t1e - t1 ).count();
#endif

                // done with indexed scoring; refine scores and create PSMs
                if (!idsOfPeptidesPossiblyObserved.empty())
                {
                    if (CrosslinkSearchTopN)
                    {
                        // take top N hits for this scan
                        std::sort(idsOfPeptidesPossiblyObserved.begin(), idsOfPeptidesPossiblyObserved.end(), [&] (int l, int r) {
                                return scoringTable[l] > scoringTable[r];}
                            );
                        if (  (int)idsOfPeptidesPossiblyObserved.size() > TopN) {
                            idsOfPeptidesPossiblyObserved.resize(TopN);
                        }
                    }
                       
#ifdef TIMING_INFO
                    auto t2 = std::chrono::high_resolution_clock::now();
#endif
                    for (auto id : idsOfPeptidesPossiblyObserved)
                    {
                        PeptideWithSetModifications *peptide = PeptideIndex[id];
                        
                        int notch = massDiffAcceptor->Accepts(scan->getPrecursorMass(), peptide->getMonoisotopicMass());
                        auto tempVar3 = new BestPeptideScoreNotch (peptide, scoringTable[id], notch);
                        bestPeptideScoreNotchList.push_back(tempVar3);
                    }
                    
                    // combine individual peptide hits with crosslinker mass to find best crosslink PSM hit
#ifdef TIMING_INFO
                    auto t2e = std::chrono::high_resolution_clock::now();
                    bestPeptideScoreNotchTime += std::chrono::duration_cast<std::chrono::microseconds>( t2e - t2 ).count();

#endif
                    auto csm = FindCrosslinkedPeptide(scan, bestPeptideScoreNotchList, scanIndex);
#ifdef TIMING_INFO
                    auto t3e = std::chrono::high_resolution_clock::now();
                    findCrosslinkedPeptideTime += std::chrono::duration_cast<std::chrono::microseconds>( t3e - t2e ).count();
#endif
                    
                    if (csm == nullptr)
                    {
                        progress++;
                        continue;
                    }
                    
                    // this scan might already have a hit from a different database partition; check to see if the score improves
                    if (GlobalCsms[scanIndex] == nullptr )
                    {
                        localinfo[scanIndex] = csm;
                    }
                    else if ( localinfo[scanIndex]->getXLTotalScore() < csm->getXLTotalScore() )
                    {
                        delete localinfo[scanIndex];
                        localinfo[scanIndex] = csm;
                    }
                }
                
                if ( tid == 0 ) 
                {
                    // report search progress
                    progress++;
                    auto percentProgress = static_cast<int>(((progress * omp_get_num_threads()) / ListOfSortedMs2Scanssize) * 100);                
                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportEngineProgress("Performing crosslink search... " + std::to_string(CurrentPartition) + "/" +
                                             std::to_string(commonParameters->getTotalPartitions()), percentProgress);
                    }
                }
                for ( auto p:  bestPeptideScoreNotchList ) {
                    delete p;
                }                
            }
#pragma omp critical
            {
                for ( auto p = localinfo.begin(); p!= localinfo.end(); p++ ) {
                    auto key = p->first;
                    auto value = p->second;
                    GlobalCsms[key] = value;
                }
            }

         } // omp parallel

#ifdef TIMING_INFO
            auto txe = std::chrono::high_resolution_clock::now();
            std::cout << " CrosslinkedSearchEngine        : " << std::chrono::duration_cast<std::chrono::milliseconds>( txe - tx ).count() << " ms in " << iteration_count << " iterations " << std::endl;
            std::cout << "     binsSearchTime             : " << ((double)binsSearchTime)/1000 << std::endl;
            std::cout << "     indexedScoringTime         : " << ((double)indexedScoringTime)/1000 << std::endl;
            std::cout << "     bestPeptideScoreNotchTime  : " << ((double)bestPeptideScoreNotchTime)/1000 << std::endl;
            std::cout << "     findCrosslinkedPeptideTime : " << ((double)findCrosslinkedPeptideTime)/1000 << std::endl;
#endif
            
            return new MetaMorpheusEngineResults(this);
        }
        
        void CrosslinkSearchEngine::GenerateCrosslinkModifications(Crosslinker *crosslinker)
        {
            ModificationMotif *motif;
            ModificationMotif::TryGetMotif("X", &motif);
            std::string oId = "";
            std::string acc = "";
            std::string modType = "";
            std::string feaType = "";
            std::string locRestr = "Unassigned.";
            Chemistry::ChemicalFormula *chemForm = nullptr;
            TrisDeadEnd = new Modification( "Tris Dead End", acc, "Crosslink", feaType, motif, 
                                            "Anywhere.", chemForm, privateCrosslinker->getDeadendMassTris());
            H2ODeadEnd = new Modification( "H2O Dead End", acc, "Crosslink", feaType, motif,
                                           "Anywhere.", chemForm, privateCrosslinker->getDeadendMassH2O());
            NH2DeadEnd = new Modification( "NH2 Dead End", acc, "Crosslink", feaType, motif, 
                                           "Anywhere.", chemForm, privateCrosslinker->getDeadendMassNH2());
            Loop = new Modification( "Loop", acc, "Crosslink", feaType, motif, 
                                     "Anywhere.", chemForm, privateCrosslinker->getLoopMass());
            
        }
        
        CrosslinkSpectralMatch *CrosslinkSearchEngine::FindCrosslinkedPeptide(Ms2ScanWithSpecificMass *theScan,
                                                                              std::vector<BestPeptideScoreNotch*> &theScanBestPeptide,
                                                                              int scanIndex)
        {
            std::vector<CrosslinkSpectralMatch*> possibleMatches;
            for (int alphaIndex = 0; alphaIndex < (int)theScanBestPeptide.size(); alphaIndex++)
            {
                PeptideWithSetModifications *bestPeptide = theScanBestPeptide[alphaIndex]->getBestPeptide();
                
                //Single Peptide
                if (XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), bestPeptide->getMonoisotopicMass()) >= 0)
                {
                    std::vector<Product*> products = bestPeptide->Fragment(commonParameters->getDissociationType(),
                                                                           FragmentationTerminus::Both);//.ToList();
                    auto matchedFragmentIons = MatchFragmentIons(theScan, products, commonParameters);
                    double score = CalculatePeptideScore(theScan->getTheScan(), matchedFragmentIons, 0);
                    
                    auto psmCrossSingle = new CrosslinkSpectralMatch(bestPeptide, theScanBestPeptide[alphaIndex]->getBestNotch(),
                                                                     score, scanIndex, theScan, commonParameters->getDigestionParams(),
                                                                     matchedFragmentIons);
                    psmCrossSingle->setCrossType(PsmCrossType::Single);
                    psmCrossSingle->setXlRank(std::vector<int> {alphaIndex});
                    
                    possibleMatches.push_back(psmCrossSingle);                    
                }
                // Deadend Peptide
                else if (QuenchTris && XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), bestPeptide->getMonoisotopicMass() +
                                                                     privateCrosslinker->getDeadendMassTris()) >= 0)
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    
                    if (!possibleCrosslinkLocations.empty())
                    {
                        // tris deadend
                        auto csm = LocalizeDeadEndSite(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, TrisDeadEnd,
                                                       theScanBestPeptide[alphaIndex]->getBestNotch(), scanIndex, alphaIndex);
                        possibleMatches.push_back(csm);
                    }
                }
                else if (QuenchH2O && XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), bestPeptide->getMonoisotopicMass() +
                                                                    privateCrosslinker->getDeadendMassH2O()) >= 0)
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    
                    if (!possibleCrosslinkLocations.empty())
                    {
                        // H2O deadend
                        auto csm = LocalizeDeadEndSite(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, H2ODeadEnd,
                                                       theScanBestPeptide[alphaIndex]->getBestNotch(), scanIndex, alphaIndex);
                        possibleMatches.push_back(csm);
                    }
                }
                else if (QuenchNH2 && XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), bestPeptide->getMonoisotopicMass() +
                                                                    privateCrosslinker->getDeadendMassNH2()) >= 0)
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    
                    if (!possibleCrosslinkLocations.empty())
                    {
                        // NH2 deadend
                        auto csm = LocalizeDeadEndSite(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, NH2DeadEnd,
                                                       theScanBestPeptide[alphaIndex]->getBestNotch(), scanIndex, alphaIndex);
                        possibleMatches.push_back(csm);
                    }
                }
                // loop peptide
                else if (privateCrosslinker->getLoopMass() != 0 && XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(),
                                                                                                 bestPeptide->getMonoisotopicMass() + privateCrosslinker->getLoopMass()) >= 0)
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    
                    if (possibleCrosslinkLocations.size() >= 2)
                    {
                        possibleMatches.push_back(LocalizeLoopSites(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, Loop,
                                                                    theScanBestPeptide[alphaIndex]->getBestNotch(), scanIndex, alphaIndex));
                    }
                }
                // Cross-linked peptide
                else if (theScan->getPrecursorMass() - bestPeptide->getMonoisotopicMass() >= (commonParameters->getDigestionParams()->getMinPeptideLength() * 50))
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    if (possibleCrosslinkLocations.empty())
                    {
                        continue;
                    }
                    
                    PeptideWithSetModifications *alphaPeptide = bestPeptide;
                    
                    for (int betaIndex = 0; betaIndex < (int)theScanBestPeptide.size(); betaIndex++)
                    {
                        PeptideWithSetModifications *betaPeptide = theScanBestPeptide[betaIndex]->getBestPeptide();
                        
                        if (XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), alphaPeptide->getMonoisotopicMass() +
                                                          betaPeptide->getMonoisotopicMass() + privateCrosslinker->getTotalMass()) >= 0)
                        {
                            std::vector<int> possibleBetaCrosslinkSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, betaPeptide);
                            
                            if (possibleBetaCrosslinkSites.empty())
                            {
                                continue;
                            }
                            
                            CrosslinkSpectralMatch *csm = LocalizeCrosslinkSites(theScan, theScanBestPeptide[alphaIndex],
                                                                                 theScanBestPeptide[betaIndex], privateCrosslinker,
                                                                                 alphaIndex, betaIndex);
                            
                            possibleMatches.push_back(csm);
                        }
                    }
                }
            }
            
            // get the best match for this spectrum
            // bestPsmCross will be null if there are no valid hits
            for ( auto v = possibleMatches.begin(); v != possibleMatches.end(); ) {
                if ( *v == nullptr ) {
                    possibleMatches.erase(v);
                }
                else {
                    v++;
                }
            }
            
            std::sort(possibleMatches.begin(), possibleMatches.end(), [&] (CrosslinkSpectralMatch *l, CrosslinkSpectralMatch *r) {
                    return l->getXLTotalScore() > r->getXLTotalScore(); }
                );
            //auto bestPsmCross = possibleMatches.FirstOrDefault();
            CrosslinkSpectralMatch* bestPsmCross = nullptr;
            if ( possibleMatches.size() > 0 ) {
                bestPsmCross = possibleMatches.front();
            }
            
            // resolve ambiguities
            if (bestPsmCross != nullptr)
            {
                bestPsmCross->ResolveAllAmbiguities();
                
                if (bestPsmCross->getBetaPeptide() != nullptr)
                {
                    bestPsmCross->getBetaPeptide()->ResolveAllAmbiguities();
                }
            }
            
            // calculate delta score
            if (possibleMatches.size() > 1)
            {
                bestPsmCross->setDeltaScore( possibleMatches[0]->getXLTotalScore() - possibleMatches[1]->getXLTotalScore());
            }

            //Some memory management to reduce the leaks
            for ( auto pM : possibleMatches ) {
                if ( pM != nullptr && pM != bestPsmCross ) {
                    delete pM;
                }
            }
            
            return bestPsmCross;
        }
        
        CrosslinkSpectralMatch *CrosslinkSearchEngine::LocalizeCrosslinkSites(Ms2ScanWithSpecificMass *theScan, BestPeptideScoreNotch *alphaPeptide,
                                                                              BestPeptideScoreNotch *betaPeptide, Crosslinker *crosslinker,
                                                                              int ind, int inx)
        {
            CrosslinkSpectralMatch *localizedCrosslinkedSpectralMatch = nullptr;
            
            std::vector<std::tuple<std::vector<int>, std::vector<int>>> pairs;
            std::string s = crosslinker->getCrosslinkerModSites();
            std::vector<char> vs (s.begin(), s.end() );
            std::string s2 = crosslinker->getCrosslinkerModSites2();
            std::vector<char> vs2 (s2.begin(), s2.end() );
            
            if (crosslinker->getCrosslinkerModSites() == crosslinker->getCrosslinkerModSites2())
            {
                std::vector<int> possibleAlphaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(
                    vs, alphaPeptide->getBestPeptide());
                std::vector<int> possibleBetaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(
                    vs, betaPeptide->getBestPeptide());
                
                pairs.push_back(std::tuple<std::vector<int>, std::vector<int>>(possibleAlphaXlSites, possibleBetaXlSites));
            }
            else
            {
                std::vector<int> possibleAlphaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(
                    vs, alphaPeptide->getBestPeptide());
                std::vector<int> possibleBetaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(
                    vs2, betaPeptide->getBestPeptide());
                
                pairs.push_back(std::tuple<std::vector<int>, std::vector<int>>(possibleAlphaXlSites, possibleBetaXlSites));
                
                possibleAlphaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(
                    vs2, alphaPeptide->getBestPeptide());
                possibleBetaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(
                    vs, betaPeptide->getBestPeptide());
                
                pairs.push_back(std::tuple<std::vector<int>, std::vector<int>>(possibleAlphaXlSites, possibleBetaXlSites));
            }
            
            for (auto pair : pairs)
            {
                std::vector<int> possibleAlphaXlSites = std::get<0>(pair);
                std::vector<int> possibleBetaXlSites = std::get<1>(pair);
                
                if (!possibleAlphaXlSites.empty() && !possibleBetaXlSites.empty())
                {
                    int bestAlphaSite = 0;
                    int bestBetaSite = 0;
                    std::vector<MatchedFragmentIon*> bestMatchedAlphaIons;
                    std::vector<MatchedFragmentIon*> bestMatchedBetaIons;
                    double bestAlphaLocalizedScore = 0;
                    double bestBetaLocalizedScore = 0;
                    
                    auto fragmentsForEachAlphaLocalizedPossibility = CrosslinkedPeptide::XlGetTheoreticalFragments(
                        commonParameters->getDissociationType(), privateCrosslinker, possibleAlphaXlSites,
                        betaPeptide->getBestPeptide()->getMonoisotopicMass(), alphaPeptide->getBestPeptide());//.ToList();
                    
                    for (auto possibleSite : possibleAlphaXlSites)
                    {
                        for (auto setOfFragments : fragmentsForEachAlphaLocalizedPossibility )
                        {
                            if ( std::get<0>(setOfFragments) != possibleSite ) {
                                continue;
                            }

                            auto matchedIons = MatchFragmentIons(theScan, std::get<1>(setOfFragments), commonParameters);
                            double score = CalculatePeptideScore(theScan->getTheScan(), matchedIons, 0);
                            
                            if (score > bestAlphaLocalizedScore)
                            {
                                bestAlphaLocalizedScore = score;
                                bestAlphaSite = possibleSite;
                                bestMatchedAlphaIons = matchedIons;
                            }
                        }
                    }
                    
                    std::vector<std::tuple<int, std::vector<Product*>>> fragmentsForEachBetaLocalizedPossibility =
                        CrosslinkedPeptide::XlGetTheoreticalFragments(
                            commonParameters->getDissociationType(), privateCrosslinker, possibleBetaXlSites,
                            alphaPeptide->getBestPeptide()->getMonoisotopicMass(), betaPeptide->getBestPeptide());
                    
                    std::unordered_set<double> alphaMz;
                    for ( auto p: bestMatchedAlphaIons ) {
                        alphaMz.emplace(p->Mz);
                    }
                    
                    
                    for (auto possibleSite : possibleBetaXlSites)
                    {
                        for ( auto setOfFragments : fragmentsForEachBetaLocalizedPossibility ) 
                        {
                            if ( std::get<0>(setOfFragments) != possibleSite ) {
                                continue;
                            }
                            
                            auto matchedIons = MatchFragmentIons(theScan, std::get<1>(setOfFragments), commonParameters);
                            
                            // remove any matched beta ions that also matched to the alpha peptide
#ifdef ORIG
                            matchedIons.RemoveAll([&] (std::any p) {
                                    std::find(alphaMz.begin(), alphaMz.end(), p::Mz) != alphaMz.end();
                                });
#endif
                            for ( auto p= matchedIons.begin(); p != matchedIons.end(); ) {
                                if ( std::find(alphaMz.begin(), alphaMz.end(), (*p)->Mz) != alphaMz.end() ) {
                                    matchedIons.erase(p);
                                }
                                else {
                                    p++;
                                }
                            }

                            double score = CalculatePeptideScore(theScan->getTheScan(), matchedIons, 0);
                            
                            if (score > bestBetaLocalizedScore)
                            {
                                bestBetaLocalizedScore = score;
                                bestBetaSite = possibleSite;
                                bestMatchedBetaIons = matchedIons;
                            }
                        }
                    }
                    
                    if (bestAlphaLocalizedScore < commonParameters->getScoreCutoff() || bestBetaLocalizedScore < commonParameters->getScoreCutoff())
                    {
                        return nullptr;
                    }
                    
                    auto localizedAlpha = new CrosslinkSpectralMatch(alphaPeptide->getBestPeptide(),
                                                                     alphaPeptide->getBestNotch(),
                                                                     bestAlphaLocalizedScore, 0, theScan,
                                                                     alphaPeptide->getBestPeptide()->getDigestionParams(),
                                                                     bestMatchedAlphaIons);
                    auto localizedBeta = new CrosslinkSpectralMatch(betaPeptide->getBestPeptide(),
                                                                    betaPeptide->getBestNotch(),
                                                                    bestBetaLocalizedScore, 0, theScan,
                                                                    betaPeptide->getBestPeptide()->getDigestionParams(),
                                                                    bestMatchedBetaIons);
                    
                    localizedAlpha->setXlRank(std::vector<int> {ind, inx});
                    localizedAlpha->setXLTotalScore(localizedAlpha->getScore() + localizedBeta->getScore());
                    localizedAlpha->setBetaPeptide(localizedBeta);
                    
                    if (crosslinker->getCleavable())
                    {
                        //TODO: re-enable intensity ranks
                        //psmCrossAlpha.ParentIonMaxIntensityRanks = psmCrossAlpha.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).Select(p => p.IntensityRank).ToList();
                        //localizedAlpha.ParentIonMaxIntensityRanks = new List<int>();
                        
                        //localizedAlpha.ParentIonExistNum = psmCrossAlpha.ParentIonMaxIntensityRanks.Count;
                    }
                    
                    localizedAlpha->setCrossType(PsmCrossType::Cross);
                    localizedCrosslinkedSpectralMatch = localizedAlpha;
                    localizedCrosslinkedSpectralMatch->setLinkPositions(std::vector<int> {bestAlphaSite});
                    localizedCrosslinkedSpectralMatch->getBetaPeptide()->setLinkPositions(std::vector<int> {bestBetaSite});
                    
                    //C# TO C++ CONVERTER TODO TASK: A 'delete localizedBeta' statement was not added since
                    //localizedBeta was assigned to another object. Handle memory management manually.
                    //C# TO C++ CONVERTER TODO TASK: A 'delete localizedAlpha' statement was not added since
                    //localizedAlpha was assigned to an outer scope variable. Handle memory management manually.
                }
            }
            
            return localizedCrosslinkedSpectralMatch;
        }
        
        CrosslinkSpectralMatch *CrosslinkSearchEngine::LocalizeDeadEndSite(PeptideWithSetModifications *originalPeptide,
                                                                           Ms2ScanWithSpecificMass *theScan,
                                                                           CommonParameters *commonParameters,
                                                                           std::vector<int> &possiblePositions,
                                                                           Modification *deadEndMod,
                                                                           int notch, int scanIndex, int peptideIndex)
        {
            double bestScore = 0;
            std::vector<MatchedFragmentIon*> bestMatchingFragments;
            PeptideWithSetModifications *bestLocalizedPeptide = nullptr;
            int bestPosition = 0;
            
            for (auto location : possiblePositions)
            {
                std::unordered_map<int, Modification*> mods = originalPeptide->getAllModsOneIsNterminus();
                
                if (mods.find(location + 1) != mods.end())
                {
                    auto alreadyAnnotatedMod = mods[location + 1];
                    double combinedMass = mods[location + 1]->getMonoisotopicMass().value() +
                        deadEndMod->getMonoisotopicMass().value();

                    std::string acc = "";
                    std::string modType = "";
                    std::string feaType = "";
                    std::string locRestr = "Unassigned.";
                    Chemistry::ChemicalFormula *chemForm = nullptr;
                    Modification *combinedMod = new Modification(alreadyAnnotatedMod->getOriginalId() + "+" + deadEndMod->getOriginalId(),
                                                                 acc, "Crosslink", feaType, alreadyAnnotatedMod->getTarget(),
                                                                 "Anywhere.", chemForm, combinedMass);
                    mods[location + 1] = combinedMod;                    
                }
                else
                {
                    mods.emplace(location + 1, deadEndMod);
                }
                
                auto localizedPeptide = new PeptideWithSetModifications(originalPeptide->getProtein(),
                                                                        originalPeptide->getDigestionParams(),
                                                                        originalPeptide->getOneBasedStartResidueInProtein(),
                                                                        originalPeptide->getOneBasedEndResidueInProtein(),
                                                                        originalPeptide->getCleavageSpecificityForFdrCategory(),
                                                                        originalPeptide->getPeptideDescription(),
                                                                        originalPeptide->getMissedCleavages(),
                                                                        mods,
                                                                        originalPeptide->NumFixedMods);
                
                auto products = localizedPeptide->Fragment(commonParameters->getDissociationType(), FragmentationTerminus::Both);//.ToList();
                auto matchedFragmentIons = MatchFragmentIons(theScan, products, commonParameters);
                
                double score = CalculatePeptideScore(theScan->getTheScan(), matchedFragmentIons, 0);
                
                if (score > bestScore)
                {
                    bestMatchingFragments = matchedFragmentIons;
                    bestScore = score;
                    if (  bestLocalizedPeptide != nullptr ) {
                        delete  bestLocalizedPeptide;
                    }
                    bestLocalizedPeptide = localizedPeptide;
                    bestPosition = location;
                }
            }
            
            if (bestScore < commonParameters->getScoreCutoff())
            {
                return nullptr;
            }
            
            auto csm = new CrosslinkSpectralMatch(bestLocalizedPeptide, notch, bestScore, scanIndex, theScan,
                                                  originalPeptide->getDigestionParams(), bestMatchingFragments);
            
            if (deadEndMod == TrisDeadEnd)
            {
                csm->setCrossType(PsmCrossType::DeadEndTris);
            }
            else if (deadEndMod == H2ODeadEnd)
            {
                csm->setCrossType(PsmCrossType::DeadEndH2O);
            }
            else if (deadEndMod == NH2DeadEnd)
            {
                csm->setCrossType(PsmCrossType::DeadEndNH2);
            }
            
            csm->setLinkPositions(std::vector<int> {bestPosition});
            csm->setXlRank(std::vector<int> {peptideIndex});
            
            return csm;
        }

        CrosslinkSpectralMatch *CrosslinkSearchEngine::LocalizeLoopSites(PeptideWithSetModifications *originalPeptide, Ms2ScanWithSpecificMass *theScan,
                                                                         CommonParameters *commonParameters, std::vector<int> &possiblePositions,
                                                                         Modification *loopMod, int notch, int scanIndex, int peptideIndex)
        {
            XLumap possibleFragmentSets = CrosslinkedPeptide::XlLoopGetTheoreticalFragments(
                commonParameters->getDissociationType(),
                Loop, possiblePositions, originalPeptide);
            double bestScore = 0;
            std::tuple<int, int> bestModPositionSites;
            std::vector<MatchedFragmentIon*> bestMatchingFragments;
            
            for (auto setOfPositions : possibleFragmentSets)
            {
                auto matchedFragmentIons = MatchFragmentIons(theScan, std::get<1>(setOfPositions), commonParameters);
                
                double score = CalculatePeptideScore(theScan->getTheScan(), matchedFragmentIons, 0);
                
                if (score > bestScore)
                {
                    bestMatchingFragments = matchedFragmentIons;
                    bestScore = score;
                    bestModPositionSites = std::get<0>(setOfPositions);
                }
            }
            
            if (bestScore < commonParameters->getScoreCutoff())
            {
                return nullptr;
            }
            
            auto csm = new CrosslinkSpectralMatch(originalPeptide, notch, bestScore, scanIndex, theScan,
                                                  originalPeptide->getDigestionParams(), bestMatchingFragments);
            csm->setCrossType(PsmCrossType::Loop);
            csm->setXlRank(std::vector<int> {peptideIndex});
            csm->setLinkPositions(std::vector<int> {std::get<0>(bestModPositionSites), std::get<1>(bestModPositionSites)});
            
            return csm;
        }
    }
}
