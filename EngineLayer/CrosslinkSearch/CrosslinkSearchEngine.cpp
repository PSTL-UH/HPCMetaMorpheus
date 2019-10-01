#include "CrosslinkSearchEngine.h"
#include "CrosslinkSpectralMatch.h"
#include "Crosslinker.h"
#include "../PrecursorSearchModes/MassDiffAcceptor.h"
#include "../CommonParameters.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../PrecursorSearchModes/SingleAbsoluteAroundZeroMassDiffAcceptor.h"
#include "../MetaMorpheusEngineResults.h"
#include "../EventArgs/ProgressEventArgs.h"
#include "BestPeptideScoreNotch.h"
#include "PsmCrossType.h"
#include "CrosslinkedPeptides.h"

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
                                                     std::vector<std::vector<int>&> &fragmentIndex,
                                                     int currentPartition, CommonParameters *commonParameters,
                                                     Crosslinker *crosslinker,
                                                     bool CrosslinkSearchTop,
                                                     int CrosslinkSearchTopNum,
                                                     bool quench_H2O,
                                                     bool quench_NH2,
                                                     bool quench_Tris,
                                                     std::vector<std::string> &nestedIds) : ModernSearchEngine(),
                                                                                            GlobalCsms(globalCsms),
                                                                                            Crosslinker(crosslinker),
                                                                                            CrosslinkSearchTopN(CrosslinkSearchTop),
                                                                                            TopN(CrosslinkSearchTopNum),
                                                                                            QuenchH2O(quench_H2O),
                                                                                            QuenchNH2(quench_NH2),
                                                                                            QuenchTris(quench_Tris)
        {
            GenerateCrosslinkModifications(crosslinker);
            AllCrosslinkerSites = Crosslinker->getCrosslinkerModSites().ToCharArray().Concat(Crosslinker->getCrosslinkerModSites2().ToCharArray())->Distinct()->ToArray();
            
            if (dynamic_cast<PpmTolerance*>(commonParameters->getPrecursorMassTolerance()) != nullptr)
            {
                XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(commonParameters->getPrecursorMassTolerance().value());
            }
            else
            {
                XLPrecusorSearchMode = new SingleAbsoluteAroundZeroSearchMode(commonParameters->getPrecursorMassTolerance().value());
            }
        }
        
        MetaMorpheusEngineResults *CrosslinkSearchEngine::RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ProgressEventArgs tempVar(oldPercentProgress, "Performing crosslink search... " + std::to_string(CurrentPartition) + "/" + std::to_string(commonParameters->getTotalPartitions()), nestedIds);
            ReportProgress(&tempVar);
            
            unsigned char byteScoreCutoff = static_cast<unsigned char>(commonParameters->getScoreCutoff());
            
            ParallelOptions *tempVar2 = new ParallelOptions();
            tempVar2->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
            Parallel::ForEach(Partitioner::Create(0, ListOfSortedMs2Scans.size()), tempVar2, [&] (range, loopState)  {
                    std::vector<unsigned char> scoringTable(PeptideIndex.size());
                    std::vector<int> idsOfPeptidesPossiblyObserved;
                    
                    for (int scanIndex = range::Item1; scanIndex < range::Item2; scanIndex++)
                    {
                        // Stop loop if canceled
                        if (GlobalVariables::getStopLoops())
                        {
                            loopState::Stop();
                            return;
                        }
                        
                        // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                        Array::Clear(scoringTable, 0, scoringTable.Length);
                        idsOfPeptidesPossiblyObserved.Clear();
                        auto scan = ListOfSortedMs2Scans[scanIndex];
                        
                        // get fragment bins for this scan
                        std::vector<int> allBinsToSearch = GetBinsToSearch(scan);
                        std::vector<BestPeptideScoreNotch*> bestPeptideScoreNotchList;
                        
                        // first-pass scoring
                        IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan->getPrecursorMass(), -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), PeptideIndex, MassDiffAcceptor, 0);
                        
                        // done with indexed scoring; refine scores and create PSMs
                        if (idsOfPeptidesPossiblyObserved.Any())
                        {
                            if (CrosslinkSearchTopN)
                            {
                                // take top N hits for this scan
                                idsOfPeptidesPossiblyObserved = idsOfPeptidesPossiblyObserved.OrderByDescending([&] (std::any p) {
                                        scoringTable[p];
                                    }).Take(TopN)->ToList();
                            }
                            
                            for (auto id : idsOfPeptidesPossiblyObserved)
                            {
                                PeptideWithSetModifications *peptide = PeptideIndex[id];
                                
                                int notch = MassDiffAcceptor->Accepts(scan->getPrecursorMass(), peptide->MonoisotopicMass);
                                BestPeptideScoreNotch tempVar3(peptide, scoringTable[id], notch);
                                bestPeptideScoreNotchList.Add(&tempVar3);
                            }
                            
                            // combine individual peptide hits with crosslinker mass to find best crosslink PSM hit
                            auto csm = FindCrosslinkedPeptide(scan, bestPeptideScoreNotchList, scanIndex);
                            
                            if (csm == nullptr)
                            {
                                progress++;
                                continue;
                            }
                            
                            // this scan might already have a hit from a different database partition; check to see if the score improves
                            if (GlobalCsms[scanIndex] == nullptr || GlobalCsms[scanIndex]->getXLTotalScore() < csm->getXLTotalScore())
                            {
                                GlobalCsms[scanIndex] = csm;
                            }
                        }
                        
                        // report search progress
                        progress++;
                        auto percentProgress = static_cast<int>((progress / ListOfSortedMs2Scans.size()) * 100);
                        
                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ProgressEventArgs tempVar4(percentProgress, "Performing crosslink search... " + std::to_string(CurrentPartition) + "/" + std::to_string(commonParameters->getTotalPartitions()), nestedIds);
                            ReportProgress(&tempVar4);
                        }
                    }
                });
            
            return new MetaMorpheusEngineResults(this);
        }
        
        void CrosslinkSearchEngine::GenerateCrosslinkModifications(Crosslinker *crosslinker)
        {
            std::any motif;
            ModificationMotif::TryGetMotif("X", motif);
            TrisDeadEnd = new Modification(_originalId: "Tris Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker->getDeadendMassTris());
            H2ODeadEnd = new Modification(_originalId: "H2O Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker->getDeadendMassH2O());
            NH2DeadEnd = new Modification(_originalId: "NH2 Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker->getDeadendMassNH2());
            Loop = new Modification(_originalId: "Loop", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker->getLoopMass());
        }
        
        CrosslinkSpectralMatch *CrosslinkSearchEngine::FindCrosslinkedPeptide(Ms2ScanWithSpecificMass *theScan, std::vector<BestPeptideScoreNotch*> &theScanBestPeptide, int scanIndex)
        {
            std::vector<CrosslinkSpectralMatch*> possibleMatches;
            
            for (int alphaIndex = 0; alphaIndex < theScanBestPeptide.size(); alphaIndex++)
            {
                PeptideWithSetModifications *bestPeptide = theScanBestPeptide[alphaIndex]->getBestPeptide();
                
                //Single Peptide
                if (XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), bestPeptide->MonoisotopicMass) >= 0)
                {
                    std::vector<Product*> products = bestPeptide->Fragment(commonParameters->getDissociationType(), FragmentationTerminus::Both).ToList();
                    auto matchedFragmentIons = MatchFragmentIons(theScan, products, commonParameters);
                    double score = CalculatePeptideScore(theScan->getTheScan(), matchedFragmentIons, 0);
                    
                    auto psmCrossSingle = new CrosslinkSpectralMatch(bestPeptide, theScanBestPeptide[alphaIndex]->getBestNotch(), score, scanIndex, theScan, commonParameters->getDigestionParams(), matchedFragmentIons);
                    psmCrossSingle->setCrossType(PsmCrossType::Single);
                    psmCrossSingle->setXlRank(std::vector<int> {alphaIndex});
                    
                    possibleMatches.push_back(psmCrossSingle);
                    
                    //C# TO C++ CONVERTER TODO TASK: A 'delete psmCrossSingle' statement was not added since
                    //psmCrossSingle was passed to a method or constructor. Handle memory management manually.
                }
                // Deadend Peptide
                else if (QuenchTris && XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), bestPeptide->MonoisotopicMass + Crosslinker->getDeadendMassTris()) >= 0)
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    
                    if (possibleCrosslinkLocations.Any())
                    {
                        // tris deadend
                        possibleMatches.push_back(LocalizeDeadEndSite(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, TrisDeadEnd, theScanBestPeptide[alphaIndex]->getBestNotch(), scanIndex, alphaIndex));
                    }
                }
                else if (QuenchH2O && XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), bestPeptide->MonoisotopicMass + Crosslinker->getDeadendMassH2O()) >= 0)
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    
                    if (possibleCrosslinkLocations.Any())
                    {
                        // H2O deadend
                        possibleMatches.push_back(LocalizeDeadEndSite(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, H2ODeadEnd, theScanBestPeptide[alphaIndex]->getBestNotch(), scanIndex, alphaIndex));
                    }
                }
                else if (QuenchNH2 && XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), bestPeptide->MonoisotopicMass + Crosslinker->getDeadendMassNH2()) >= 0)
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    
                    if (possibleCrosslinkLocations.Any())
                    {
                        // NH2 deadend
                        possibleMatches.push_back(LocalizeDeadEndSite(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, NH2DeadEnd, theScanBestPeptide[alphaIndex]->getBestNotch(), scanIndex, alphaIndex));
                    }
                }
                // loop peptide
                else if (Crosslinker->getLoopMass() != 0 && XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), bestPeptide->MonoisotopicMass + Crosslinker->getLoopMass()) >= 0)
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    
                    if (possibleCrosslinkLocations.size() >= 2)
                    {
                        possibleMatches.push_back(LocalizeLoopSites(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, Loop, theScanBestPeptide[alphaIndex]->getBestNotch(), scanIndex, alphaIndex));
                    }
                }
                // Cross-linked peptide
                else if (theScan->getPrecursorMass() - bestPeptide->MonoisotopicMass >= (commonParameters->getDigestionParams()->MinPeptideLength * 50))
                {
                    std::vector<int> possibleCrosslinkLocations = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    if (!possibleCrosslinkLocations.Any())
                    {
                        continue;
                    }
                    
                    PeptideWithSetModifications *alphaPeptide = bestPeptide;
                    
                    for (int betaIndex = 0; betaIndex < theScanBestPeptide.size(); betaIndex++)
                    {
                        PeptideWithSetModifications *betaPeptide = theScanBestPeptide[betaIndex]->getBestPeptide();
                        
                        if (XLPrecusorSearchMode->Accepts(theScan->getPrecursorMass(), alphaPeptide->MonoisotopicMass + betaPeptide->MonoisotopicMass + Crosslinker->getTotalMass()) >= 0)
                        {
                            std::vector<int> possibleBetaCrosslinkSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(AllCrosslinkerSites, betaPeptide);
                            
                            if (!possibleBetaCrosslinkSites.Any())
                            {
                                continue;
                            }
                            
                            CrosslinkSpectralMatch *csm = LocalizeCrosslinkSites(theScan, theScanBestPeptide[alphaIndex], theScanBestPeptide[betaIndex], Crosslinker, alphaIndex, betaIndex);
                            
                            possibleMatches.push_back(csm);
                        }
                    }
                }
            }
            
            // get the best match for this spectrum
            // bestPsmCross will be null if there are no valid hits
            possibleMatches.RemoveAll([&] (std::any v)  {
                    return v == nullptr;
                });
            possibleMatches = possibleMatches.OrderByDescending([&] (std::any p)  {
                    p::XLTotalScore;
                }).ToList();
            auto bestPsmCross = possibleMatches.FirstOrDefault();
            
            // resolve ambiguities
            if (bestPsmCross != nullptr)
            {
                bestPsmCross->ResolveAllAmbiguities();
                
                if (bestPsmCross->BetaPeptide != nullptr)
                {
                    bestPsmCross->BetaPeptide.ResolveAllAmbiguities();
                }
            }
            
            // calculate delta score
            if (possibleMatches.size() > 1)
            {
                bestPsmCross->DeltaScore = possibleMatches[0]->getXLTotalScore() - possibleMatches[1]->getXLTotalScore();
            }
            
            return bestPsmCross;
        }
        
        CrosslinkSpectralMatch *CrosslinkSearchEngine::LocalizeCrosslinkSites(Ms2ScanWithSpecificMass *theScan, BestPeptideScoreNotch *alphaPeptide, BestPeptideScoreNotch *betaPeptide, Crosslinker *crosslinker, int ind, int inx)
        {
            CrosslinkSpectralMatch *localizedCrosslinkedSpectralMatch = nullptr;
            
            std::vector<std::tuple<std::vector<int>, std::vector<int>>> pairs;
            
            if (crosslinker->getCrosslinkerModSites() == crosslinker->getCrosslinkerModSites2())
            {
                std::vector<int> possibleAlphaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(crosslinker->getCrosslinkerModSites().ToCharArray(), alphaPeptide->getBestPeptide());
                std::vector<int> possibleBetaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(crosslinker->getCrosslinkerModSites().ToCharArray(), betaPeptide->getBestPeptide());
                
                pairs.push_back(std::tuple<std::vector<int>, std::vector<int>>(possibleAlphaXlSites, possibleBetaXlSites));
            }
            else
            {
                std::vector<int> possibleAlphaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(crosslinker->getCrosslinkerModSites().ToCharArray(), alphaPeptide->getBestPeptide());
                std::vector<int> possibleBetaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(crosslinker->getCrosslinkerModSites2().ToCharArray(), betaPeptide->getBestPeptide());
                
                pairs.push_back(std::tuple<std::vector<int>, std::vector<int>>(possibleAlphaXlSites, possibleBetaXlSites));
                
                possibleAlphaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(crosslinker->getCrosslinkerModSites2().ToCharArray(), alphaPeptide->getBestPeptide());
                possibleBetaXlSites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(crosslinker->getCrosslinkerModSites().ToCharArray(), betaPeptide->getBestPeptide());
                
                pairs.push_back(std::tuple<std::vector<int>, std::vector<int>>(possibleAlphaXlSites, possibleBetaXlSites));
            }
            
            for (auto pair : pairs)
            {
                std::vector<int> possibleAlphaXlSites = pair.Item1;
                std::vector<int> possibleBetaXlSites = pair.Item2;
                
                if (possibleAlphaXlSites.Any() && possibleBetaXlSites.Any())
                {
                    int bestAlphaSite = 0;
                    int bestBetaSite = 0;
                    std::vector<MatchedFragmentIon*> bestMatchedAlphaIons;
                    std::vector<MatchedFragmentIon*> bestMatchedBetaIons;
                    double bestAlphaLocalizedScore = 0;
                    double bestBetaLocalizedScore = 0;
                    
                    auto fragmentsForEachAlphaLocalizedPossibility = CrosslinkedPeptide::XlGetTheoreticalFragments(commonParameters->getDissociationType(), Crosslinker, possibleAlphaXlSites, betaPeptide->getBestPeptide()->MonoisotopicMass, alphaPeptide->getBestPeptide()).ToList();
                    
                    for (auto possibleSite : possibleAlphaXlSites)
                    {
                        for (auto setOfFragments : fragmentsForEachAlphaLocalizedPossibility.Where([&] (std::any v) {
                                    return v->Item1 == possibleSite;
                                }))
                        {
                            auto matchedIons = MatchFragmentIons(theScan, setOfFragments::Item2, commonParameters);
                            double score = CalculatePeptideScore(theScan->getTheScan(), matchedIons, 0);
                            
                            if (score > bestAlphaLocalizedScore)
                            {
                                bestAlphaLocalizedScore = score;
                                bestAlphaSite = possibleSite;
                                bestMatchedAlphaIons = matchedIons;
                            }
                        }
                    }
                    
                    auto fragmentsForEachBetaLocalizedPossibility = CrosslinkedPeptide::XlGetTheoreticalFragments(commonParameters->getDissociationType(), Crosslinker, possibleBetaXlSites, alphaPeptide->getBestPeptide()->MonoisotopicMass, betaPeptide->getBestPeptide()).ToList();
                    
                    auto alphaMz = std::unordered_set<double>(bestMatchedAlphaIons.Select([&] (std::any p)
                    {
                        p::Mz;
                    }));
                    
                    for (auto possibleSite : possibleBetaXlSites)
                    {
                        for (auto setOfFragments : fragmentsForEachBetaLocalizedPossibility.Where([&] (std::any v) {
                                    return v->Item1 == possibleSite;
                                }))
                        {
                            auto matchedIons = MatchFragmentIons(theScan, setOfFragments::Item2, commonParameters);
                            
                            // remove any matched beta ions that also matched to the alpha peptide
                            matchedIons.RemoveAll([&] (std::any p) {
                                    std::find(alphaMz.begin(), alphaMz.end(), p::Mz) != alphaMz.end();
                                });
                            
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
                    
                    auto localizedAlpha = new CrosslinkSpectralMatch(alphaPeptide->getBestPeptide(), alphaPeptide->getBestNotch(), bestAlphaLocalizedScore, 0, theScan, alphaPeptide->getBestPeptide()->DigestionParams, bestMatchedAlphaIons);
                    auto localizedBeta = new CrosslinkSpectralMatch(betaPeptide->getBestPeptide(), betaPeptide->getBestNotch(), bestBetaLocalizedScore, 0, theScan, betaPeptide->getBestPeptide()->DigestionParams, bestMatchedBetaIons);
                    
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
        
        CrosslinkSpectralMatch *CrosslinkSearchEngine::LocalizeDeadEndSite(PeptideWithSetModifications *originalPeptide, Ms2ScanWithSpecificMass *theScan, CommonParameters *commonParameters, std::vector<int> &possiblePositions, Modification *deadEndMod, int notch, int scanIndex, int peptideIndex)
        {
            double bestScore = 0;
            std::vector<MatchedFragmentIon*> bestMatchingFragments;
            PeptideWithSetModifications *bestLocalizedPeptide = nullptr;
            int bestPosition = 0;
            
            for (auto location : possiblePositions)
            {
                std::unordered_map<int, Modification*> mods = originalPeptide->AllModsOneIsNterminus.ToDictionary([&] (std::any p) {
                        p::Key;
                    }, [&] (std::any p)
                    {
                        p->Value;
                    });
                if (mods.find(location + 1) != mods.end())
                {
                    auto alreadyAnnotatedMod = mods[location + 1];
                    double combinedMass = mods[location + 1]->MonoisotopicMass->Value + deadEndMod->MonoisotopicMass->Value;
                    Modification *combinedMod = new Modification(_originalId: alreadyAnnotatedMod->OriginalId + "+" + deadEndMod->OriginalId, _modificationType: "Crosslink", _target: alreadyAnnotatedMod->Target, _locationRestriction: "Anywhere.", _monoisotopicMass: combinedMass);
                    mods[location + 1] = combinedMod;
                    
                    delete combinedMod;
                }
                else
                {
                    mods.emplace(location + 1, deadEndMod);
                }
                
                auto localizedPeptide = new PeptideWithSetModifications(originalPeptide->Protein, originalPeptide->DigestionParams, originalPeptide->OneBasedStartResidueInProtein, originalPeptide->OneBasedEndResidueInProtein, originalPeptide->CleavageSpecificityForFdrCategory, originalPeptide->PeptideDescription, originalPeptide->MissedCleavages, mods, originalPeptide->NumFixedMods);
                
                auto products = localizedPeptide->Fragment(commonParameters->getDissociationType(), FragmentationTerminus::Both).ToList();
                auto matchedFragmentIons = MatchFragmentIons(theScan, products, commonParameters);
                
                double score = CalculatePeptideScore(theScan->getTheScan(), matchedFragmentIons, 0);
                
                if (score > bestScore)
                {
                    bestMatchingFragments = matchedFragmentIons;
                    bestScore = score;
                    bestLocalizedPeptide = localizedPeptide;
                    bestPosition = location;
                }
                
                //C# TO C++ CONVERTER TODO TASK: A 'delete localizedPeptide' statement was not added since
                //localizedPeptide was assigned to an outer scope variable. Handle memory management manually.
            }
            
            if (bestScore < commonParameters->getScoreCutoff())
            {
                return nullptr;
            }
            
            auto csm = new CrosslinkSpectralMatch(bestLocalizedPeptide, notch, bestScore, scanIndex, theScan, originalPeptide->DigestionParams, bestMatchingFragments);
            
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
            
            //C# TO C++ CONVERTER TODO TASK: A 'delete csm' statement was not added since
            //csm was used in a 'return' or 'throw' statement.
            return csm;
        }

        CrosslinkSpectralMatch *CrosslinkSearchEngine::LocalizeLoopSites(PeptideWithSetModifications *originalPeptide, Ms2ScanWithSpecificMass *theScan, CommonParameters *commonParameters, std::vector<int> &possiblePositions, Modification *loopMod, int notch, int scanIndex, int peptideIndex)
        {
            auto possibleFragmentSets = CrosslinkedPeptide::XlLoopGetTheoreticalFragments(commonParameters->getDissociationType(), Loop, possiblePositions, originalPeptide);
            double bestScore = 0;
            std::tuple<int, int> bestModPositionSites = nullptr;
            std::vector<MatchedFragmentIon*> bestMatchingFragments;
            
            for (auto setOfPositions : possibleFragmentSets)
            {
                auto matchedFragmentIons = MatchFragmentIons(theScan, setOfPositions.Value, commonParameters);
                
                double score = CalculatePeptideScore(theScan->getTheScan(), matchedFragmentIons, 0);
                
                if (score > bestScore)
                {
                    bestMatchingFragments = matchedFragmentIons;
                    bestScore = score;
                    bestModPositionSites = setOfPositions.Key;
                }
            }
            
            if (bestScore < commonParameters->getScoreCutoff())
            {
                return nullptr;
            }
            
            auto csm = new CrosslinkSpectralMatch(originalPeptide, notch, bestScore, scanIndex, theScan, originalPeptide->DigestionParams, bestMatchingFragments);
            csm->setCrossType(PsmCrossType::Loop);
            csm->setXlRank(std::vector<int> {peptideIndex});
            csm->setLinkPositions(std::vector<int> {std::get<0>(bestModPositionSites), std::get<1>(bestModPositionSites)});
            
            //C# TO C++ CONVERTER TODO TASK: A 'delete csm' statement was not added since csm was used in a 'return' or 'throw' statement.
            return csm;
        }
    }
}
