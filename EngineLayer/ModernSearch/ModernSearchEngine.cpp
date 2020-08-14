#include "ModernSearchEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../PrecursorSearchModes/MassDiffAcceptor.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "../EventArgs/ProgressEventArgs.h"
#include "../GlobalVariables.h"
#include "bankersrounding.h"

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
    namespace ModernSearch
    {
        
        ModernSearchEngine::ModernSearchEngine(std::vector<PeptideSpectralMatch*> globalPsms,
                                               std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans,
                                               std::vector<PeptideWithSetModifications*> &peptideIndex,
                                               std::vector<std::vector<int>> &fragmentIndex,
                                               int currentPartition,
                                               CommonParameters *commonParameters,
                                               MassDiffAcceptor *massDiffAcceptor,
                                               double maximumMassThatFragmentIonScoreIsDoubled,
                                               std::vector<std::string> &nestedIds) :
            MetaMorpheusEngine(commonParameters, nestedIds),
            FragmentIndex(fragmentIndex),
            PeptideSpectralMatches(globalPsms),
            ListOfSortedMs2Scans(listOfSortedms2Scans),
            PeptideIndex(peptideIndex),
            CurrentPartition(currentPartition + 1),
            massDiffAcceptor(massDiffAcceptor),
            dissociationType(commonParameters->getDissociationType()),
            MaxMassThatFragmentIonScoreIsDoubled(maximumMassThatFragmentIonScoreIsDoubled)
        {
        }
        
        MetaMorpheusEngineResults *ModernSearchEngine::RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ProgressEventArgs tempVar(oldPercentProgress, "Performing modern search... " + std::to_string(CurrentPartition) + "/" +
                                      std::to_string(commonParameters->getTotalPartitions()),
                                      const_cast<std::vector<std::string>&>(nestedIds));
            ReportProgress(&tempVar);
            
            unsigned char byteScoreCutoff = static_cast<unsigned char>(commonParameters->getScoreCutoff());
            if (commonParameters->getCalculateEValue())
            {
                byteScoreCutoff = 1;
            }
            
            //ParallelOptions *tempVar2 = new ParallelOptions();
            //tempVar2->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
            //Parallel::ForEach(Partitioner::Create(0, ListOfSortedMs2Scans.size()), tempVar2, [&] (range, loopState) {

            std::vector<unsigned char> scoringTable(PeptideIndex.size());
            std::vector<int> idsOfPeptidesPossiblyObserved;
            for (int i = 0; i<(int)ListOfSortedMs2Scans.size(); i++ ) {
                
                //for (int i = range::Item1; i < range::Item2; i++)
                //{
                // Stop loop if canceled
                if (GlobalVariables::getStopLoops())
                {
                    //loopState::Stop();
                    return nullptr;
                }
                        
                //empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                //Array::Clear(scoringTable, 0, scoringTable.Length);
                scoringTable.clear();
                idsOfPeptidesPossiblyObserved.clear();
                Ms2ScanWithSpecificMass *scan = ListOfSortedMs2Scans[i];
                
                // get fragment bins for this scan
                std::vector<int> allBinsToSearch = GetBinsToSearch(scan);
                        
                // get allowed theoretical masses from the known experimental mass
                // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)
                // this is just PRELIMINARY precursor-mass filtering
                // additional checks are made later to ensure that the theoretical precursor mass is acceptable
                std::vector<AllowedIntervalWithNotch*> notches = massDiffAcceptor->GetAllowedPrecursorMassIntervalsFromObservedMass(scan->getPrecursorMass());
                        
#ifdef ORIG
                double lowestMassPeptideToLookFor = notches.Min([&] (std::any p)  {
                        p::AllowedInterval::Minimum;
                    });
                double highestMassPeptideToLookFor = notches.Max([&] (std::any p) {
                        p::AllowedInterval::Maximum;
                    });
#endif
                double lowestMassPeptideToLookFor = std::numeric_limits<double>::max();
                double highestMassPeptideToLookFor = std::numeric_limits<double>::min();
                for ( auto p: notches ) {
                    if ( p->AllowedInterval->getMinimum() < lowestMassPeptideToLookFor ) {
                        lowestMassPeptideToLookFor = p->AllowedInterval->getMinimum();
                    }
                    if ( p->AllowedInterval->getMaximum() < highestMassPeptideToLookFor ) {
                        highestMassPeptideToLookFor = p->AllowedInterval->getMaximum();
                    }
                }
                
                // first-pass scoring
                IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved,
                               scan->getPrecursorMass(), lowestMassPeptideToLookFor, highestMassPeptideToLookFor,
                               PeptideIndex, massDiffAcceptor, MaxMassThatFragmentIonScoreIsDoubled);
                
                // done with indexed scoring; refine scores and create PSMs
                for (auto id : idsOfPeptidesPossiblyObserved)
                {
                    PeptideWithSetModifications *peptide = PeptideIndex[id];
                    
                    std::vector<Product*> peptideTheorProducts = peptide->Fragment(commonParameters->getDissociationType(),
                                                                                   FragmentationTerminus::Both); //.ToList();
                    
                    std::vector<MatchedFragmentIon*> matchedIons = MatchFragmentIons(scan, peptideTheorProducts, commonParameters);
                            
                    double thisScore = CalculatePeptideScore(scan->getTheScan(), matchedIons, 0);
                    int notch = massDiffAcceptor->Accepts(scan->getPrecursorMass(), peptide->getMonoisotopicMass());
                    
                    bool meetsScoreCutoff = thisScore >= commonParameters->getScoreCutoff();
                    bool scoreImprovement = PeptideSpectralMatches[i] == nullptr ||
                        (thisScore - PeptideSpectralMatches[i]->getRunnerUpScore()) > -PeptideSpectralMatch::ToleranceForScoreDifferentiation;
                    
                    if ( (meetsScoreCutoff && scoreImprovement) || commonParameters->getCalculateEValue())
                    {
                        if (PeptideSpectralMatches[i] == nullptr)
                        {
                            PeptideSpectralMatches[i] = new PeptideSpectralMatch(peptide, notch, thisScore, i, scan,
                                                                                 commonParameters->getDigestionParams(), matchedIons);
                        }
                        else
                        {
                            PeptideSpectralMatches[i]->AddOrReplace(peptide, thisScore, notch,
                                                                    commonParameters->getReportAllAmbiguity(), matchedIons);
                        }
                        
                        if (commonParameters->getCalculateEValue())
                        {
                            PeptideSpectralMatches[i]->getAllScores().push_back(thisScore);
                        }
                    }
                }
                
                // report search progress
                progress++;
                auto percentProgress = static_cast<int>((progress / ListOfSortedMs2Scans.size()) * 100);
                
                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ProgressEventArgs tempVar3(percentProgress, "Performing modern search... " + std::to_string(CurrentPartition) + "/" +
                                               std::to_string(commonParameters->getTotalPartitions()),
                                               const_cast<std::vector<std::string>&> (nestedIds));
                    ReportProgress(&tempVar3);
                }
            }
            //});
            
            // remove peptides below the score cutoff that were stored to calculate expectation values
            if (commonParameters->getCalculateEValue())
            {
                for (int i = 0; i < (int)PeptideSpectralMatches.size(); i++)
                {
                    if (PeptideSpectralMatches[i] != nullptr && PeptideSpectralMatches[i]->getScore() < commonParameters->getScoreCutoff())
                    {
                        PeptideSpectralMatches[i] = nullptr;
                    }
                }
            }
            
#ifdef ORIG
            for (PeptideSpectralMatch *psm : PeptideSpectralMatches.Where([&] (std::any p) {
                        return p != nullptr;
                    }))
            {
                psm::ResolveAllAmbiguities();
            }
#endif
            for (PeptideSpectralMatch *psm : PeptideSpectralMatches ) {
                if ( psm != nullptr ) {
                    psm->ResolveAllAmbiguities();
                }
            }
            
            return new MetaMorpheusEngineResults(this);
        }
        
        std::vector<int> ModernSearchEngine::GetBinsToSearch(Ms2ScanWithSpecificMass *scan)
        {
            int obsPreviousFragmentCeilingMz = 0;
            std::vector<int> binsToSearch;

            for (auto envelope : scan->getExperimentalFragments())
            {
                // assume charge state 1 to calculate mass tolerance
                double experimentalFragmentMass = envelope->monoisotopicMass;
                
                // get theoretical fragment bins within mass tolerance
                int obsFragmentFloorMass = static_cast<int>(std::floor((commonParameters->getProductMassTolerance()->GetMinimumValue(experimentalFragmentMass)) * FragmentBinsPerDalton));
                int obsFragmentCeilingMass = static_cast<int>(std::ceil((commonParameters->getProductMassTolerance()->GetMaximumValue(experimentalFragmentMass)) * FragmentBinsPerDalton));
                
                // prevents double-counting peaks close in m/z and lower-bound out of range exceptions
                if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)
                {
                    obsFragmentFloorMass = obsPreviousFragmentCeilingMz;
                }
                obsPreviousFragmentCeilingMz = obsFragmentCeilingMass + 1;
                
                // prevent upper-bound index out of bounds errors;
                // lower-bound is handled by the previous "if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)" statement
                if (obsFragmentCeilingMass >= (int)FragmentIndex.size())
                {
                    obsFragmentCeilingMass = (int)FragmentIndex.size() - 1;
                    
                    if (obsFragmentFloorMass >= (int)FragmentIndex.size())
                    {
                        obsFragmentFloorMass = FragmentIndex.size() - 1;
                    }
                }
                
                // search mass bins within a tolerance
                for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                {
                    if (FragmentIndex[fragmentBin].size() > 0)
                    {
                        binsToSearch.push_back(fragmentBin);
                    }
                }
                
                // add complementary ions
                if (commonParameters->getAddCompIons())
                {
                    //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut
                    //and just straight up adding the bins assuming that they're z=1
                    
                    double protonMassShift;
                    std::unordered_map<DissociationType, double>::const_iterator complementaryIonConversionDictionary_iterator = complementaryIonConversionDictionary.find(commonParameters->getDissociationType());
                    if (complementaryIonConversionDictionary_iterator != complementaryIonConversionDictionary.end()) //TODO: this is broken for EThcD because that method needs two conversions
                    {
                        protonMassShift = complementaryIonConversionDictionary_iterator->second;
                        protonMassShift = Chemistry::ClassExtensions::ToMass(protonMassShift, 1);
                        int compFragmentFloorMass = static_cast<int>(BankersRounding::round(((scan->getPrecursorMass() + protonMassShift) * FragmentBinsPerDalton))) - obsFragmentCeilingMass;
                        int compFragmentCeilingMass = static_cast<int>(BankersRounding::round(((scan->getPrecursorMass() + protonMassShift) * FragmentBinsPerDalton))) - obsFragmentFloorMass;
                        
                        // prevent index out of bounds errors
                        if (compFragmentCeilingMass >= (int)FragmentIndex.size())
                        {
                            compFragmentCeilingMass = FragmentIndex.size() - 1;
                            
                            if (compFragmentFloorMass >= (int)FragmentIndex.size())
                            {
                                compFragmentFloorMass = FragmentIndex.size() - 1;
                            }
                        }
                        if (compFragmentFloorMass < 0)
                        {
                            compFragmentFloorMass = 0;
                        }
                        
                        for (int fragmentBin = compFragmentFloorMass; fragmentBin <= compFragmentCeilingMass; fragmentBin++)
                        {
                            if (FragmentIndex[fragmentBin].size() > 0)
                            {
                                binsToSearch.push_back(fragmentBin);
                            }
                        }
                    }
                    else
                    {
                        protonMassShift = complementaryIonConversionDictionary_iterator->second;
                        throw NotImplementedException();
                    }
                }
            }
            return binsToSearch;
        }
        
        int ModernSearchEngine::BinarySearchBinForPrecursorIndex(std::vector<int> &peptideIdsInThisBin, double peptideMassToLookFor, std::vector<PeptideWithSetModifications*> &peptideIndex)
        {
            int m = 0;
            int l = 0;
            int r = peptideIdsInThisBin.size() - 1;
            
            // binary search in the fragment bin for precursor mass
            while (l <= r)
            {
                m = l + ((r - l) / 2);
                
                if (r - l < 2)
                {
                    break;
                }
                if (peptideIndex[peptideIdsInThisBin[m]]->getMonoisotopicMass() < peptideMassToLookFor)
                {
                    l = m + 1;
                }
                else
                {
                    r = m - 1;
                }
            }
            if (m > 0)
            {
                m--;
            }
            return m;
        }
        
        void ModernSearchEngine::IndexedScoring(std::vector<int> &binsToSearch, std::vector<unsigned char> &scoringTable, unsigned char byteScoreCutoff, std::vector<int> &idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor, std::vector<PeptideWithSetModifications*> &peptideIndex, MassDiffAcceptor *massDiffAcceptor, double maxMassThatFragmentIonScoreIsDoubled)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < (int)binsToSearch.size(); i++)
            {
                std::vector<int> peptideIdsInThisBin = FragmentIndex[binsToSearch[i]];
                
                //get index for minimum monoisotopic allowed
                int lowestPeptideMassIndex = std::isinf(lowestMassPeptideToLookFor) ? 0 : BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, lowestMassPeptideToLookFor, peptideIndex);
                
                // get index for highest mass allowed
                int highestPeptideMassIndex = (int)peptideIdsInThisBin.size() - 1;
                
                if (!std::isinf(highestMassPeptideToLookFor))
                {
                    highestPeptideMassIndex = BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, highestMassPeptideToLookFor,
                                                                               peptideIndex);
                    
                    for (int j = highestPeptideMassIndex; j < (int)peptideIdsInThisBin.size(); j++)
                    {
                        int nextId = peptideIdsInThisBin[j];
                        auto nextPep = peptideIndex[nextId];
                        if (nextPep->getMonoisotopicMass() < highestMassPeptideToLookFor)
                        {
                            highestPeptideMassIndex = j;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                
                // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                {
                    int id = peptideIdsInThisBin[j];
                    scoringTable[id]++;
                    
                    // add possible search results to the hashset of id's
                    if (scoringTable[id] == byteScoreCutoff && massDiffAcceptor->Accepts(scanPrecursorMass,
                                                                                         peptideIndex[id]->getMonoisotopicMass()) >= 0)
                    {
                        idsOfPeptidesPossiblyObserved.push_back(id);
                    }
                }
                
                if (maxMassThatFragmentIonScoreIsDoubled > 0)
                {
                    for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                    {
                        if (j < maxMassThatFragmentIonScoreIsDoubled * FragmentBinsPerDalton)
                        {
                            int id = peptideIdsInThisBin[j];
                            scoringTable[id]++;
                            
                            // add possible search results to the hashset of id's
                            if (scoringTable[id] == byteScoreCutoff && massDiffAcceptor->Accepts(scanPrecursorMass,
                                                                                                 peptideIndex[id]->getMonoisotopicMass()) >= 0)
                            {
                                idsOfPeptidesPossiblyObserved.push_back(id);
                            }
                        }
                    }
                }
            }
        }
    }
}
