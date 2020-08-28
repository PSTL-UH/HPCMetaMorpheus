#include <optional>

#include "NonSpecificEnzymeSearchEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../PrecursorSearchModes/MassDiffAcceptor.h"
#include "../MetaMorpheusEngineResults.h"
#include "../EventArgs/ProgressEventArgs.h"
#include "../ProteinScoringAndFdr/FdrClassifier.h"

#include "bankersrounding.h"
#include "Group.h"

using namespace Chemistry;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::ModernSearch;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
    namespace NonSpecificEnzymeSearch
    {

        //const double NonSpecificEnzymeSearchEngine::WaterMonoisotopicMass = PeriodicTable::GetElement("H")->getPrincipalIsotope()->getAtomicMass() * 2 + PeriodicTable::GetElement("O")->getPrincipalIsotope()->getAtomicMass();
        // Edgar: replacing by the actual value obtained by running
        const double NonSpecificEnzymeSearchEngine::WaterMonoisotopicMass = 18.01056468403;
        
        NonSpecificEnzymeSearchEngine::NonSpecificEnzymeSearchEngine(std::vector<std::vector<PeptideSpectralMatch*>> &globalPsms,
                                                                     std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans,
                                                                     std::vector<PeptideWithSetModifications*> &peptideIndex,
                                                                     std::vector<std::vector<int>> &fragmentIndex,
                                                                     std::vector<std::vector<int>> &precursorIndex,
                                                                     int currentPartition,
                                                                     CommonParameters *commonParams,
                                                                     MassDiffAcceptor *massDiffAcceptor,
                                                                     double maximumMassThatFragmentIonScoreIsDoubled,
                                                                     std::vector<std::string> &nestedIds) :
            ModernSearchEngine(unusedPsms,
                               listOfSortedms2Scans,
                               peptideIndex,
                               fragmentIndex,
                               currentPartition,
                               commonParams,
                               massDiffAcceptor,
                               maximumMassThatFragmentIonScoreIsDoubled,
                               nestedIds),
            PrecursorIndex(precursorIndex),
            MinimumPeptideLength(commonParameters->getDigestionParams()->getMinPeptideLength())
        {
            GlobalCategorySpecificPsms = globalPsms;
            //ModifiedParametersNoComp = commonParameters->CloneWithNewTerminus(, std::make_optional(false));
            std::optional<FragmentationTerminus> t1;
            std::optional<bool> t2;
            ModifiedParametersNoComp = new CommonParameters(commonParams, t1, t2);

#ifdef ORIG
            ProductTypesToSearch = DissociationTypeCollection::ProductsFromDissociationType[commonParameters->getDissociationType()].
                Intersect(
                    TerminusSpecificProductTypes::ProductIonTypesFromSpecifiedTerminus[commonParameters->getDigestionParams()->getFragmentationTerminus()]);
#endif
            std::vector<ProductType> t3 = DissociationTypeCollection::ProductsFromDissociationType[commonParameters->getDissociationType()];
            std::vector<ProductType> t4 = TerminusSpecificProductTypes::ProductIonTypesFromSpecifiedTerminus[commonParameters->getDigestionParams()->getFragmentationTerminus()];
            for ( auto p: t3 ) {
                bool found = false;
                for ( auto q : t4 ) {
                    if ( p == q ) {
                        found = true;
                        break;
                    }
                }
                if ( found ) {
                    ProductTypesToSearch.push_back(p);
                }
            }
        }
        
        MetaMorpheusEngineResults *NonSpecificEnzymeSearchEngine::RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ProgressEventArgs tempVar(oldPercentProgress, "Performing nonspecific search... " +
                                      std::to_string(CurrentPartition) + "/" +
                                      std::to_string(commonParameters->getTotalPartitions()),
                                      const_cast<std::vector<std::string>&>(nestedIds));
            ReportProgress(&tempVar);
            
            unsigned char byteScoreCutoff = static_cast<unsigned char>(commonParameters->getScoreCutoff());
            
#ifdef ORIG
           //ParallelOptions *tempVar2 = new ParallelOptions();
           //tempVar2->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
           // Parallel::ForEach(Partitioner::Create(0, ListOfSortedMs2Scans.size()), tempVar2, [&] (std::any range)  {
#endif
            std::vector<unsigned char> scoringTable(PeptideIndex.size());
            std::unordered_set<int> idsOfPeptidesPossiblyObserved;
                    
#ifdef ORIG
            //for (int i = range::Item1; i < range::Item2; i++)
#endif
            for (int i = 0; i < (int)ListOfSortedMs2Scans.size(); i++)
            {
                // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                //Array::Clear(scoringTable, 0, scoringTable.Length);
                scoringTable.clear();
                idsOfPeptidesPossiblyObserved.clear();
                Ms2ScanWithSpecificMass *scan = ListOfSortedMs2Scans[i];
                
                //get bins to add points to
                std::vector<int> allBinsToSearch = GetBinsToSearch(scan);
                
                //the entire indexed scoring is done here
                for (int j = 0; j < (int)allBinsToSearch.size(); j++)
                {
                    for (auto id = FragmentIndex[allBinsToSearch[j]].begin(); id!=FragmentIndex[allBinsToSearch[j]].end(); id++)    {
                        scoringTable[*id]++;
                    }
                }
                        
                //populate ids of possibly observed with those containing allowed precursor masses
                std::vector<AllowedIntervalWithNotch*> validIntervals = massDiffAcceptor->GetAllowedPrecursorMassIntervalsFromObservedMass(scan->getPrecursorMass()); //.ToList(); //get all valid notches
                for (auto interval : validIntervals)
                {
                    int obsPrecursorFloorMz = static_cast<int>(std::floor(interval->AllowedInterval->getMinimum() * FragmentBinsPerDalton));
                    int obsPrecursorCeilingMz = static_cast<int>(std::ceil(interval->AllowedInterval->getMaximum() * FragmentBinsPerDalton));
                    
                    for (auto pt : ProductTypesToSearch)
                    {
                        int dissociationBinShift = static_cast<int>(BankersRounding::round((WaterMonoisotopicMass -
                                                           DissociationTypeCollection::GetMassShiftFromProductType(pt))
                                                                                           * FragmentBinsPerDalton));
                        int lowestBin = obsPrecursorFloorMz - dissociationBinShift;
                        int highestBin = obsPrecursorCeilingMz - dissociationBinShift;
                        for (int bin = lowestBin; bin <= highestBin; bin++)
                        {
                            if (bin < (int)FragmentIndex.size() && FragmentIndex[bin].size() > 0)
                            {
                                for (auto id = FragmentIndex[bin].begin(); id != FragmentIndex[bin].end(); id++ ) {
                                    idsOfPeptidesPossiblyObserved.insert(*id);
                                }
                            }
                        }
                    }
                    
                    for (int bin = obsPrecursorFloorMz; bin <= obsPrecursorCeilingMz; bin++) //no bin shift, since they're precursor masses
                    {
                        if (bin < (int)PrecursorIndex.size() && PrecursorIndex[bin].size() > 0)
                        {
                            for (auto id = PrecursorIndex[bin].begin(); id != PrecursorIndex[bin].end();  id++) {
                                idsOfPeptidesPossiblyObserved.insert(*id);
                            }
                        }
                    }
                }
                
                // done with initial scoring; refine scores and create PSMs
                if (!idsOfPeptidesPossiblyObserved.empty())
                {
                    int maxInitialScore= scoringTable[0];
                    for ( auto id : idsOfPeptidesPossiblyObserved ) {
                        if ( maxInitialScore < scoringTable[id] ) {
                            maxInitialScore = scoringTable[id];
                        }
                    };
                    maxInitialScore += 1;
                     
                    while (maxInitialScore > commonParameters->getScoreCutoff()) //go through all until we hit the end
                    {
                        maxInitialScore--;
                        for (int id : idsOfPeptidesPossiblyObserved ) {
                            if ( scoringTable[id] != maxInitialScore ) {
                                continue;
                            }
                            
                            PeptideWithSetModifications *peptide = PeptideIndex[id];
                            std::vector<Product*> peptideTheorProducts = peptide->Fragment(
                                commonParameters->getDissociationType(),
                                commonParameters->getDigestionParams()->getFragmentationTerminus());
                            
                            std::tuple<int, PeptideWithSetModifications*> notchAndUpdatedPeptide = Accepts(peptideTheorProducts,
                                                              scan->getPrecursorMass(),
                                                              peptide,
                                                              commonParameters->getDigestionParams()->getFragmentationTerminus(),
                                                              massDiffAcceptor);
                            int notch = std::get<0>(notchAndUpdatedPeptide);
                            if (notch >= 0)
                            {
                                peptide = std::get<1>(notchAndUpdatedPeptide);
                                peptideTheorProducts = peptide->Fragment(commonParameters->getDissociationType(),
                                                                         FragmentationTerminus::Both);//.ToList();
                                std::vector<MatchedFragmentIon*> matchedIons = MatchFragmentIons(scan, peptideTheorProducts,
                                                                                                 ModifiedParametersNoComp);
                                
                                double thisScore = CalculatePeptideScore(scan->getTheScan(), matchedIons,
                                                                         MaxMassThatFragmentIonScoreIsDoubled);
                                if (thisScore > commonParameters->getScoreCutoff())
                                {
                                    auto tmp  = peptide->getCleavageSpecificityForFdrCategory();
                                    std::vector<PeptideSpectralMatch*> localPeptideSpectralMatches =
                                        GlobalCategorySpecificPsms[static_cast<int>(FdrClassifier::GetCleavageSpecificityCategory(&tmp))];
                                    if (localPeptideSpectralMatches[i] == nullptr)
                                    {
                                        localPeptideSpectralMatches[i] = new PeptideSpectralMatch(peptide, notch, thisScore,
                                                                                                  i, scan,
                                                                                                  commonParameters->getDigestionParams(),
                                                                                                  matchedIons);
                                    }
                                    else
                                    {
                                        localPeptideSpectralMatches[i]->AddOrReplace(peptide, thisScore, notch,
                                                                                     commonParameters->getReportAllAmbiguity(),
                                                                                     matchedIons);
                                    }
                                }
                            }
                        }
                    }
                }
                // report search progress
                progress++;
                int percentProgress = static_cast<int>((progress / ListOfSortedMs2Scans.size()) * 100);
                
                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ProgressEventArgs tempVar3(percentProgress, "Performing nonspecific search... " +
                                               std::to_string(CurrentPartition) + "/" +
                                               std::to_string(commonParameters->getTotalPartitions()),
                                               const_cast<std::vector<std::string>&>(nestedIds));
                    ReportProgress(&tempVar3);
                }
            }
            // });
            
            return new MetaMorpheusEngineResults(this);
        }
        
        std::tuple<int, PeptideWithSetModifications*> NonSpecificEnzymeSearchEngine::Accepts(std::vector<Product*> &fragments,
                                                                                             double scanPrecursorMass,
                                                                                             PeptideWithSetModifications *peptide,
                                                                                             FragmentationTerminus fragmentationTerminus,
                                                                                             MassDiffAcceptor *searchMode)
        {
            //all masses in N and CTerminalMasses are b-ion masses, which are one water away from a full peptide
            int localminPeptideLength = commonParameters->getDigestionParams()->getMinPeptideLength();
            
            for (int i = localminPeptideLength - 1; i < (int)fragments.size(); i++) //minus one start, because fragment 1 is at index 0
            {
                Product *fragment = fragments[i];
                double theoMass = fragment->NeutralMass -
                    DissociationTypeCollection::GetMassShiftFromProductType(fragment->productType) + WaterMonoisotopicMass;
                int notch = searchMode->Accepts(scanPrecursorMass, theoMass);
                if (notch >= 0)
                {
                    PeptideWithSetModifications *updatedPwsm = nullptr;
                    if ( fragmentationTerminus == FragmentationTerminus::N)
                    {
                        //-1 for one based index
                        int endResidue = peptide->getOneBasedStartResidueInProtein() + fragment->TerminusFragment->FragmentNumber - 1;
                        std::unordered_map<int, Modification*> updatedMods;
                        for (auto mod : peptide->getAllModsOneIsNterminus())
                        {
                            //check if we cleaved it off, +1 for N-terminus being mod 1 and first residue being mod 2,
                            //+1 again for the -1 on end residue for one based index, +1 (again) for the one-based start residue
                            if (mod.first < endResidue - peptide->getOneBasedStartResidueInProtein() + 3) 
                            {
                                updatedMods.emplace(mod.first, mod.second);
                            }
                        }
                        updatedPwsm = new PeptideWithSetModifications(peptide->getProtein(), peptide->getDigestionParams(),
                                                                      peptide->getOneBasedStartResidueInProtein(), endResidue,
                                                                      CleavageSpecificity::Unknown, "", 0, updatedMods, 0);
                    }
                    else
                    {
                        //plus one for one based index
                        int startResidue = peptide->getOneBasedEndResidueInProtein() - fragment->TerminusFragment->FragmentNumber + 1; 
                        std::unordered_map<int, Modification*> updatedMods; //updateMods
                        int indexShift = startResidue - peptide->getOneBasedStartResidueInProtein();
                        for (auto mod : peptide->getAllModsOneIsNterminus())
                        {
                            if (mod.first > indexShift + 1) //check if we cleaved it off, +1 for N-terminus being mod 1 and first residue being 2
                            {
                                int key = mod.first - indexShift;
                                updatedMods.emplace(key, mod.second);
                            }
                        }
                        updatedPwsm = new PeptideWithSetModifications(peptide->getProtein(), peptide->getDigestionParams(),
                                                                      startResidue, peptide->getOneBasedEndResidueInProtein(),
                                                                      CleavageSpecificity::Unknown, "", 0, updatedMods, 0);
                    }
                    
                    delete updatedPwsm;
                    return std::tuple<int, PeptideWithSetModifications*>(notch, updatedPwsm);
                }
                else if (theoMass > scanPrecursorMass)
                {
                    break;
                }
            }
            //if the theoretical and experimental have the same mass
            if ( (int)fragments.size() > localminPeptideLength)
            {
                double totalMass = peptide->getMonoisotopicMass(); // + Constants.ProtonMass;
                int notch = searchMode->Accepts(scanPrecursorMass, totalMass);
                if (notch >= 0)
                {
                    auto tempv = peptide->getAllModsOneIsNterminus();
                    //need to update so that the cleavage specificity is recorded
                    PeptideWithSetModifications *updatedPwsm = new PeptideWithSetModifications(peptide->getProtein(),
                                                                                               peptide->getDigestionParams(),
                                                                                               peptide->getOneBasedStartResidueInProtein(),
                                                                                               peptide->getOneBasedEndResidueInProtein(),
                                                                                               CleavageSpecificity::Unknown, "", 0,
                                                                                               tempv,
                                                                                               peptide->NumFixedMods);
                    
                    return std::make_tuple(notch, updatedPwsm);
                }
            }
            return std::make_tuple(-1, nullptr);
        }
        
        std::vector<PeptideSpectralMatch*> NonSpecificEnzymeSearchEngine::ResolveFdrCategorySpecificPsms( std::vector<std::vector<PeptideSpectralMatch*>> &AllPsms, int numNotches, const std::string &taskId, CommonParameters *commonParameters )
        {
            //update all psms with peptide info
#ifdef ORIG
            AllPsms.ToList()->Where([&] (std::any psmArray) {
                    return psmArray != nullptr;
                }).ToList()->ForEach([&] (std::any psmArray) {
                        psmArray::Where([&] (std::any psm) {
                                return psm != nullptr;
                            }).ToList()->ForEach([&] (std::any psm){
                                    psm::ResolveAllAmbiguities();
				});
                    });
#endif
            std::vector<std::vector<PeptideSpectralMatch*>> tmpPsms;
            int current=0;
            for ( auto psmArray : AllPsms ) {
                auto newvec = new std::vector<PeptideSpectralMatch*> ();

                if ( !psmArray.empty() ) {
                    tmpPsms.push_back(*newvec);
                    for ( auto psm : psmArray ) {
                        if ( psm != nullptr ) {
                            psm->ResolveAllAmbiguities();
                            tmpPsms[current].push_back(psm);                            
                        }
                    }
                    current++;
                }
            }
                            
            // Now clear AllPsms and copy everything over.
            for ( auto psmArray: AllPsms ) {
                if ( !psmArray.empty() ) {
                    psmArray.clear();
                }
            }
            AllPsms.clear();

            current=0;
            for ( auto psmArray : tmpPsms ) {
                auto newvec = new std::vector<PeptideSpectralMatch*> ();
                AllPsms.push_back(*newvec);
                for ( auto psm: psmArray ) {
                    AllPsms[current].push_back(psm);                    
                }
                current++;
            }
            // Clear tmpPsms
            for ( auto psmArray: tmpPsms ) {
                if ( !psmArray.empty() ) {
                    psmArray.clear();
                }
            }
            tmpPsms.clear();
            
            for (auto psmsArray : AllPsms)
            {
                if (psmsArray.size() > 0)
                {
#ifdef ORIG
                    std::vector<PeptideSpectralMatch*> cleanedPsmsArray = psmsArray.Where([&] (std::any b)    {
                            return b != nullptr;
                        }).OrderByDescending([&] (std::any b)   {
                                b::Score;
                            }).ThenBy([&] (std::any b) {
                                    b::PeptideMonisotopicMass.HasValue ?
                                    std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) :
                                    std::numeric_limits<double>::max();
                                }).GroupBy([&] (std::any b)  {
                                        (b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass);
                                    })->Select([&] (std::any b) {
                                            b::First();
					}).ToList();
#endif
                    std::vector<PeptideSpectralMatch*> tmpPsmsArray;
                    for ( auto b: psmsArray ) {
                        if ( b != nullptr ) {
                            tmpPsmsArray.push_back(b);
                        }
                    }
                    std::sort(tmpPsmsArray.begin(), tmpPsmsArray.end(), [&] (PeptideSpectralMatch* l , PeptideSpectralMatch* r) {
                            if ( l->getScore() > r->getScore() ) return true;
                            if ( l->getScore() < r->getScore() ) return false;

                            double lval = std::numeric_limits<double>::max();
                            double rval = std::numeric_limits<double>::max();
                            if ( l->getPeptideMonisotopicMass().has_value() ) {
                                lval = std::abs(l->getScanPrecursorMass() - l->getPeptideMonisotopicMass().value());
                            }
                            if ( r->getPeptideMonisotopicMass().has_value() ) {
                                rval = std::abs(r->getScanPrecursorMass() - r->getPeptideMonisotopicMass().value());
                            }
                            if ( lval < rval ) return true;
                            return false;
                        });

                    std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f1 = [&](PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                        if ( l->getFullFilePath() < r->getFullFilePath() ) return true;
                        if ( l->getFullFilePath() > r->getFullFilePath() ) return true;

                        if ( l->getScanNumber() < r->getScanNumber() ) return true;
                        if ( l->getScanNumber() > r->getScanNumber() ) return true;

                        if ( l->getPeptideMonisotopicMass() < r->getPeptideMonisotopicMass() ) return true;
                        return true;
                    } ;
                    std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f2 = [&](PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                        return l->getFullFilePath() ==  r->getFullFilePath() &&
                               l->getScanNumber()   == r->getScanNumber() &&
                               l->getPeptideMonisotopicMass() == r->getPeptideMonisotopicMass();
                    };
                    std::vector<std::vector<PeptideSpectralMatch*>> cond = Group::GroupBy ( tmpPsmsArray, f1, f2);
                    std::vector<PeptideSpectralMatch*> cleanedPsmsArray;
                    for ( auto c: cond ) {
                        cleanedPsmsArray.push_back(c.front() );
                    }

                    std::vector<std::string> vs1 = {taskId};
                    FdrAnalysisEngine tempVar(cleanedPsmsArray, numNotches, commonParameters, vs1);
                    (&tempVar)->Run();
                    
                    for (int i = 0; i < (int)psmsArray.size(); i++)
                    {
                        if (psmsArray[i] != nullptr)
                        {
                            if (psmsArray[i]->getFdrInfo() == nullptr) //if it was grouped in the cleanedPsmsArray
                            {
                                psmsArray[i] = nullptr;
                            }
                        }
                    }
                }
            }
            
            std::vector<int> ranking(AllPsms.size()); //high int is good ranking
            std::vector<int> indexesOfInterest;
            for (int i = 0; i < (int)ranking.size(); i++)
            {
                if (AllPsms[i].size() > 0)
                {
#ifdef ORIG
                    ranking[i] = AllPsms[i].Where([&] (std::any x)   {
                            return x != nullptr;
                        })->Count([&] (std::any x)  {
                                return x::FdrInfo::QValue <= 0.01;
                            }); //set ranking as number of psms above 1% FDR
#endif
                    ranking[i]=0;
                    for ( auto x: AllPsms[i] ) {
                        if ( x != nullptr && x->getFdrInfo()->getQValue() <= 0.01 ) {
                            ranking[i]++;
                        }
                    }
                    indexesOfInterest.push_back(i);
                }
            }
            
            //get the index of the category with the highest ranking
            int majorCategoryIndex = indexesOfInterest[0];
            for (int i = 1; i < (int)indexesOfInterest.size(); i++)
            {
                int currentCategoryIndex = indexesOfInterest[i];
                if (ranking[currentCategoryIndex] > ranking[majorCategoryIndex])
                {
                    majorCategoryIndex = currentCategoryIndex;
                }
            }
            
            //update other category q-values
            //There's a chance of weird categories getting a random decoy before a random target,
            //but we don't want to give that target a q value of zero.
            //We can't just take the q of the first decoy, because if the target wasn't random (score = 40),
            //but there are no other targets before the decoy (score = 5), then we're incorrectly dinging the target
            //The current solution is such that if a minor category has a lower q value than it's corresponding
            //score in the major category, then its q-value is changed to what it would be in the major category
#ifdef ORIG
            std::vector<PeptideSpectralMatch*> majorCategoryPsms = AllPsms[majorCategoryIndex].Where([&] (std::any x)  {
                    return x != nullptr;
                }).OrderByDescending([&] (std::any x)   {
                        x::Score;
                    }).ToList(); //get sorted major category
#endif
            std::vector<PeptideSpectralMatch*> majorCategoryPsms;
            for ( auto x : AllPsms[majorCategoryIndex] ) {
                if ( x != nullptr ) {
                    majorCategoryPsms.push_back(x);
                }
            }
            std::sort(majorCategoryPsms.begin(), majorCategoryPsms.end(), [&] (PeptideSpectralMatch* l , PeptideSpectralMatch* r) {
                    return l->getScore() > r->getScore();
                });
                      
            for (int i = 0; i < (int)indexesOfInterest.size(); i++)
            {
                int minorCategoryIndex = indexesOfInterest[i];
                if (minorCategoryIndex != majorCategoryIndex)
                {
#ifdef ORIG
                    std::vector<PeptideSpectralMatch*> minorCategoryPsms = AllPsms[minorCategoryIndex].Where([&] (std::any x){
                            return x != nullptr;
                        }).OrderByDescending([&] (std::any x){
                                x::Score;
                            }).ToList(); //get sorted minor category
#endif
                    std::vector<PeptideSpectralMatch*> minorCategoryPsms;
                    for ( auto x : AllPsms[minorCategoryIndex] ) {
                        if ( x != nullptr ) {
                            minorCategoryPsms.push_back(x);
                        }
                    }
                    std::sort(minorCategoryPsms.begin(), minorCategoryPsms.end(), [&] (PeptideSpectralMatch* l , PeptideSpectralMatch* r) {
                            return l->getScore() > r->getScore();
                        });
                    
                    int minorPsmIndex = 0;
                    int majorPsmIndex = 0;
                    while (minorPsmIndex < (int)minorCategoryPsms.size() &&
                           majorPsmIndex < (int)majorCategoryPsms.size()) //while in the lists
                    {
                        PeptideSpectralMatch *majorPsm = majorCategoryPsms[majorPsmIndex];
                        PeptideSpectralMatch *minorPsm = minorCategoryPsms[minorPsmIndex];
                        //major needs to be a lower score than the minor
                        if (majorPsm->getScore() > minorPsm->getScore())
                        {
                            majorPsmIndex++;
                        }
                        else
                        {
                            if (majorPsm->getFdrInfo()->getQValue() > minorPsm->getFdrInfo()->getQValue())
                            {
                                minorPsm->getFdrInfo()->setQValue(majorPsm->getFdrInfo()->getQValue());
                            }
                            minorPsmIndex++;
                        }
                    }
                    //wrap up if we hit the end of the major category
                    while (minorPsmIndex < (int)minorCategoryPsms.size())
                    {
                        PeptideSpectralMatch *majorPsm = majorCategoryPsms[majorPsmIndex - 1]; //-1 because it's out of index right now
                        PeptideSpectralMatch *minorPsm = minorCategoryPsms[minorPsmIndex];
                        if (majorPsm->getFdrInfo()->getQValue() > minorPsm->getFdrInfo()->getQValue())
                        {
                            minorPsm->getFdrInfo()->setQValue(majorPsm->getFdrInfo()->getQValue());
                        }
                        minorPsmIndex++;
                    }
                }
            }
            
            int numTotalSpectraWithPrecursors = AllPsms[indexesOfInterest[0]].size();
            std::vector<PeptideSpectralMatch*> bestPsmsList;
            for (int i = 0; i < numTotalSpectraWithPrecursors; i++)
            {
                PeptideSpectralMatch *bestPsm = nullptr;
                double lowestQ = std::numeric_limits<double>::max();
                int bestIndex = -1;
                for (auto index : indexesOfInterest) //foreach category
                {
                    PeptideSpectralMatch *currentPsm = AllPsms[index][i];
                    if (currentPsm != nullptr)
                    {
                        double currentQValue = currentPsm->getFdrInfo()->getQValue();
                        if (currentQValue < lowestQ || (currentQValue == lowestQ && currentPsm->getScore() > bestPsm->getScore()))
                        {
                            if (bestIndex != -1)
                            {
                                //remove the old one so we don't use it for fdr later
                                AllPsms[bestIndex][i] = nullptr;
                            }
                            bestPsm = currentPsm;
                            lowestQ = currentQValue;
                            bestIndex = index;
                        }
                        else //remove the old one so we don't use it for fdr later
                        {
                            AllPsms[index][i] = nullptr;
                        }
                    }
                }
                if (bestPsm != nullptr)
                {
                    bestPsmsList.push_back(bestPsm);
                }
            }
            
            //It's probable that psms from some categories were removed by psms from other categories.
            //however, the fdr is still affected by their presence, since it was calculated before their removal.
            for (auto psmsArray : AllPsms)
            {
                if (psmsArray.size() > 0)
                {
#ifdef ORIG
                    std::vector<PeptideSpectralMatch*> cleanedPsmsArray = psmsArray.Where([&] (std::any b)  {
                            return b != nullptr;
                        }).OrderByDescending([&] (std::any b) {
                                b::Score;
                            }).ThenBy([&] (std::any b)	{
                                    b::PeptideMonisotopicMass.HasValue ?
                                    std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) :
                                    std::numeric_limits<double>::max();
                                }).ToList();
#endif
                    std::vector<PeptideSpectralMatch*> cleanedPsmsArray;
                    for ( auto x : psmsArray ) {
                        if ( x != nullptr ) {
                            cleanedPsmsArray.push_back(x);
                        }
                    }
                    std::sort(cleanedPsmsArray.begin(), cleanedPsmsArray.end(), [&] (PeptideSpectralMatch* l , PeptideSpectralMatch* r) {
                            if ( l->getScore() > r->getScore() ) return true;
                            if ( l->getScore() < r->getScore() ) return false;

                            double lval = std::numeric_limits<double>::max();
                            double rval = std::numeric_limits<double>::max();
                            if ( l->getPeptideMonisotopicMass().has_value() ) {
                                lval = std::abs(l->getScanPrecursorMass() - l->getPeptideMonisotopicMass().value());
                            }
                            if ( r->getPeptideMonisotopicMass().has_value() ) {
                                rval = std::abs(r->getScanPrecursorMass() - r->getPeptideMonisotopicMass().value());
                            }
                            if ( lval < rval ) return true;
                            return false;
                        });
                    
                    
                    std::vector<std::string> vs1 = {taskId};
                    FdrAnalysisEngine tempVar2(cleanedPsmsArray, numNotches, commonParameters, vs1);
                    (&tempVar2)->Run();
                }
            }
           
#ifdef ORIG
            return bestPsmsList.OrderBy([&] (std::any b) {
                    b::FdrInfo::QValue;
                }).ThenByDescending([&] (std::any b) {
                        b::Score;
                    }).ToList();
#endif
            std::sort(bestPsmsList.begin(), bestPsmsList.end(), [&] (PeptideSpectralMatch* l , PeptideSpectralMatch* r) {
                    if ( l->getFdrInfo()->getQValue() < r->getFdrInfo()->getQValue() ) return true;
                    if ( l->getFdrInfo()->getQValue() > r->getFdrInfo()->getQValue() ) return true;
                    
                    if ( l->getScore() > r->getScore() ) return true;
                    return false;
                });
            return bestPsmsList;
        }
    }
}
