#include "ProteinParsimonyEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "ProteinGroup.h"

#include "Group.h"

using namespace EngineLayer;
using namespace EngineLayer::ProteinParsimony;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

    ProteinParsimonyEngine::ProteinParsimonyEngine(std::vector<PeptideSpectralMatch*> &allPsms, bool modPeptidesAreDifferent,
                                                   CommonParameters *commonParameters, std::vector<std::string> &nestedIds,
                                                   int verbosityLevel ) :
        MetaMorpheusEngine(commonParameters, nestedIds, verbosityLevel ),
        _fdrFilteredPeptides(std::unordered_set<PeptideWithSetModifications*>()), _allPsms(allPsms),
        _treatModPeptidesAsDifferentPeptides(modPeptidesAreDifferent)
    {
        
        // parsimony will only use non-ambiguous, high-confidence PSMs
        // KEEP decoys and contaminants for use in parsimony!
        if (modPeptidesAreDifferent)
        {
#ifdef ORIG
            _fdrFilteredPsms = allPsms.Where([&] (std::any p) {
                    return p::FullSequence != nullptr && p::FdrInfo::QValue <= FdrCutoffForParsimony &&
                    p::FdrInfo::QValueNotch <= FdrCutoffForParsimony;
                }).ToList();
#endif
            for ( auto p: allPsms ) {
                if (p->getFullSequence().length() != 0                      &&
                    p->getFdrInfo()->getQValue() <= FdrCutoffForParsimony &&
                    p->getFdrInfo()->getQValueNotch() <= FdrCutoffForParsimony ) {
                    _fdrFilteredPsms.push_back(p);
                }
            }
        }
        else
        {
#ifdef ORIG
            _fdrFilteredPsms = allPsms.Where([&] (std::any p) {
                    return p::BaseSequence != nullptr && p::FdrInfo::QValue <= FdrCutoffForParsimony &&
                    p::FdrInfo::QValueNotch <= FdrCutoffForParsimony;
                }).ToList();
#endif
            for ( auto p: allPsms ) {
                if (p->getBaseSequence().length() != 0                    &&
                    p->getFdrInfo()->getQValue() <= FdrCutoffForParsimony &&
                    p->getFdrInfo()->getQValueNotch() <= FdrCutoffForParsimony ) {
                    _fdrFilteredPsms.push_back(p);
                }
            }
        }
        
        // if PSM is a decoy, add only decoy sequences; same for contaminants
        // peptides to use in parsimony = peptides observed in high-confidence PSMs
        
        for (auto psm : _fdrFilteredPsms)
        {
            if (psm->getIsDecoy())
            {
#ifdef ORIG
                for (auto peptide : psm->BestMatchingPeptides->Select([&] (std::any p)    {
                            p::Peptide;
                        }).Where([&] (std::any p)  {
                                p::Protein::IsDecoy;
                            }))
                {
                    _fdrFilteredPeptides.Add(peptide);
                }
#endif
                for ( auto p : psm->getBestMatchingPeptides() ) {
                    auto peptide = std::get<1>(p);
                    if ( peptide->getProtein()->getIsDecoy() ) {
                        _fdrFilteredPeptides.insert(peptide);
                    }
                }
            }
            else if (psm->getIsContaminant() )
            {
                for ( auto p : psm->getBestMatchingPeptides() ) {
                    auto peptide = std::get<1>(p);
                    if ( peptide->getProtein()->getIsContaminant() ) {                        
                        _fdrFilteredPeptides.insert(peptide);
                    }
                }
            }
            else // PSM is target
            {
                for ( auto p : psm->getBestMatchingPeptides() ) {
                    auto peptide = std::get<1>(p);
                    if (!peptide->getProtein()->getIsDecoy() &&
                        !peptide->getProtein()->getIsContaminant() ) {
                        _fdrFilteredPeptides.insert(peptide);
                    }
                }
            }
        }
        // we're storing all PSMs (not just FDR-filtered ones) here because we will remove some protein associations 
        // from low-confidence PSMs if they can be explained by a parsimonious protein
        // allready done at the top of the constructor.
        // _allPsms = allPsms;
    }
    
    MetaMorpheusEngineResults* ProteinParsimonyEngine::RunSpecific()
    {
        auto myAnalysisResults = new ProteinParsimonyResults(this);
        myAnalysisResults->setProteinGroups(RunProteinParsimonyEngine());
        
        return myAnalysisResults;
    }
    
    
    std::vector<ProteinGroup*> ProteinParsimonyEngine::RunProteinParsimonyEngine()
    {
        // parsimonious list of proteins built by this protein parsimony engine
        // EDGAR: changing from unordered_set to vector, because that's what is really required later
        std::vector<Protein*> parsimoniousProteinList;
        
        // list of peptides that can only be digestion products of one protein in the proteome (considering different
        // protease digestion rules)
        std::unordered_set<PeptideWithSetModifications*> uniquePeptides;
        
        // if there are no peptides observed, there are no proteins; return an empty list of protein groups
        if (_fdrFilteredPeptides.size()  == 0) {
            return std::vector<ProteinGroup*>();
        }
        
        // Parsimony stage 0: create peptide-protein associations if needed because the user wants a
        // modification-agnostic parsimony
        if (!_treatModPeptidesAsDifferentPeptides) {
#ifdef ORIG
            //for (auto protease : _fdrFilteredPsms::GroupBy([&] (std::any p) {
            //            p::DigestionParams::Protease;
            //        })) {
#endif
            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f1 = [&](PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                return l->digestionParams->getProtease() < r->digestionParams->getProtease(); } ;
            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f2 = [&](PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                return l->digestionParams->getProtease() != r->digestionParams->getProtease(); } ;
            std::vector<std::vector<PeptideSpectralMatch*>> tmpPsms = Group::GroupBy ( _fdrFilteredPsms, f1, f2);
            
            for ( auto protease : tmpPsms ) {
                std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>> sequenceWithPsms;

                // for each protease, match the base sequence of each peptide to its PSMs
                for (PeptideSpectralMatch *psm : protease) {
                    std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>>::const_iterator sequenceWithPsms_iterator = sequenceWithPsms.find(psm->getBaseSequence());
                    if (sequenceWithPsms_iterator != sequenceWithPsms.end()) {
                        sequenceWithPsms[psm->getBaseSequence()].push_back(psm);
                    }
                    else {
                        sequenceWithPsms.emplace(psm->getBaseSequence(), std::vector<PeptideSpectralMatch*>{psm} );
                    }
                }

                // create new peptide-protein associations
                for (auto baseSequence : sequenceWithPsms) {
#ifdef ORIG
                    auto peptidesWithNotchInfo = baseSequence.Value->SelectMany([&] (std::any p) {
                            p::BestMatchingPeptides;
                        })->Distinct()->ToList();
#endif
                    std::vector<std::tuple<int, PeptideWithSetModifications*>> peptidesWithNotchInfo;
                    for ( auto p : std::get<1>(baseSequence) ) {
                        for ( auto b  : p->getBestMatchingPeptides() ) {
                            bool  found = false;
                            for ( auto z : peptidesWithNotchInfo ) {
                                if ( std::get<0>(b) == std::get<0>(z) &&
                                     std::get<1>(b)->Equals(std::get<1>(z) )) {
                                    found = true;
                                    break;
                                }
                            }
                            if (!found ) {
                                peptidesWithNotchInfo.push_back(b);
                            }
                        }
                    }
                    // if the base seq has >1 PeptideWithSetMods object and has >0 mods, it might need to be matched
                    //to new proteins
#ifdef ORIG
                    //if (peptidesWithNotchInfo.size() > 1 && peptidesWithNotchInfo.Any([&] (std::any p) {
                    //            return p::Peptide::NumMods > 0;
                    //        })) {
#endif
                    bool cond=false;
                    for ( auto p: peptidesWithNotchInfo ) {
                        if ( std::get<1>(p)->getNumMods() > 0 ) {
                            cond = true;
                            break;
                        }
                    }
                    if (peptidesWithNotchInfo.size() > 1 && cond ) {
                        // list of proteins along with start/end residue in protein and the # missed cleavages
                        // this is needed to create new PeptideWithSetModification objects
                        auto peptideInProteinInfo = std::vector<std::tuple<Protein*, DigestionParams*, int, int, int, int>>();
                        for (auto peptide : peptidesWithNotchInfo) {
                            peptideInProteinInfo.push_back(
                                std::make_tuple(std::get<1>(peptide)->getProtein(),
                                                std::get<1>(peptide)->getDigestionParams(),
                                                std::get<1>(peptide)->getOneBasedStartResidueInProtein(),
                                                std::get<1>(peptide)->getOneBasedEndResidueInProtein(),
                                                std::get<1>(peptide)->getMissedCleavages(),
                                                std::get<0>(peptide) ));
                        }
                        
                        // add the protein associations to the PSM
                        for (PeptideSpectralMatch *psm : std::get<1>(baseSequence)) {
                            for (auto proteinInfo : peptideInProteinInfo) {
                                auto originalPep = std::get<1>(psm->getBestMatchingPeptides().front());
                                auto telem  = originalPep->getAllModsOneIsNterminus();
                                auto pep = new PeptideWithSetModifications(std::get<0>(proteinInfo),
                                                                           std::get<1>(proteinInfo),
                                                                           std::get<2>(proteinInfo),
                                                                           std::get<3>(proteinInfo),
                                                                           originalPep->getCleavageSpecificityForFdrCategory(),
                                                                           originalPep->getPeptideDescription(),
                                                                           std::get<4>(proteinInfo),
                                                                           telem,
                                                                           originalPep->NumFixedMods);
                                _fdrFilteredPeptides.insert(pep);
                                psm->AddProteinMatch(std::make_tuple(std::get<5>(proteinInfo), pep));

                                //C# TO C++ CONVERTER TODO TASK: A 'delete pep' statement was not added since pep was
                                //passed to a method or constructor. Handle memory management manually.
                            }
                        }
                    }
                }
            }
        }

        // Parsimony stage 1: add proteins with unique peptides (for each protease)
#ifdef ORIG
        auto peptidesGroupedByProtease = _fdrFilteredPeptides::GroupBy([&] (std::any p) {
                p::DigestionParams::Protease;
            });
#endif
        std::vector<PeptideWithSetModifications *> tmpFilteredPeptides;
        for ( auto p = _fdrFilteredPeptides.begin(); p != _fdrFilteredPeptides.end(); p++ ) {
            tmpFilteredPeptides.push_back(*p);
        }
        std::function<bool(PeptideWithSetModifications*,PeptideWithSetModifications*)> f3 = [&](PeptideWithSetModifications *l, PeptideWithSetModifications *r) {
            return l->getDigestionParams()->getProtease() < r->getDigestionParams()->getProtease(); } ;
        std::function<bool(PeptideWithSetModifications*,PeptideWithSetModifications*)> f4 = [&](PeptideWithSetModifications *l, PeptideWithSetModifications *r) {
            return l->getDigestionParams()->getProtease() != r->getDigestionParams()->getProtease(); } ;
        std::vector<std::vector<PeptideWithSetModifications*>> peptidesGroupedByProtease = Group::GroupBy ( tmpFilteredPeptides, f3, f4);
        
        for (auto peptidesForThisProtease : peptidesGroupedByProtease) {
            std::unordered_map<std::string, std::vector<Protein*>> peptideSequenceToProteinsForThisProtease;
            std::unordered_map<std::string, std::vector<PeptideWithSetModifications*>> sequenceToPwsm;

            for (auto peptide : peptidesForThisProtease) {
                std::string sequence = peptide->getBaseSequence();
                if (_treatModPeptidesAsDifferentPeptides) {
                    //these and next set to full sequence but might be base sequence. treat modified as unique makes
                    //sense to use full
                    sequence = peptide->getFullSequence();
                }

                std::unordered_map<std::string, std::vector<Protein*>>::const_iterator peptideSequenceToProteinsForThisProtease_iterator = peptideSequenceToProteinsForThisProtease.find(sequence);
                if (peptideSequenceToProteinsForThisProtease_iterator != peptideSequenceToProteinsForThisProtease.end()) {
                    peptideSequenceToProteinsForThisProtease[sequence].push_back(peptide->getProtein());
                }
                else {
                    peptideSequenceToProteinsForThisProtease.emplace(sequence, std::vector<Protein*> {peptide->getProtein()});
                }

                std::unordered_map<std::string, std::vector<PeptideWithSetModifications*>>::const_iterator sequenceToPwsm_iterator = sequenceToPwsm.find(sequence);
                if (sequenceToPwsm_iterator != sequenceToPwsm.end()) {
                    sequenceToPwsm[sequence].push_back(peptide);
                }
                else {
                    sequenceToPwsm.emplace(sequence, std::vector<PeptideWithSetModifications*> {peptide} );
                }
            }

#ifdef ORIG
            //for (auto uniquePeptide : peptideSequenceToProteinsForThisProtease.Where([&] (std::any p) {
            //            return p->Value->Count == 1;
            //        })) {
#endif
            for (auto uniquePeptide : peptideSequenceToProteinsForThisProtease ) {
                if ( std::get<1>(uniquePeptide).size() != 1 )  {
                    continue;
                }
                // add the protein with the unique peptide to the parsimonious protein list
                Protein *proteinWithUniquePeptideSequence = std::get<1>(uniquePeptide).front();
                parsimoniousProteinList.push_back(proteinWithUniquePeptideSequence);

                // add the unique peptide to the list of unique peptides
                PeptideWithSetModifications *uniquePwsm = sequenceToPwsm[std::get<0>(uniquePeptide)].front();
                uniquePeptides.insert(uniquePwsm);
            }
        }
        
        // Parsimony stage 2: build the peptide-protein matching structure for the parsimony greedy algorithm
        // and remove all peptides observed by proteins with unique peptides
        std::unordered_map<ParsimonySequence*, std::vector<Protein*>> peptideSequenceToProteins = std::unordered_map<ParsimonySequence*, std::vector<Protein*>>();

        // this dictionary associates proteins w/ all peptide sequences (list will NOT shrink over time)
        // this is used in case of greedy algorithm ties to figure out which protein has more total peptides observed
        std::unordered_map<Protein*, std::unordered_set<ParsimonySequence*>> proteinToPepSeqMatch;
        for (auto peptide : _fdrFilteredPeptides) {
            ParsimonySequence *sequence = new ParsimonySequence(peptide, _treatModPeptidesAsDifferentPeptides);

            std::unordered_map<ParsimonySequence*, std::vector<Protein*>>::const_iterator peptideSequenceToProteins_iterator = peptideSequenceToProteins.find(sequence);
            if (peptideSequenceToProteins_iterator != peptideSequenceToProteins.end()) {
                peptideSequenceToProteins[sequence].push_back(peptide->getProtein());
            }
            else {
                peptideSequenceToProteins.emplace(sequence,  std::vector<Protein*> {peptide->getProtein()} );
            }

            std::unordered_map<Protein*, std::unordered_set<ParsimonySequence*>>::const_iterator proteinToPepSeqMatch_iterator = proteinToPepSeqMatch.find(peptide->getProtein());
            if (proteinToPepSeqMatch_iterator != proteinToPepSeqMatch.end()) {
                proteinToPepSeqMatch[peptide->getProtein()].insert(sequence);
            }
            else {
                proteinToPepSeqMatch.emplace(peptide->getProtein(), std::unordered_set<ParsimonySequence*> {sequence});
            }

            //C# TO C++ CONVERTER TODO TASK: A 'delete sequence' statement was not added since sequence was passed
            //to a method or constructor. Handle memory management manually.
        }

        // remove the peptides observed by proteins with unique peptides
        std::unordered_set<ParsimonySequence*> toRemove;
        for (auto seq : peptideSequenceToProteins) {
#ifdef ORIG
            bool observedAlready = seq.Value->Any([&] (std::any p) {
                    std::find(parsimoniousProteinList.begin(), parsimoniousProteinList.end(), p) != parsimoniousProteinList.end();
                });
#endif
            bool observedAlready = false;
            for ( auto p: std::get<1>(seq) ) {
                if ( std::find(parsimoniousProteinList.begin(), parsimoniousProteinList.end(), p) !=
                     parsimoniousProteinList.end() ) {
                    observedAlready = true;
                    break;
                }
            }

            if (observedAlready) {
                toRemove.insert(std::get<0>(seq) );
            }
        }
        for (auto sequence : toRemove) {
            peptideSequenceToProteins.erase(sequence);
        }

        if (!peptideSequenceToProteins.empty()) {
            // Parsimony stage 3: greedy algorithm

            // dictionary with proteins as keys and list of associated peptide sequences as the values.
            // this data structure makes parsimony easier because the algorithm can look up a protein's peptides
            // to remove them from the list of available peptides. this list will shrink as the algorithm progresses
            std::unordered_map<Protein*, std::unordered_set<std::string>> algDictionary;
            std::unordered_map<Protein*, std::unordered_set<ParsimonySequence*>> algDictionaryProtease;
            for (auto kvp : peptideSequenceToProteins) {
                for (auto protein : std::get<1>(kvp) ) {
                    std::unordered_map<Protein*, std::unordered_set<ParsimonySequence*>>::const_iterator algDictionaryProtease_iterator = algDictionaryProtease.find(protein);
                    if (algDictionaryProtease_iterator != algDictionaryProtease.end()) {
                        algDictionaryProtease[protein].insert(std::get<0>(kvp) );
                    }
                    else {
                        algDictionaryProtease.emplace(protein, std::unordered_set<ParsimonySequence*> {std::get<0>(kvp)});
                    }

                    std::unordered_map<Protein*, std::unordered_set<std::string>>::const_iterator algDictionary_iterator = algDictionary.find(protein);
                    if (algDictionary_iterator != algDictionary.end()) {
                        algDictionary[protein].insert(std::get<0>(kvp)->getSequence());
                    }
                    else {
                        std::unordered_set<std::string> usp4 = {std::get<0>(kvp)->getSequence()};
                        algDictionary.emplace(protein, usp4);
                    }
                }
            }

            // *** greedy algorithm loop
#ifdef ORIG
            int numNewSeqs = algDictionary.Max([&] (std::any p) {
                    p->Value->Count;
                });
#endif
            int numNewSeqs = -1;
            for ( auto p = algDictionary.begin(); p!= algDictionary.end(); p++ ) {
                if ( numNewSeqs == -1 ) {
                    numNewSeqs = (int)std::get<1>(*p).size();
                }
                if ( numNewSeqs < (int)std::get<1>(*p).size() ) {
                    numNewSeqs = (int)std::get<1>(*p).size();
                }
            }

            while (numNewSeqs != 0) {
                // gets list of proteins with the most unaccounted-for peptide sequences
#ifdef ORIG
                auto possibleBestProteinList = algDictionary.Where([&] (std::any p) {
                        return p->Value->Count == numNewSeqs;
                    })->ToList();
#endif
                std::vector<std::tuple<Protein*, std::unordered_set<std::string>>> possibleBestProteinList;
                for ( auto p= algDictionary.begin(); p!= algDictionary.end(); p++ ) {
                    if ( (int)std::get<1>(*p).size() == numNewSeqs ) {
                        possibleBestProteinList.push_back(*p);
                    }
                }

                Protein *bestProtein = std::get<0>(possibleBestProteinList.front());

                // may need to select different protein in case of a greedy algorithm tie
                // the protein with the most total peptide sequences wins in this case (doesn't matter if parsimony
                // has grabbed them or not)
                if (possibleBestProteinList.size() > 1) {
                    int highestNumTotalPep = proteinToPepSeqMatch[bestProtein].size();
                    for (auto kvp : possibleBestProteinList) {
                        if ( (int)proteinToPepSeqMatch[std::get<0>(kvp)].size() > highestNumTotalPep) {
                            highestNumTotalPep = proteinToPepSeqMatch[std::get<0>(kvp)].size();
                            bestProtein = std::get<0>(kvp);
                        }
                    }
                }

                parsimoniousProteinList.push_back(bestProtein);
                
                // remove observed peptide seqs
                std::vector<ParsimonySequence*> temp;
                for ( auto p = algDictionaryProtease[bestProtein].begin();
                      p != algDictionaryProtease[bestProtein].end(); p++ ) {
                    temp.push_back(*p);
                }
                
                for (auto peptideSequence : temp) {
                    std::vector<Protein*> proteinsWithThisPeptide = peptideSequenceToProteins[peptideSequence];

                    for (auto protein : proteinsWithThisPeptide) {
                        algDictionary[protein].erase(peptideSequence->getSequence());
                        algDictionaryProtease[protein].erase(peptideSequence);
                    }
                }

                algDictionary.erase(bestProtein);
                algDictionaryProtease.erase(bestProtein);
#ifdef ORIG
                numNewSeqs = algDictionary.Any() ? algDictionary.Max([&] (std::any p) {
                        p->Value->Count;
                    }) : 0;
#endif
                numNewSeqs = 0;
                if ( !algDictionary.empty() ) {
                    for ( auto p = algDictionary.begin(); p!= algDictionary.end(); p++ ) {
                        if ( numNewSeqs < (int)std::get<1>(*p).size() ) {
                            numNewSeqs = (int)std::get<1>(*p).size();
                        }
                    }
                }
            }

            // *** done with greedy algorithm
            // Parsimony stage 4: add back indistinguishable proteins (proteins that have identical peptide sets as
            // parsimonious proteins)
#ifdef ORIG
            auto allProteinsGroupedByNumPeptides = proteinToPepSeqMatch.GroupBy([&] (std::any p) {
                    p->Value->Count;
                });
#endif

            std::vector <std::tuple<Protein*, std::unordered_set<ParsimonySequence*>>> tmpSeqMatch;
            for ( auto p = proteinToPepSeqMatch.begin(); p!= proteinToPepSeqMatch.end(); p++ ) {
                tmpSeqMatch.push_back(std::make_tuple(std::get<0>(*p), std::get<1>(*p)));
            }
            std::function<bool(std::tuple<Protein*, std::unordered_set<ParsimonySequence*>>,std::tuple<Protein*, std::unordered_set<ParsimonySequence*>>)> f5 = [&]
                (std::tuple<Protein*, std::unordered_set<ParsimonySequence*>> l, std::tuple<Protein*, std::unordered_set<ParsimonySequence*>> r) {
                return std::get<1>(l).size() < std::get<1>(r).size(); } ;
            std::function<bool(std::tuple<Protein*, std::unordered_set<ParsimonySequence*>>,std::tuple<Protein*, std::unordered_set<ParsimonySequence*>>)> f6 = [&]
                (std::tuple<Protein*, std::unordered_set<ParsimonySequence*>> l, std::tuple<Protein*, std::unordered_set<ParsimonySequence*>> r) {
                return std::get<1>(l).size() < std::get<1>(r).size(); };
            std::vector<std::vector<std::tuple<Protein*, std::unordered_set<ParsimonySequence*>>>> allProteinsGroupedByNumPeptides = Group::GroupBy ( tmpSeqMatch, f5, f6);
            
#ifdef ORIG            
            auto parsimonyProteinsGroupedByNumPeptides = parsimoniousProteinList.GroupBy([&] (std::any p) {
                    proteinToPepSeqMatch[p].size();
                });
#endif
            std::function<bool(Protein*, Protein*)> f7 = [&] (Protein* l, Protein* r) {
                return  proteinToPepSeqMatch[l].size() <  proteinToPepSeqMatch[r].size(); } ;
            std::function<bool(Protein*, Protein*)> f8 = [&] (Protein *l, Protein *r) {
                return  proteinToPepSeqMatch[l].size() <  proteinToPepSeqMatch[r].size(); };
            std::vector<std::vector<Protein*>> parsimonyProteinsGroupedByNumPeptides = Group::GroupBy ( parsimoniousProteinList, f7, f8);
            
#ifdef ORIG
            auto indistinguishableProteins = new ConcurrentBag<Protein*>();
#endif
            std::unordered_set<Protein*> indistinguishableProteins;
            
            for (auto group : allProteinsGroupedByNumPeptides) {
#ifdef ORIG
                auto parsimonyProteinsWithSameNumPeptides = parsimonyProteinsGroupedByNumPeptides->FirstOrDefault([&] (std::any p) {
                        return p->Key == group->Key;
                    });
#endif
                auto parsimonyProteinsWithSameNumPeptides = *(parsimonyProteinsGroupedByNumPeptides.begin());
                auto list = group;

                if ( !parsimonyProteinsWithSameNumPeptides.empty() ) {
                    //ParallelOptions *tempVar = new ParallelOptions();
                    //tempVar->MaxDegreeOfParallelism = commonParameters::MaxThreadsToUsePerFile;
                    //Parallel::ForEach(Partitioner::Create(0, list.size()), tempVar, [&] (range, loopState) {
                    //        for (int i = range::Item1; i < range::Item2; i++) {
                    for ( int i = 0; i < (int) list.size(); i++ ) {
                        Protein *otherProtein = std::get<0>(list[i]);

                        for (auto parsimonyProtein : parsimonyProteinsWithSameNumPeptides) {
                            // if the two proteins have the same set of peptide sequences, they're indistinguishable
                            if (parsimonyProtein != otherProtein   &&
                                proteinToPepSeqMatch[parsimonyProtein] == proteinToPepSeqMatch[otherProtein] ) {
                                indistinguishableProteins.insert(otherProtein);
                            }
                        }
                    }
                    
                    //C# TO C++ CONVERTER TODO TASK: A 'delete tempVar' statement was not added since tempVar was
                    //passed to a method or constructor. Handle memory management manually.
                }
            }

            for (auto protein : indistinguishableProteins) {
                parsimoniousProteinList.push_back(protein);
            }

            //C# TO C++ CONVERTER TODO TASK: A 'delete indistinguishableProteins' statement was not added since
            //indistinguishableProteins was passed to a method or constructor. Handle memory management manually.
        }

        // Parsimony stage 5: remove peptide objects that do not have proteins in the parsimonious list
        for (PeptideSpectralMatch *psm : _allPsms) {
            // if this PSM has a protein in the parsimonious list, it removes the proteins NOT in the parsimonious list
            // otherwise, no proteins are removed (i.e., for PSMs that cannot be explained by a parsimonious protein,
            // no protein associations are removed)
#ifdef ORIG
            // if (psm->BestMatchingPeptides.Any([&] (std::any p) {
            //            std::find(parsimoniousProteinList.begin(), parsimoniousProteinList.end(), p::Peptide::Protein) != parsimoniousProteinList.end();
            //        })) {
#endif
            for ( auto p : psm->getBestMatchingPeptides() ) {
                if ( std::find(parsimoniousProteinList.begin(), parsimoniousProteinList.end(), std::get<1>(p)->getProtein() ) !=
                     parsimoniousProteinList.end()  ) {
                    psm->TrimProteinMatches(parsimoniousProteinList);
                }
            }
        }

        // construct protein groups
        std::vector<ProteinGroup*> proteinGroups = ConstructProteinGroups(uniquePeptides);

        // finished with parsimony
        return proteinGroups;
    }

    std::vector<ProteinGroup*> ProteinParsimonyEngine::ConstructProteinGroups(std::unordered_set<PeptideWithSetModifications*> &uniquePeptides)
    {
        std::vector<ProteinGroup*> proteinGroups;
        std::unordered_map<Protein*, std::unordered_set<PeptideWithSetModifications*>> proteinToPeptidesMatching;

        for (auto peptide : _fdrFilteredPeptides) {
            std::unordered_map<Protein*, std::unordered_set<PeptideWithSetModifications*>>::const_iterator proteinToPeptidesMatching_iterator = proteinToPeptidesMatching.find(peptide->getProtein());
            if (proteinToPeptidesMatching_iterator != proteinToPeptidesMatching.end()) {
                proteinToPeptidesMatching[peptide->getProtein()].insert(peptide);
            }
            else {
                std::unordered_set<PeptideWithSetModifications*> usp1 = {peptide};
                proteinToPeptidesMatching.emplace(peptide->getProtein(), usp1);
            }
        }

        for (auto kvp : proteinToPeptidesMatching) {
            auto allPeptidesHere = proteinToPeptidesMatching[std::get<0>(kvp)];
#ifdef ORIG
            auto uniquePeptidesHere = std::unordered_set<PeptideWithSetModifications*>(allPeptidesHere.Where([&] (std::any p) {
                        std::find(uniquePeptides.begin(), uniquePeptides.end(), p) != uniquePeptides.end();
                    }));
#endif
            std::unordered_set<PeptideWithSetModifications*> uniquePeptidesHere;
            for ( auto p : allPeptidesHere ) {
                if (std::find(uniquePeptides.begin(), uniquePeptides.end(), p) != uniquePeptides.end() ) {
                    uniquePeptidesHere.insert(p);
                }
            }

            std::unordered_set<Proteomics::Protein*> tmp1 = {std::get<0>(kvp)};
            
            auto tempVar = new ProteinGroup (tmp1, allPeptidesHere, uniquePeptidesHere);
            proteinGroups.push_back(tempVar);
        }

        for (auto proteinGroup : proteinGroups) {
#ifdef ORIG
            proteinGroup->AllPeptides.RemoveWhere([&] (std::any p) {
                    !proteinGroup->Proteins->Contains(p::Protein);
                });
#endif
            std::unordered_set<PeptideWithSetModifications*> newAllPeptides;// = proteinGroup->getAllPeptides();
            for ( auto p : proteinGroup->getAllPeptides() ){
                if ( proteinGroup->getProteins().find(p->getProtein()) != proteinGroup->getProteins().end() ) {
                    newAllPeptides.insert(p);
                }
            }
            proteinGroup->setAllPeptides(newAllPeptides);
            
            proteinGroup->setDisplayModsOnPeptides(_treatModPeptidesAsDifferentPeptides);
        }

        return proteinGroups;
    }
}
