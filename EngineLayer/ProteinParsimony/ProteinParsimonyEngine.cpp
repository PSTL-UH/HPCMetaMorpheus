#include "ProteinParsimonyEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"

using namespace EngineLayer::ProteinParsimony;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

    ProteinParsimonyEngine::ProteinParsimonyEngine(std::vector<PeptideSpectralMatch*> &allPsms, bool modPeptidesAreDifferent,
                                                   CommonParameters *commonParameters, std::vector<std::string> &nestedIds) :
        MetaMorpheusEngine(commonParameters, nestedIds),
        _fdrFilteredPeptides(std::unordered_set<PeptideWithSetModifications*>()), _allPsms(allPsms),
        _treatModPeptidesAsDifferentPeptides(modPeptidesAreDifferent)
    {
        
        if (!allPsms.Any()) 	{
            _fdrFilteredPsms = std::vector<PeptideSpectralMatch*>();
        }
        
        // parsimony will only use non-ambiguous, high-confidence PSMs
        // KEEP decoys and contaminants for use in parsimony!
        if (modPeptidesAreDifferent)
        {
                _fdrFilteredPsms = allPsms.Where([&] (std::any p) {
                        return p::FullSequence != nullptr && p::FdrInfo::QValue <= FdrCutoffForParsimony &&
                        p::FdrInfo::QValueNotch <= FdrCutoffForParsimony;
                    }).ToList();
        }
        else
        {
            _fdrFilteredPsms = allPsms.Where([&] (std::any p) {
                    return p::BaseSequence != nullptr && p::FdrInfo::QValue <= FdrCutoffForParsimony &&
                    p::FdrInfo::QValueNotch <= FdrCutoffForParsimony;
                }).ToList();
        }
        
        // if PSM is a decoy, add only decoy sequences; same for contaminants
        // peptides to use in parsimony = peptides observed in high-confidence PSMs
        
        for (auto psm : _fdrFilteredPsms)
        {
            if (psm->getIsDecoy())
            {
                for (auto peptide : psm->BestMatchingPeptides->Select([&] (std::any p)    {
                            p::Peptide;
                        }).Where([&] (std::any p)  {
                                p::Protein::IsDecoy;
                            }))
                {
                }
            }
        }
        else if (psm::IsContaminant)
        {
            for (auto peptide : psm::BestMatchingPeptides->Select([&] (std::any p)     {
                        p::Peptide;
                    }).Where([&] (std::any p) {
                            p::Protein::IsContaminant;
                        }))
            {
            }
        }
        else // PSM is target
        {
            for (auto peptide : psm::BestMatchingPeptides->Select([&] (std::any p)  {
                        p::Peptide;
                    }).Where([&] (std::any p)  {
                            return !p::Protein::IsDecoy && !p::Protein::IsContaminant;
                        }))
            {
            }
        }
        // we're storing all PSMs (not just FDR-filtered ones) here because we will remove some protein associations 
        // from low-confidence PSMs if they can be explained by a parsimonious protein
        _allPsms = allPsms;
    }
    
    ProteinParsimonyEngine::MetaMorpheusEngineResults* RunSpecific();
    {
        auto myAnalysisResults = new ProteinParsimonyResults(this);
        myAnalysisResults.ProteinGroups = RunProteinParsimonyEngine();
        
        return myAnalysisResults;
    }
    
    
    std::vector<ProteinGroup*> ProteinParsimonyEngine::RunProteinParsimonyEngine()
    {
        // parsimonious list of proteins built by this protein parsimony engine
        std::unordered_set<Protein*> parsimoniousProteinList = std::unordered_set<Protein*>();
        
        // list of peptides that can only be digestion products of one protein in the proteome (considering different
        // protease digestion rules)
        std::unordered_set<PeptideWithSetModifications*> uniquePeptides = std::unordered_set<PeptideWithSetModifications*>();
        
        // if there are no peptides observed, there are no proteins; return an empty list of protein groups
        if (_fdrFilteredPeptides->Count == 0) {
            return std::vector<ProteinGroup*>();
        }
        
        // Parsimony stage 0: create peptide-protein associations if needed because the user wants a
        // modification-agnostic parsimony
        if (!_treatModPeptidesAsDifferentPeptides) {
            for (auto protease : _fdrFilteredPsms::GroupBy([&] (std::any p) {
                        p::DigestionParams::Protease;
                    })) {
                std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>> sequenceWithPsms = std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>>();

                // for each protease, match the base sequence of each peptide to its PSMs
                for (PeptideSpectralMatch *psm : protease) {
                    List<PeptideSpectralMatch*> peptidesForThisBaseSequence;
                    std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>>::const_iterator sequenceWithPsms_iterator = sequenceWithPsms.find(psm.BaseSequence);
                    if (sequenceWithPsms_iterator != sequenceWithPsms.end()) {
                        peptidesForThisBaseSequence = sequenceWithPsms_iterator->second;
                        peptidesForThisBaseSequence->Add(psm);
                    }
                    else {
                        peptidesForThisBaseSequence = sequenceWithPsms_iterator->second;
                        sequenceWithPsms[psm->BaseSequence] = {psm};
                    }
                }

                // create new peptide-protein associations
                for (auto baseSequence : sequenceWithPsms) {
                    auto peptidesWithNotchInfo = baseSequence.Value->SelectMany([&] (std::any p) {
                            p::BestMatchingPeptides;
                        })->Distinct()->ToList();

                    // if the base seq has >1 PeptideWithSetMods object and has >0 mods, it might need to be matched to new proteins
                    if (peptidesWithNotchInfo.size() > 1 && peptidesWithNotchInfo.Any([&] (std::any p) {
                                return p::Peptide::NumMods > 0;
                            })) {
                        // list of proteins along with start/end residue in protein and the # missed cleavages
                        // this is needed to create new PeptideWithSetModification objects
                        auto peptideInProteinInfo = std::vector<std::tuple<Protein*, DigestionParams*, int, int, int, int>>();
                        for (auto peptide : peptidesWithNotchInfo) {
                            peptideInProteinInfo.push_back(std::tuple<Protein*, DigestionParams*, int, int, int, int>(peptide.Peptide::Protein, peptide.Peptide::DigestionParams, peptide.Peptide::OneBasedStartResidueInProtein, peptide.Peptide::OneBasedEndResidueInProtein, peptide.Peptide::MissedCleavages, peptide.Notch));
                        }
                        
                        // add the protein associations to the PSM
                        for (PeptideSpectralMatch *psm : baseSequence.Value) {
                            for (auto proteinInfo : peptideInProteinInfo) {
                                auto originalPep = psm->BestMatchingPeptides.First()->Peptide;
                                auto pep = new PeptideWithSetModifications(proteinInfo.Item1, proteinInfo.Item2, proteinInfo.Item3, proteinInfo.Item4, originalPep->CleavageSpecificityForFdrCategory, originalPep->PeptideDescription, proteinInfo.Item5, originalPep->AllModsOneIsNterminus, originalPep->NumFixedMods);
                                _fdrFilteredPeptides->Add(pep);
                                psm->AddProteinMatch((proteinInfo.Item6, pep));

                                //C# TO C++ CONVERTER TODO TASK: A 'delete pep' statement was not added since pep was
                                //passed to a method or constructor. Handle memory management manually.
                            }
                        }
                    }
                }
            }
        }

        // Parsimony stage 1: add proteins with unique peptides (for each protease)
        auto peptidesGroupedByProtease = _fdrFilteredPeptides::GroupBy([&] (std::any p) {
                p::DigestionParams::Protease;
            });
        for (auto peptidesForThisProtease : *peptidesGroupedByProtease) {
            std::unordered_map<std::string, std::vector<Protein*>> peptideSequenceToProteinsForThisProtease = std::unordered_map<std::string, std::vector<Protein*>>();
            std::unordered_map<std::string, std::vector<PeptideWithSetModifications*>> sequenceToPwsm = std::unordered_map<std::string, std::vector<PeptideWithSetModifications*>>();

            for (auto peptide : *peptidesForThisProtease) {
                std::string sequence = peptide->BaseSequence;
                if (_treatModPeptidesAsDifferentPeptides) {
                    //these and next set to full sequence but might be base sequence. treat modified as unique makes
                    //sense to use full
                    sequence = peptide->FullSequence;
                }

                List<Protein*> proteinsForThisPeptideSequence;
                std::unordered_map<std::string, std::vector<Protein*>>::const_iterator peptideSequenceToProteinsForThisProtease_iterator = peptideSequenceToProteinsForThisProtease.find(sequence);
                if (peptideSequenceToProteinsForThisProtease_iterator != peptideSequenceToProteinsForThisProtease.end()) {
                    proteinsForThisPeptideSequence = peptideSequenceToProteinsForThisProtease_iterator->second;
                    proteinsForThisPeptideSequence->Add(peptide->Protein);
                }
                else {
                    proteinsForThisPeptideSequence = peptideSequenceToProteinsForThisProtease_iterator->second;
                    peptideSequenceToProteinsForThisProtease.emplace(sequence, std::vector<Protein*> {peptide->Protein});
                }

                List<PeptideWithSetModifications*> peptidesForThisSequence;
                std::unordered_map<std::string, std::vector<PeptideWithSetModifications*>>::const_iterator sequenceToPwsm_iterator = sequenceToPwsm.find(sequence);
                if (sequenceToPwsm_iterator != sequenceToPwsm.end()) {
                    peptidesForThisSequence = sequenceToPwsm_iterator->second;
                    peptidesForThisSequence->Add(peptide);
                }
                else {
                    peptidesForThisSequence = sequenceToPwsm_iterator->second;
                    sequenceToPwsm.emplace(sequence, std::vector<PeptideWithSetModifications*> {peptide});
                }
            }

            for (auto uniquePeptide : peptideSequenceToProteinsForThisProtease.Where([&] (std::any p) {
                        return p->Value->Count == 1;
                    })) {
                // add the protein with the unique peptide to the parsimonious protein list
                Protein *proteinWithUniquePeptideSequence = uniquePeptide->Value->First();
                parsimoniousProteinList.insert(proteinWithUniquePeptideSequence);

                // add the unique peptide to the list of unique peptides
                PeptideWithSetModifications *uniquePwsm = sequenceToPwsm[uniquePeptide::Key].front();
                uniquePeptides.insert(uniquePwsm);
            }
        }

        // Parsimony stage 2: build the peptide-protein matching structure for the parsimony greedy algorithm
        // and remove all peptides observed by proteins with unique peptides
        std::unordered_map<ParsimonySequence*, std::vector<Protein*>> peptideSequenceToProteins = std::unordered_map<ParsimonySequence*, std::vector<Protein*>>();

        // this dictionary associates proteins w/ all peptide sequences (list will NOT shrink over time)
        // this is used in case of greedy algorithm ties to figure out which protein has more total peptides observed
        std::unordered_map<Protein*, std::unordered_set<ParsimonySequence*>> proteinToPepSeqMatch = std::unordered_map<Protein*, std::unordered_set<ParsimonySequence*>>();
        for (auto peptide : _fdrFilteredPeptides) {
            ParsimonySequence *sequence = new ParsimonySequence(peptide, _treatModPeptidesAsDifferentPeptides);

            List<Protein*> proteinsForThisPeptideSequence;
            std::unordered_map<ParsimonySequence*, std::vector<Protein*>>::const_iterator peptideSequenceToProteins_iterator = peptideSequenceToProteins.find(sequence);
            if (peptideSequenceToProteins_iterator != peptideSequenceToProteins.end()) {
                proteinsForThisPeptideSequence = peptideSequenceToProteins_iterator->second;
                proteinsForThisPeptideSequence->Add(peptide->Protein);
            }
            else {
                proteinsForThisPeptideSequence = peptideSequenceToProteins_iterator->second;
                peptideSequenceToProteins.emplace(sequence, std::vector<Protein*> {peptide->Protein});
            }

            TValue peptideSequences;
            std::unordered_map<Protein*, std::unordered_set<ParsimonySequence*>>::const_iterator proteinToPepSeqMatch_iterator = proteinToPepSeqMatch.find(peptide.Protein);
            if (proteinToPepSeqMatch_iterator != proteinToPepSeqMatch.end()) {
                peptideSequences = proteinToPepSeqMatch_iterator->second;
                peptideSequences->Add(sequence);
            }
            else {
                peptideSequences = proteinToPepSeqMatch_iterator->second;
                proteinToPepSeqMatch.emplace(peptide->Protein, std::unordered_set<ParsimonySequence*> {sequence});
            }

            //C# TO C++ CONVERTER TODO TASK: A 'delete sequence' statement was not added since sequence was passed
            //to a method or constructor. Handle memory management manually.
        }

        // remove the peptides observed by proteins with unique peptides
        std::unordered_set<ParsimonySequence*> toRemove = std::unordered_set<ParsimonySequence*>();
        for (auto seq : peptideSequenceToProteins) {
            bool observedAlready = seq.Value->Any([&] (std::any p) {
                    std::find(parsimoniousProteinList.begin(), parsimoniousProteinList.end(), p) != parsimoniousProteinList.end();
                });

            if (observedAlready) {
                toRemove.insert(seq.Key);
            }
        }
        for (auto sequence : toRemove) {
            peptideSequenceToProteins.erase(sequence);
        }

        if (peptideSequenceToProteins.Any()) {
            // Parsimony stage 3: greedy algorithm

            // dictionary with proteins as keys and list of associated peptide sequences as the values.
            // this data structure makes parsimony easier because the algorithm can look up a protein's peptides
            // to remove them from the list of available peptides. this list will shrink as the algorithm progresses
            auto algDictionary = std::unordered_map<Protein*, std::unordered_set<std::string>>();
            auto algDictionaryProtease = std::unordered_map<Protein*, std::unordered_set<ParsimonySequence*>>();
            for (auto kvp : peptideSequenceToProteins) {
                for (auto protein : kvp.Value) {
                    HashSet<ParsimonySequence*> peptideSequencesWithProtease;
                    std::unordered_map<Protein*, std::unordered_set<ParsimonySequence*>>::const_iterator algDictionaryProtease_iterator = algDictionaryProtease.find(protein);
                    if (algDictionaryProtease_iterator != algDictionaryProtease.end()) {
                        peptideSequencesWithProtease = algDictionaryProtease_iterator->second;
                        peptideSequencesWithProtease->Add(kvp.Key);
                    }
                    else {
                        peptideSequencesWithProtease = algDictionaryProtease_iterator->second;
                        algDictionaryProtease.emplace(protein, std::unordered_set<ParsimonySequence*> {kvp.Key});
                    }

                    HashSet<std::string> peptideSequences;
                    std::unordered_map<Protein*, std::unordered_set<std::string>>::const_iterator algDictionary_iterator = algDictionary.find(protein);
                    if (algDictionary_iterator != algDictionary.end()) {
                        peptideSequences = algDictionary_iterator->second;
                        peptideSequences->Add(kvp.Key->Sequence);
                    }
                    else {
                        peptideSequences = algDictionary_iterator->second;
                        algDictionary.emplace(protein, std::unordered_set<std::string> {kvp.Key->Sequence});
                    }
                }
            }

            // *** greedy algorithm loop
            int numNewSeqs = algDictionary.Max([&] (std::any p) {
                    p->Value->Count;
                });
            while (numNewSeqs != 0) {
                // gets list of proteins with the most unaccounted-for peptide sequences
                auto possibleBestProteinList = algDictionary.Where([&] (std::any p) {
                        return p->Value->Count == numNewSeqs;
                    })->ToList();

                Protein *bestProtein = possibleBestProteinList.front().Key;

                // may need to select different protein in case of a greedy algorithm tie
                // the protein with the most total peptide sequences wins in this case (doesn't matter if parsimony
                // has grabbed them or not)
                if (possibleBestProteinList.size() > 1) {
                    int highestNumTotalPep = proteinToPepSeqMatch[bestProtein].size();
                    for (auto kvp : possibleBestProteinList) {
                        if (proteinToPepSeqMatch[kvp.Key].size() > highestNumTotalPep) {
                            highestNumTotalPep = proteinToPepSeqMatch[kvp.Key].size();
                            bestProtein = kvp.Key;
                        }
                    }
                }

                parsimoniousProteinList.insert(bestProtein);

                // remove observed peptide seqs
                std::vector<ParsimonySequence*> temp = algDictionaryProtease[bestProtein].ToList();
                for (auto peptideSequence : temp) {
                    std::vector<Protein*> proteinsWithThisPeptide = peptideSequenceToProteins[peptideSequence];

                    for (auto protein : proteinsWithThisPeptide) {
                        algDictionary[protein].Remove(peptideSequence->Sequence);
                        algDictionaryProtease[protein].Remove(peptideSequence);
                    }
                }

                algDictionary.erase(bestProtein);
                algDictionaryProtease.erase(bestProtein);
                numNewSeqs = algDictionary.Any() ? algDictionary.Max([&] (std::any p) {
                        p->Value->Count;
                    }) : 0;
            }

            // *** done with greedy algorithm
            // Parsimony stage 4: add back indistinguishable proteins (proteins that have identical peptide sets as
            // parsimonious proteins)
            auto allProteinsGroupedByNumPeptides = proteinToPepSeqMatch.GroupBy([&] (std::any p) {
                    p->Value->Count;
                });
            auto parsimonyProteinsGroupedByNumPeptides = parsimoniousProteinList.GroupBy([&] (std::any p) {
                    proteinToPepSeqMatch[p].size();
                });
            auto indistinguishableProteins = new ConcurrentBag<Protein*>();

            for (auto group : *allProteinsGroupedByNumPeptides) {
                auto parsimonyProteinsWithSameNumPeptides = parsimonyProteinsGroupedByNumPeptides->FirstOrDefault([&] (std::any p) {
                        delete indistinguishableProteins;
                        return p->Key == group->Key;
                    });
                auto list = group->ToList();

                if (parsimonyProteinsWithSameNumPeptides != nullptr) {
                    ParallelOptions *tempVar = new ParallelOptions();
                    tempVar->MaxDegreeOfParallelism = commonParameters::MaxThreadsToUsePerFile;
                    Parallel::ForEach(Partitioner::Create(0, list.size()), tempVar, [&] (range, loopState) {
                            for (int i = range::Item1; i < range::Item2; i++) {
                                Protein *otherProtein = list[i].Key;

                                for (auto parsimonyProtein : *parsimonyProteinsWithSameNumPeptides) {
                                    // if the two proteins have the same set of peptide sequences, they're indistinguishable
                                    if (parsimonyProtein != otherProtein && proteinToPepSeqMatch[parsimonyProtein].SetEquals(proteinToPepSeqMatch[otherProtein])) {
                                        indistinguishableProteins->Add(otherProtein);
                                    }
                                }
                            }
                        });

                    //C# TO C++ CONVERTER TODO TASK: A 'delete tempVar' statement was not added since tempVar was
                    //passed to a method or constructor. Handle memory management manually.
                }
            }

            for (auto protein : *indistinguishableProteins) {
                parsimoniousProteinList.insert(protein);
            }

            //C# TO C++ CONVERTER TODO TASK: A 'delete indistinguishableProteins' statement was not added since
            //indistinguishableProteins was passed to a method or constructor. Handle memory management manually.
        }

        // Parsimony stage 5: remove peptide objects that do not have proteins in the parsimonious list
        for (PeptideSpectralMatch *psm : _allPsms) {
            // if this PSM has a protein in the parsimonious list, it removes the proteins NOT in the parsimonious list
            // otherwise, no proteins are removed (i.e., for PSMs that cannot be explained by a parsimonious protein,
            // no protein associations are removed)
            if (psm->BestMatchingPeptides.Any([&] (std::any p) {
                        std::find(parsimoniousProteinList.begin(), parsimoniousProteinList.end(), p::Peptide::Protein) != parsimoniousProteinList.end();
                    })) {
                psm->TrimProteinMatches(parsimoniousProteinList);
            }
        }

        // construct protein groups
        std::vector<ProteinGroup*> proteinGroups = ConstructProteinGroups(uniquePeptides);

        // finished with parsimony
        return proteinGroups;
    }

    std::vector<ProteinGroup*> ProteinParsimonyEngine::ConstructProteinGroups(std::unordered_set<PeptideWithSetModifications*> &uniquePeptides)
    {
        std::vector<ProteinGroup*> proteinGroups = std::vector<ProteinGroup*>();
        auto proteinToPeptidesMatching = std::unordered_map<Protein*, std::unordered_set<PeptideWithSetModifications*>>();

        for (auto peptide : _fdrFilteredPeptides) {
            HashSet<PeptideWithSetModifications*> peptidesHere;
            std::unordered_map<Protein*, std::unordered_set<PeptideWithSetModifications*>>::const_iterator proteinToPeptidesMatching_iterator = proteinToPeptidesMatching.find(peptide.Protein);
            if (proteinToPeptidesMatching_iterator != proteinToPeptidesMatching.end()) {
                peptidesHere = proteinToPeptidesMatching_iterator->second;
                peptidesHere->Add(peptide);
            }
            else {
                peptidesHere = proteinToPeptidesMatching_iterator->second;
                proteinToPeptidesMatching.emplace(peptide->Protein, std::unordered_set<PeptideWithSetModifications*> {peptide});
            }
        }

        for (auto kvp : proteinToPeptidesMatching) {
            auto allPeptidesHere = proteinToPeptidesMatching[kvp.Key];
            auto uniquePeptidesHere = std::unordered_set<PeptideWithSetModifications*>(allPeptidesHere.Where([&] (std::any p) {
                        std::find(uniquePeptides.begin(), uniquePeptides.end(), p) != uniquePeptides.end();
                    }));

            ProteinGroup tempVar({kvp.Key}, allPeptidesHere, uniquePeptidesHere);
            proteinGroups.push_back(&tempVar);
        }

        for (auto proteinGroup : proteinGroups) {
            proteinGroup->AllPeptides.RemoveWhere([&] (std::any p) {
                    !proteinGroup->Proteins->Contains(p::Protein);
                });
            proteinGroup->DisplayModsOnPeptides = _treatModPeptidesAsDifferentPeptides;
        }

        return proteinGroups;
    }
}
