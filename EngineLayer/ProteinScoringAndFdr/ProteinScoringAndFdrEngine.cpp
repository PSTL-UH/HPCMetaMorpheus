#include "ProteinScoringAndFdrEngine.h"
#include "ProteinScoringAndFdrResults.h"
#include "../PeptideSpectralMatch.h"
#include "../ProteinParsimony/ProteinGroup.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"

using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

    ProteinScoringAndFdrEngine::ProteinScoringAndFdrEngine(std::vector<ProteinGroup*> &proteinGroups,
                                                           std::vector<PeptideSpectralMatch*> &newPsms,
                                                           bool noOneHitWonders,
                                                           bool treatModPeptidesAsDifferentPeptides,
                                                           bool mergeIndistinguishableProteinGroups,
                                                           CommonParameters *commonParameters,
                                                           std::vector<std::string> &nestedIds) :
        MetaMorpheusEngine(commonParameters, nestedIds), NewPsms(newPsms),
        NoOneHitWonders(noOneHitWonders),
        TreatModPeptidesAsDifferentPeptides(treatModPeptidesAsDifferentPeptides),
        MergeIndistinguishableProteinGroups(mergeIndistinguishableProteinGroups),
        ProteinGroups(proteinGroups)
    {
    }
    
    MetaMorpheusEngineResults *ProteinScoringAndFdrEngine::RunSpecific()
    {
        ProteinScoringAndFdrResults *myAnalysisResults = new ProteinScoringAndFdrResults(this);
        
        ScoreProteinGroups(ProteinGroups, NewPsms);
        myAnalysisResults->SortedAndScoredProteinGroups = DoProteinFdr(ProteinGroups);
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete myAnalysisResults' statement was not added since
        //myAnalysisResults was used in a 'return' or 'throw' statement.
        return myAnalysisResults;
    }
    
    std::string ProteinScoringAndFdrEngine::StripDecoyIdentifier(const std::string &proteinGroupName)
    {
        return proteinGroupName.find("DECOY_") != std::string::npos ?
            StringHelper::replace(proteinGroupName, "DECOY_", "") : proteinGroupName;
    }
    
    void ProteinScoringAndFdrEngine::ScoreProteinGroups(std::vector<ProteinGroup*> &proteinGroups,
                                                        std::vector<PeptideSpectralMatch*> &psmList)
    {
        // add each protein groups PSMs
        auto peptideToPsmMatching = std::unordered_map<PeptideWithSetModifications*,
                                                       std::unordered_set<PeptideSpectralMatch*>>();
        for (auto psm : psmList)
        {
            if (psm->getFdrInfo()->getQValueNotch() <= 0.01 && psm->getFdrInfo()->getQValue() <= 0.01)
            {
                if ((TreatModPeptidesAsDifferentPeptides && psm->getFullSequence() != "") ||
                    (!TreatModPeptidesAsDifferentPeptides && psm->getBaseSequence() != ""))
                {
#ifdef ORIG
                    //for (auto pepWithSetMods : psm->BestMatchingPeptides->Select([&] (std::any p)  {
                    //            p::Peptide;
                    //        }))
#endif
                    std::vector<PeptideWithSetModifications *> tmpPep;
                    for ( auto p: psm->getBestMatchingPeptides() ) {
                        tmpPep.push_back(std::get<1>(p));
                    }
                    for ( auto pepWithSetMods : tmpPep )
                    {
                        std::unordered_set<PeptideSpectralMatch*> psmsForThisPeptide;
                        std::unordered_map<PeptideWithSetModifications*, std::unordered_set<PeptideSpectralMatch*>>::const_iterator peptideToPsmMatching_iterator = peptideToPsmMatching.find(pepWithSetMods);
                        if (peptideToPsmMatching_iterator == peptideToPsmMatching.end())
                        {
                            //psmsForThisPeptide = peptideToPsmMatching_iterator->second;
                            peptideToPsmMatching.emplace(pepWithSetMods, std::unordered_set<PeptideSpectralMatch*> {psm});
                        }
                        else
                        {
                            psmsForThisPeptide = peptideToPsmMatching_iterator->second;
                            psmsForThisPeptide.insert(psm);
                        }
                    }
                }
            }
        }
        
        for (auto proteinGroup : proteinGroups)
        {
            std::vector<PeptideWithSetModifications*> pepsToRemove;
            for (auto peptide : proteinGroup->getAllPeptides())
            {
                // build PSM list for scoring
                std::unordered_set<PeptideSpectralMatch*> psms;
                std::unordered_map<PeptideWithSetModifications*, std::unordered_set<PeptideSpectralMatch*>>::const_iterator peptideToPsmMatching_iterator = peptideToPsmMatching.find(peptide);
                if (peptideToPsmMatching_iterator != peptideToPsmMatching.end())
                {
                    psms = peptideToPsmMatching_iterator->second;
#ifdef ORIG
                    proteinGroup->getAllPsmsBelowOnePercentFDR().UnionWith(psms);
#endif
                    auto tmp = proteinGroup->getAllPsmsBelowOnePercentFDR();
                    for ( auto p: psms ){
                        bool found = false;
                        for ( auto q: tmp ) {
                            if ( p == q ) {
                                found = true;
                                break;
                            }
                        }
                        if ( !found ) {
                            tmp.insert(p);
                        }
                    }
                    proteinGroup->setAllPsmsBelowOnePercentFDR(tmp);
                }
                else
                {
                    //psms = peptideToPsmMatching_iterator->second;
                    pepsToRemove.push_back(peptide);
                }
            }
#ifdef ORIG            
            proteinGroup->getAllPeptides().ExceptWith(pepsToRemove);
#endif
            auto tmp1 =  proteinGroup->getAllPeptides();
            for ( auto p: tmp1 ) {
                for ( auto q: pepsToRemove) {
                    if ( p == q ) {
                        tmp1.erase (p);
                        break;
                    }
                }
            }
            proteinGroup->setAllPeptides(tmp1);
#ifdef ORIG            
            proteinGroup->getUniquePeptides().ExceptWith(pepsToRemove);
#endif
            auto tmp2 = proteinGroup->getUniquePeptides();
            for ( auto p: tmp2 ) {
                for ( auto q: pepsToRemove) {
                    if ( p == q ) {
                        tmp1.erase (p);
                        break;
                    }
                }
            }
            proteinGroup->setUniquePeptides(tmp2);
        }
        
        // score the group
        for (auto proteinGroup : proteinGroups)
        {
            proteinGroup->Score();
        }
        
        if (MergeIndistinguishableProteinGroups)
        {
            // merge protein groups that are indistinguishable after scoring
#ifdef ORIG
            auto pg = proteinGroups.OrderByDescending([&] (std::any p) {
                    p::ProteinGroupScore;
                }).ToList();
#endif
            auto pg = proteinGroups;
            std::sort(pg.begin(), pg.end(), [&] (ProteinGroup *l, ProteinGroup *r) {
                    return l->getProteinGroupScore() > r->getProteinGroupScore();
                });

            for (int i = 0; i < ((int)pg.size() - 1); i++)
            {
                if (pg[i]->getProteinGroupScore() == pg[i + 1]->getProteinGroupScore() &&
                    pg[i]->getProteinGroupScore() != 0)
                {
#ifdef ORIG
                    auto pgsWithThisScore = pg.Where([&] (std::any p)  {
                            return p->ProteinGroupScore == pg[i].ProteinGroupScore;
                        }).ToList();
#endif
                    std::vector<ProteinGroup*> pgsWithThisScore;
                    for ( auto p: pg ) {
                        if ( p->getProteinGroupScore() == pg[i]->getProteinGroupScore() ) {
                            pgsWithThisScore.push_back(p);
                        }
                    }
                    
                    // check to make sure they have the same peptides, then merge them
                    for (auto p : pgsWithThisScore)
                    {
#ifdef ORIG
                        auto seqs1 = std::unordered_set<std::string>(p.AllPeptides->Select([&] (std::any x)  {
                                    return x::FullSequence + x::DigestionParams::Protease;
                                }));
                        auto seqs2 = std::unordered_set<std::string>(pg[i].AllPeptides->Select([&] (std::any x) {
                                    return x::FullSequence + x::DigestionParams::Protease;
                                }));
#endif
                        std::set<std::string> seqs1;
                        for ( auto x: p->getAllPeptides() ) {
                            std::string s = x->getFullSequence() + x->getDigestionParams()->getProtease()->ToString();
                            seqs1.insert(s);
                        }
                        
                        std::set<std::string> seqs2;
                        for ( auto x: pg[i]->getAllPeptides()) {
                            std::string s = x->getFullSequence() + x->getDigestionParams()->getProtease()->ToString();
                            seqs1.insert(s);
                        }

#ifdef ORIG
                        //if (p != pg[i] && seqs1.SetEquals(seqs2))
#endif
                        if( p != pg[i] && seqs1 == seqs2 ) 
                        {
                            pg[i]->MergeProteinGroupWith(p);
                        }
                    }
                }
            }
        }
        
        // remove empty protein groups (peptides were too poor quality or group was merged)
#ifdef ORIG
        proteinGroups.RemoveAll([&] (std::any p)  {
                return p->ProteinGroupScore == 0;
            });
#endif
        for ( auto  p=0; p < (int)proteinGroups.size(); p++ ) {
            if ( proteinGroups[p]->getProteinGroupScore() == 0 ) {
                proteinGroups.erase(proteinGroups.begin() + p);
            }
        }
        
        // calculate sequence coverage
        for (auto proteinGroup : proteinGroups)
        {
            proteinGroup->CalculateSequenceCoverage();
        }
    }
    
    std::vector<ProteinGroup*> ProteinScoringAndFdrEngine::DoProteinFdr(std::vector<ProteinGroup*> &proteinGroups)
    {
        if (NoOneHitWonders)
        {
            if (TreatModPeptidesAsDifferentPeptides)
            {
#ifdef ORIG
                proteinGroups = proteinGroups.Where([&] (std::any p)   {
                        return p::IsDecoy ||
                        (std::unordered_set<std::string>(p::AllPeptides->Select([&] (std::any x)
                {
                    x::FullSequence;
                })))->size() > 1;
                    }).ToList();
#endif
                std::vector<ProteinGroup*> pG = proteinGroups;
                proteinGroups.clear();
                for ( auto p:  proteinGroups ){
                    int size=0;
                    for ( auto x: p->getAllPeptides() ) {
                        if ( x->getFullSequence().length() != 0 ) {
                            size++;
                        }
                    }
                    if (p->getIsDecoy() || size > 1 ) {
                        proteinGroups.push_back(p);
                    }
                }
            }
            else
            {
#ifdef ORIG
                //    proteinGroups = proteinGroups.Where([&] (std::any p)     {
                //        return p::IsDecoy ||
                //        (std::unordered_set<std::string>(p::AllPeptides->Select([&] (std::any x)
                // {
                //    x::BaseSequence;
                //})))->size() > 1;
                //    }).ToList();
            
#endif                
                std::vector<ProteinGroup*> pG = proteinGroups;
                proteinGroups.clear();
                for ( auto p:  proteinGroups ){
                    int size=0;
                    for ( auto x: p->getAllPeptides() ) {
                        if ( x->getBaseSequence().length() != 0  ) {
                            size++;
                        }
                    }
                    if (p->getIsDecoy() || size > 1 ) {
                        proteinGroups.push_back(p);
                    }
                }
            }
        }
        
        // pair decoys and targets by accession
        // then use the best peptide notch-QValue as the score for the protein group
        std::unordered_map<std::string, std::vector<ProteinGroup*>> accessionToProteinGroup;
        for (auto pg : proteinGroups)
        {
            for (auto protein : pg->getProteins())
            {
                std::string stippedAccession = StripDecoyIdentifier(protein->getAccession());
                
                std::vector<ProteinGroup*> groups;
                std::unordered_map<std::string, std::vector<ProteinGroup*>>::const_iterator accessionToProteinGroup_iterator = accessionToProteinGroup.find(stippedAccession);
                if (accessionToProteinGroup_iterator != accessionToProteinGroup.end())
                {
                    groups = accessionToProteinGroup_iterator->second;
                    groups.push_back(pg);
                }
                else
                {
                    //groups = accessionToProteinGroup_iterator->second;
                    accessionToProteinGroup.emplace(stippedAccession, std::vector<ProteinGroup*> {pg});
                }
            }
#ifdef ORIG            
            pg->setBestPeptideScore(pg->getAllPsmsBelowOnePercentFDR().Max([&] (std::any psm) {
                        psm::Score;
                    }));
            pg->setBestPeptideQValue(pg->getAllPsmsBelowOnePercentFDR().Min([&] (std::any psm) {
                        psm::FdrInfo::QValueNotch;
                    }));
#endif
            double maxv = (*(pg->getAllPsmsBelowOnePercentFDR().begin()))->getScore();
            for ( auto psm : pg->getAllPsmsBelowOnePercentFDR() ) {
                if ( psm->getScore() > maxv ) {
                    maxv = psm->getScore();
                }
            }
            pg->setBestPeptideScore(maxv);

            double minv = (*(pg->getAllPsmsBelowOnePercentFDR().begin()))->getScore();
            for ( auto psm : pg->getAllPsmsBelowOnePercentFDR() ) {
                if ( psm->getScore() < minv ) {
                    minv = psm->getScore();
                }
            }
            pg->setBestPeptideQValue(minv);
        }
        
        // pick the best notch-QValue for each paired accession
        for (auto accession : accessionToProteinGroup)
        {
            if (std::get<1>(accession).size() > 1)
            {
#ifdef ORIG
                auto pgList = accession.Value->OrderBy([&] (std::any p) {
                        p::BestPeptideQValue;
                    }).ThenByDescending([&] (std::any p) {
                            p::BestPeptideScore;
                        }).ToList();
#endif
                auto pgList = std::get<1>(accession);
                std::sort(pgList.begin(), pgList.end(), [&] (ProteinGroup *l, ProteinGroup *r) {
                        if (l->getBestPeptideQValue() < r->getBestPeptideQValue() ) return true;
                        if (l->getBestPeptideQValue() > r->getBestPeptideQValue() ) return false;
                        
                        if (l->getBestPeptideScore() > r->getBestPeptideScore() ) return true;
                        return false;
                    });
                    
                        
                auto pgToUse = pgList.front(); // pick lowest notch QValue and remove the rest
                pgList.erase(pgList.begin());
#ifdef ORIG
                proteinGroups = proteinGroups.Except(pgList).ToList();
#endif
                for ( auto q: pgList ) {
                    for ( auto p=0; p < (int)proteinGroups.size(); p++ ) {
                        if ( proteinGroups[p] == q ) {
                            proteinGroups.erase (proteinGroups.begin() + p);
                        }
                    }
                }
            }
        }
        
        // order protein groups by notch-QValue
#ifdef ORIG
            auto sortedProteinGroups = proteinGroups.OrderBy([&] (std::any b) {
                b::BestPeptideQValue;
            }).ThenByDescending([&] (std::any p) {
                    p::BestPeptideScore;
		}).ToList();
#endif
            auto sortedProteinGroups = proteinGroups;
            std::sort(sortedProteinGroups.begin(), sortedProteinGroups.end(), [&] (ProteinGroup *l, ProteinGroup *r) {
                    if (l->getBestPeptideQValue() < r->getBestPeptideQValue() ) return true;
                    if (l->getBestPeptideQValue() > r->getBestPeptideQValue() ) return false;
                    
                    if (l->getBestPeptideScore() > r->getBestPeptideScore() ) return true;
                    return false;
                });

        
        // calculate protein QValues
        int cumulativeTarget = 0;
        int cumulativeDecoy = 0;
        
        for (auto proteinGroup : sortedProteinGroups)
        {
            if (proteinGroup->getIsDecoy())
            {
                cumulativeDecoy++;
            }
            else
            {
                cumulativeTarget++;
            }
            
            proteinGroup->setCumulativeTarget(cumulativeTarget);
            proteinGroup->setCumulativeDecoy(cumulativeDecoy);
            proteinGroup->setQValue(static_cast<double>(cumulativeDecoy) / cumulativeTarget);
        }
        
        return sortedProteinGroups;
    }
}
