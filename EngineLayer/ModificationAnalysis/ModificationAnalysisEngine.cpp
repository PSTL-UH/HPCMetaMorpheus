#include "ModificationAnalysisEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "ModificationAnalysisResults.h"

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

#include "Group.h"

namespace EngineLayer
{
    namespace ModificationAnalysis
    {
        
        ModificationAnalysisEngine::ModificationAnalysisEngine(std::vector<PeptideSpectralMatch*> &newPsms,
                                                               CommonParameters *commonParameters,
                                                               std::vector<std::string> &nestedIds) :
            MetaMorpheusEngine(commonParameters, nestedIds), NewPsms(newPsms)
        {
        }
        
        MetaMorpheusEngineResults *ModificationAnalysisEngine::RunSpecific()
        {
            Status("Running modification analysis...");
            
            ModificationAnalysisResults *myAnalysisResults = new ModificationAnalysisResults(this);
            
#ifdef ORIG
            auto confidentTargetPsms = NewPsms.Where([&] (std::any b)    {
                    return b::FdrInfo::QValue <= 0.01 && !b::IsDecoy;
                }).ToList();
#endif
            std::vector<PeptideSpectralMatch*> confidentTargetPsms;
            for ( auto b: NewPsms ) {
                if ( b->getFdrInfo()->getQValue() < 0.01 && !b->getIsDecoy() ) {
                    confidentTargetPsms.push_back(b);
                }
            }

            // For the database ones, only need un-ambiguous protein and location in protein
#ifdef ORIG
            auto forObserved = confidentTargetPsms.Where([&] (std::any b){
                    return b::ProteinAccession != nullptr && b::OneBasedEndResidueInProtein != nullptr &&
                    b::OneBasedStartResidueInProtein != nullptr;
                });
#endif
            std::vector<PeptideSpectralMatch*> forObserved;
            for ( auto b: confidentTargetPsms ) {
                if ( b->getProteinAccession().length() != 0           &&
                     !b->getOneBasedEndResidueInProtein().has_value() &&
                     !b->getOneBasedStartResidueInProtein().has_value() ) {
                    forObserved.push_back(b);
                }
            }
            
            // For the unambiguously localized ones, need FullSequence and un-ambiguous protein and location in protein
#ifdef ORIG
            auto forUnambiguouslyLocalized = confidentTargetPsms.Where([&] (std::any b)	{
                    return b::FullSequence != nullptr && b::ProteinAccession != nullptr &&
                    b::OneBasedEndResidueInProtein != nullptr && b::OneBasedStartResidueInProtein != nullptr;
                });
#endif
            std::vector<PeptideSpectralMatch*> forUnambiguouslyLocalized;
            for ( auto b: confidentTargetPsms ) {
                if ( b->getFullSequence().length() != 0               &&
                     b->getProteinAccession().length() != 0           &&
                     !b->getOneBasedEndResidueInProtein().has_value() &&
                     !b->getOneBasedStartResidueInProtein().has_value() ) {
                    forUnambiguouslyLocalized.push_back(b);
                }
            }
#ifdef ORIG
            //**DEBUG
            std::vector<PeptideSpectralMatch*> toby;
            for (auto psm : confidentTargetPsms)
            {
                if (psm->getFullSequence() != "" && psm->getProteinAccession() != "" &&
                    !psm->getOneBasedEndResidueInProtein().has_value()               &&
                    !psm->getOneBasedStartResidueInProtein().has_value() )
                {
                    toby.push_back(psm);
                }
            }
                        
            int i = forUnambiguouslyLocalized->Count() + toby.size()();
            //**END DEBUG
#endif            
            // For the localized but ambiguous ones, need FullSequence
#ifdef ORIG
            auto forAmbiguousButLocalized = confidentTargetPsms.Where([&] (std::any b)     {
                    return b::FullSequence != nullptr && !(b::ProteinAccession != nullptr &&
                                                           b::OneBasedEndResidueInProtein != nullptr &&
                                                           b::OneBasedStartResidueInProtein != nullptr);
                }).GroupBy([&] (std::any b)			{
                        b::FullSequence;
                    });
#endif
            std::vector<PeptideSpectralMatch*> tmpPsms;
            for ( auto b: confidentTargetPsms ) {
                if ( b->getFullSequence().length() != 0               &&
                     b->getProteinAccession().length() != 0           &&
                     !b->getOneBasedEndResidueInProtein().has_value() &&
                     !b->getOneBasedStartResidueInProtein().has_value() ) {
                    tmpPsms.push_back(b);
                }
            }
            
            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f1 = [&]
                (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                return l->getFullSequence() < r->getFullSequence(); } ;
            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f2 = [&]
                (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                return l->getFullSequence() != r->getFullSequence(); } ;
            std::vector<std::vector<PeptideSpectralMatch*>> forAmbiguousButLocalized = Group::GroupBy ( tmpPsms, f1, f2);

            
            // For unlocalized but identified modifications, skip ones with full sequences!
#ifdef ORIG
            auto forUnlocalized = confidentTargetPsms.Where([&] (std::any b)         {
                    return b::BaseSequence != nullptr &&
                           b->FullSequence == nullptr &&
                           b::ModsIdentified != nullptr;
                }).GroupBy([&] (std::any b)  {
                        (b::BaseSequence,
                         std::string::Join(" ", b::ModsIdentified->Values->OrderBy([&] (std::any c)  {
                                     return c;
                                 })));
                    });
#endif
            tmpPsms.clear();
            for ( auto b: confidentTargetPsms ) {
                if ( b->getFullSequence().length() != 0               &&
                     b->getBaseSequence().length() != 0           &&
                     b->getModsIdentified().size() != 0 ){
                    tmpPsms.push_back(b);
                }
            }
            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f3 = [&]
                (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                std::string sl, sr;
                std::string sb = " ";
                std::vector<std::string> vl, vr;
                for ( auto p = l->getModsIdentified().begin(); p != l->getModsIdentified().end(); p++  ) {
                    vl.push_back (std::to_string(p->second));
                }
                for ( auto p = r->getModsIdentified().begin(); p != r->getModsIdentified().end(); p++  ) {
                    vr.push_back (std::to_string(p->second));
                }                
                sl = StringHelper::join(vl, sb);
                sr = StringHelper::join(vr, sb);
                return l->getBaseSequence() < r->getBaseSequence() && sl < sr; } ;

            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f4 = [&]
                (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                std::string sl, sr;
                std::string sb = " ";
                std::vector<std::string> vl, vr;
                for ( auto p = l->getModsIdentified().begin(); p != l->getModsIdentified().end(); p++  ) {
                    vl.push_back (std::to_string(p->second));
                }
                for ( auto p = r->getModsIdentified().begin(); p != r->getModsIdentified().end(); p++  ) {
                    vr.push_back (std::to_string(p->second));
                }                
                sl = StringHelper::join(vl, sb);
                sr = StringHelper::join(vr, sb);

                return l->getBaseSequence() != r->getBaseSequence() && sl != sr ; } ;
            std::vector<std::vector<PeptideSpectralMatch*>> forUnlocalized = Group::GroupBy ( tmpPsms, f3, f4);

            
            // For chemical formulas of modifications, skip ones with full sequences and identified mods!
#ifdef ORIG
            auto forChemicalFormulas = confidentTargetPsms.Where([&] (std::any b)      {
                    return b::BaseSequence != nullptr && b->FullSequence == nullptr &&
                    b->ModsIdentified == nullptr && b::ModsChemicalFormula != nullptr;
                }).GroupBy([&] (std::any b)	{
                        (b::BaseSequence, b::ModsChemicalFormula);
                    });
#endif
            tmpPsms.clear();
            for ( auto b: confidentTargetPsms ) {
                if ( b->getFullSequence().length() != 0   &&
                     b->getBaseSequence().length() == 0   &&
                     b->getModsIdentified().size() == 0   &&
                     b->getModsChemicalFormula() != nullptr ){
                    tmpPsms.push_back(b);
                }
            }
            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f5 = [&]
                (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                return l->getFullSequence() < r->getFullSequence() &&
                l->getModsChemicalFormula() < r->getModsChemicalFormula() ; } ;
            std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f6 = [&]
                (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
                return l->getFullSequence() != r->getFullSequence() &&
                l->getModsChemicalFormula() != r->getModsChemicalFormula() ; } ;
            std::vector<std::vector<PeptideSpectralMatch*>> forChemicalFormulas = Group::GroupBy ( tmpPsms, f5, f6);
            
            // We do not want to double-count modifications. Hence the HashSet!!!
            ModTuple_set modsOnProteins;
            for (auto psm : forObserved)
            {
                auto singlePeptide = std::get<1>(psm->getBestMatchingPeptides().front());
#ifdef ORIG
                //for (auto modInProtein : singlePeptide->getProtein()->getOneBasedPossibleLocalizedModifications().Where([&] (std::any b) {
                //        return b::Key >= singlePeptide->OneBasedStartResidueInProtein &&
                //                b::Key <= singlePeptide->OneBasedEndResidueInProtein;
                //        }))
#endif
                for (auto modInProtein : singlePeptide->getProtein()->getOneBasedPossibleLocalizedModifications() )
                {
                    if ( std::get<0>(modInProtein) < singlePeptide->getOneBasedStartResidueInProtein() ||
                         std::get<0>(modInProtein) > singlePeptide->getOneBasedEndResidueInProtein() ) {
                        continue;
                    }
                    for (auto huh : std::get<1>(modInProtein) )
                    {
                        modsOnProteins.insert(std::make_tuple(singlePeptide->getProtein()->getAccession(),
                                                              huh->getIdWithMotif(),
                                                              std::get<0>(modInProtein)));
                    }
                }
            }
            
            // We do not want to double-count modifications. Hence the HashSet!!!
            ModTuple_set modsSeenAndLocalized;
            for (auto psm : forUnambiguouslyLocalized)
            {
                auto singlePeptide = std::get<1>(psm->getBestMatchingPeptides().front());
                for (auto nice : singlePeptide->getAllModsOneIsNterminus())
                {
                    int locInProtein;
                    if (std::get<0>(nice) == 1)
                    {
                        locInProtein = singlePeptide->getOneBasedStartResidueInProtein();
                    }
                    else if (std::get<0>(nice) == singlePeptide->getLength() + 2)
                    {
                        locInProtein = singlePeptide->getOneBasedEndResidueInProtein();
                    }
                    else
                    {
                        locInProtein = singlePeptide->getOneBasedStartResidueInProtein() + std::get<0>(nice) - 2;
                    }
                    modsSeenAndLocalized.insert(std::make_tuple(singlePeptide->getProtein()->getAccession(),
                                                                std::get<1>(nice)->getIdWithMotif(),
                                                                locInProtein));
                }
            }
            
            // Might have some double counting...
            std::unordered_map<std::string, int> ambiguousButLocalizedModsSeen;
#ifdef ORIG
            //for (auto representativePsm : forAmbiguousButLocalized->Select([&] (std::any b)   {
            //            b::First();
            //        }))
#endif
            for (auto r : forAmbiguousButLocalized ) 
            {
                auto  representativePsm = r.front();
                
                for (auto modCountKvp : representativePsm->getModsIdentified())
                {
                    if (ambiguousButLocalizedModsSeen.find(std::get<0>(modCountKvp) ) != ambiguousButLocalizedModsSeen.end())
                    {
                        ambiguousButLocalizedModsSeen[std::get<0>(modCountKvp)] += std::get<1>(modCountKvp);
                    }
                    else
                    {
                        ambiguousButLocalizedModsSeen.emplace(std::get<0>(modCountKvp), std::get<1>(modCountKvp));
                    }
                }
            }
            
            // Might have some double counting...
            std::unordered_map<std::string, int> unlocalizedMods;
#ifdef ORIG
            // for (auto representativePsm : forUnlocalized->Select([&] (std::any b)      {
            //            b::First();
            //        }))
#endif
            for ( auto r: forUnlocalized ) 
            {
                auto representativePsm = r.front();
                
                for (auto modCountKvp : representativePsm->getModsIdentified() )
                {
                    if (unlocalizedMods.find(std::get<0>(modCountKvp)) != unlocalizedMods.end())
                    {
                        unlocalizedMods[std::get<0>(modCountKvp)] += std::get<1>(modCountKvp);
                    }
                    else
                    {
                        unlocalizedMods.emplace(std::get<0>(modCountKvp), std::get<1>(modCountKvp));
                    }
                }
            }
            
            // Might have some double counting...
            std::unordered_map<ChemicalFormula*, int> unlocalizedFormulas;
#ifdef ORIG
            //for (auto representativePsm : forChemicalFormulas->Select([&] (std::any b)   {
            //            b::First();
            //        }))
#endif
            for ( auto r :  forChemicalFormulas ) 
            {
                auto representativePsm = r.front();
                
                if (unlocalizedFormulas.find(representativePsm->getModsChemicalFormula()) != unlocalizedFormulas.end())
                {
                    unlocalizedFormulas[representativePsm->getModsChemicalFormula()] += 1;
                }
                else
                {
                    unlocalizedFormulas.emplace(representativePsm->getModsChemicalFormula(), 1);
                }
            }
            
#ifdef ORIG
            myAnalysisResults->setCountOfEachModSeenOnProteins(modsOnProteins.GroupBy([&] (std::any b)   {
                        b::Item2;
                    }).ToDictionary([&] (std::any b)    {
                            b::Key;
			}, [&] (std::any b)	{
                            b->Count();
			}));
#endif
            std::unordered_map<std::string, int> modsOP;
            for ( auto r = modsOnProteins.begin(); r != modsOnProteins.end(); r++ ) {
                if (modsOP.find(std::get<1>(*r)) != modsOP.end()) {
                    modsOP[std::get<1>(*r)] += 1;
                }
                else  {
                    modsOP.emplace(std::get<1>(*r), 1);
                }
            }
            myAnalysisResults->setCountOfEachModSeenOnProteins( modsOP);
            
#ifdef ORIG
            myAnalysisResults->setCountOfModsSeenAndLocalized(modsSeenAndLocalized.GroupBy([&] (std::any b) {
                        b::Item2;
			}).ToDictionary([&] (std::any b) {
				b::Key;
                            }, [&] (std::any b)	{
				b->Count();
                            }));
#endif
            std::unordered_map<std::string, int> modsSAL;
            for ( auto r = modsSeenAndLocalized.begin(); r != modsSeenAndLocalized.end(); r++ ) {
                if (modsSAL.find(std::get<1>(*r)) != modsSAL.end()) {
                    modsSAL[std::get<1>(*r)] += 1;
                }
                else  {
                    modsSAL.emplace(std::get<1>(*r), 1);
                }
            }
            myAnalysisResults->setCountOfModsSeenAndLocalized(modsSAL);
            
            myAnalysisResults->setCountOfAmbiguousButLocalizedModsSeen(ambiguousButLocalizedModsSeen);
            myAnalysisResults->setCountOfUnlocalizedMods(unlocalizedMods);
            myAnalysisResults->setCountOfUnlocalizedFormulas(unlocalizedFormulas);
            
            //C# TO C++ CONVERTER TODO TASK: A 'delete myAnalysisResults' statement was not added since
            //myAnalysisResults was used in a 'return' or 'throw' statement.
            return myAnalysisResults;
        }
    }
}
