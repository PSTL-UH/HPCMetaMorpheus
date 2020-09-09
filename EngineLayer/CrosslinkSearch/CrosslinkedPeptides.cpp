#include "CrosslinkedPeptides.h"
#include "Crosslinker.h"

using namespace MassSpectrometry;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
    namespace CrosslinkSearch
    {
        std::vector<std::tuple<int, std::vector<Product*>>> CrosslinkedPeptide::XlGetTheoreticalFragments(DissociationType dissociationType,
                                                                                     Crosslinker *crosslinker,
                                                                                     std::vector<int> &possibleCrosslinkerPositions,
                                                                                     double otherPeptideMass,
                                                                                     PeptideWithSetModifications *peptide)
        {
            std::vector<std::tuple<int, std::vector<Product*>>> *retvec = new std::vector<std::tuple<int, std::vector<Product*>>>;
            std::vector<double> massesToLocalize;
            if (crosslinker->getCleavable())
            {
                massesToLocalize.push_back(crosslinker->getCleaveMassShort());
                massesToLocalize.push_back(crosslinker->getCleaveMassLong());
            }
            else
            {
                massesToLocalize.push_back(crosslinker->getTotalMass() + otherPeptideMass);
            }
            
            for (auto crosslinkerPosition : possibleCrosslinkerPositions)
            {
                std::vector<Product*> theoreticalProducts;
                std::unordered_set<double> masses;
                
                for (auto massToLocalize : massesToLocalize)
                {
                    std::string oId="", accs="", modType="", featType="", locRestr="Unassigned";
                    ModificationMotif *targ=nullptr;
                    ChemicalFormula *chem=nullptr;
                    
                    std::unordered_map<int, Modification*> testMods =
                        {{crosslinkerPosition + 1, new Modification( oId, accs, modType, featType, targ, locRestr,
                                                                     chem, std::make_optional<double>(massToLocalize))}};
                    //{{crosslinkerPosition + 1, new Modification(_monoisotopicMass: massToLocalize)}};
                    auto testPeptide = new PeptideWithSetModifications(peptide->getProtein(),
                                                                       peptide->getDigestionParams(),
                                                                       peptide->getOneBasedStartResidueInProtein(),
                                                                       peptide->getOneBasedEndResidueInProtein(),
                                                                       peptide->getCleavageSpecificityForFdrCategory(),
                                                                       peptide->getPeptideDescription(),
                                                                       peptide->getMissedCleavages(),
                                                                       testMods,
                                                                       peptide->NumFixedMods);
                    
                    // add fragmentation ions for this crosslinker position guess
                    for (auto fragment : testPeptide->Fragment(dissociationType, FragmentationTerminus::Both))
                    {
                        if (std::find(masses.begin(), masses.end(), fragment->NeutralMass) == masses.end())
                        {
                            theoreticalProducts.push_back(fragment);
                            masses.insert(fragment->NeutralMass);
                        }
                    }
                    
                    // add signature ions
                    if (crosslinker->getCleavable())
                    {
                        auto tempVar = new Product(ProductType::M,
                                                   new NeutralTerminusFragment(FragmentationTerminus::None,
                                                                               peptide->getMonoisotopicMass() + massToLocalize,
                                                                               peptide->getLength(),
                                                                               peptide->getLength()), 0);
                        theoreticalProducts.push_back(tempVar);
                    }
                    
                    delete testPeptide;
                }
                
                retvec->push_back(std::make_tuple<int, std::vector<Product*>>((int)crosslinkerPosition,
                                                                              (std::vector<Product*>)theoreticalProducts));
            }
            return *retvec;
        }
        
        XLumap CrosslinkedPeptide::XlLoopGetTheoreticalFragments(DissociationType dissociationType, Modification *loopMass,
                                                                 std::vector<int> &modPos, PeptideWithSetModifications *peptide)
        {
            XLumap AllTheoreticalFragmentsLists;
            auto originalFragments = peptide->Fragment(dissociationType, FragmentationTerminus::Both);
            for (auto position1 : modPos)
            {
                for (auto position2 : modPos)
                {
                    if (position2 <= position1)
                    {
                        continue;
                    }
                    
                    // add N and C terminal fragments that do not contain the loop
                    std::vector<int> loopPositions = {position1, position2};
#ifdef ORIG
                    std::vector<Product*> loopFragments = originalFragments.Where([&] (std::any p)  {
                            return p::TerminusFragment->Terminus == FragmentationTerminus::N &&
                            p::TerminusFragment::AminoAcidPosition < position1               ||
                            p::TerminusFragment->Terminus == FragmentationTerminus::C        &&
                            p::TerminusFragment::AminoAcidPosition > position2;
                        }).ToList();
#endif
                    std::vector<Product*> loopFragments;
                    for ( auto p: originalFragments ) {
                        if (  (p->TerminusFragment->Terminus == FragmentationTerminus::N &&
                               p->TerminusFragment->AminoAcidPosition < position1)        ||
                              (p->TerminusFragment->Terminus == FragmentationTerminus::C &&
                               p->TerminusFragment->AminoAcidPosition > position2)  ) {
                            loopFragments.push_back(p);
                        }
                    }
                    // add N-terminal fragments containing the loop
                    std::unordered_map<int, Modification*> modDict;

                    if (!peptide->getAllModsOneIsNterminus().empty())
                    {
#ifdef ORIG
                        double combinedModMass = loopMass->getMonoisotopicMass().value() + peptide->getAllModsOneIsNterminus().Where([&] (std::any v) {
                                return v::Key <= position2 + 1;
                            }).Sum([&] (std::any p) {
                                    p->Value->MonoisotopicMass->Value;
                                });
#endif
                        double combinedModMass = loopMass->getMonoisotopicMass().value();                        
                        for ( auto v: peptide->getAllModsOneIsNterminus()) {
                            if ( std::get<0>(v) <= position2+1 ) {
                                combinedModMass += std::get<1>(v)->getMonoisotopicMass().value();
                            }
                        }
                        // Modification *combined = new Modification(_monoisotopicMass: combinedModMass);
                        
                        std::string oId="", accs="", modType="", featType="", locRestr="Unassigned";
                        ModificationMotif *targ=nullptr;
                        ChemicalFormula *chem=nullptr;
                        Modification *combined = new Modification( oId, accs, modType, featType, targ, locRestr,
                                                                   chem, std::make_optional<double>(combinedModMass));
                        
                        modDict.emplace(position1 + 1, combined);
                        
#ifdef ORIG
                        for (auto mod : peptide->AllModsOneIsNterminus.Where([&] (std::any m)  {
                                    //C# TO C++ CONVERTER TODO TASK: A 'delete combined' statement was not added
                                    //since combined was passed to a method or constructor. Handle memory management manually.
                                    return m::Key > position2 + 1;
                                }));
                        {
                            modDict.emplace(mod::Key, mod->Value);
                        }
#endif
                        for ( auto mod: peptide->getAllModsOneIsNterminus() ) {
                            if ( std::get<0>(mod) > position2+1 ) {
                                modDict.emplace(std::get<0>(mod), std::get<1>(mod));
                            }
                        }
                         
                        //C# TO C++ CONVERTER TODO TASK: A 'delete combined' statement was not added since
                        //combined was passed to a method or constructor. Handle memory management manually.
                    }
                    else
                    {
                        modDict.emplace(position1 + 1, loopMass);
                    }
                    PeptideWithSetModifications *peptideWithLoop = new PeptideWithSetModifications(peptide->getProtein(),
                                                                                         peptide->getDigestionParams(),
                                                                                         peptide->getOneBasedStartResidueInProtein(),
                                                                                         peptide->getOneBasedEndResidueInProtein(),
                                                                                         peptide->getCleavageSpecificityForFdrCategory(),
                                                                                         peptide->getPeptideDescription(),
                                                                                         peptide->getMissedCleavages(),
                                                                                         modDict,
                                                                                         peptide->NumFixedMods);
#ifdef ORIG
                    loopFragments.AddRange(peptideWithLoop->Fragment(dissociationType, FragmentationTerminus::Both).Where([&] (std::any p) {
                                delete peptideWithLoop;
                                return p::TerminusFragment->Terminus == FragmentationTerminus::N &&
                                    p::TerminusFragment::AminoAcidPosition >= position2;
                            }));
#endif
                    auto tmp = peptideWithLoop->Fragment(dissociationType, FragmentationTerminus::Both);
                    for ( auto p : tmp ) {
                        if ( p->TerminusFragment->Terminus == FragmentationTerminus::N &&
                             p->TerminusFragment->AminoAcidPosition >= position2 ) {
                            loopFragments.push_back(p);
                        }
                    }
                    
                    // add C-terminal fragments containing the loop
                    modDict.clear();
                    if (!peptide->getAllModsOneIsNterminus().empty())
                    {
#ifdef ORIG
                        double combinedModMass = loopMass->getMonoisotopicMass().value() + peptide->getAllModsOneIsNterminus().Where([&] (std::any v)      {
                                delete peptideWithLoop;
                                return v::Key >= position1 + 1;
                            }).Sum([&] (std::any p) {
                                    p->Value->MonoisotopicMass->Value;
                                });
#endif
                        double combinedModMass = loopMass->getMonoisotopicMass().value();                        
                        for ( auto v: peptide->getAllModsOneIsNterminus()) {
                            if ( std::get<0>(v) >= position1+1 ) {
                                combinedModMass += std::get<1>(v)->getMonoisotopicMass().value();
                            }
                        }

                        //Modification *combined = new Modification(_monoisotopicMass: combinedModMass);
                        std::string oId="", accs="", modType="", featType="", locRestr="Unassigned";
                        ModificationMotif *targ=nullptr;
                        ChemicalFormula *chem=nullptr;
                        Modification *combined = new Modification( oId, accs, modType, featType, targ, locRestr,
                                                                   chem, std::make_optional<double>(combinedModMass));
                        modDict.emplace(position2 + 1, combined);
                        
#ifdef ORIG
                        for (auto mod : peptide->AllModsOneIsNterminus.Where([&] (std::any m) {
                                    //C# TO C++ CONVERTER TODO TASK: A 'delete combined' statement was not added
                                    //since combined was passed to a method or constructor. Handle memory management manually.
                                    delete peptideWithLoop;
                                    return m::Key < position1 + 1;
                                }));
                        {
                            modDict.emplace(mod::Key, mod->Value);
                        }
#endif
                        for ( auto mod: peptide->getAllModsOneIsNterminus() ) {
                            if ( std::get<0>(mod) < position1+1 ) {
                                modDict.emplace(std::get<0>(mod), std::get<1>(mod));
                            }
                        }
                        
                        //C# TO C++ CONVERTER TODO TASK: A 'delete combined' statement was not added since
                        //combined was passed to a method or constructor. Handle memory management manually.
                    }
                    else
                    {
                        modDict.emplace(position2 + 1, loopMass);
                    }
                    peptideWithLoop = new PeptideWithSetModifications(peptide->getProtein(),
                                                                      peptide->getDigestionParams(),
                                                                      peptide->getOneBasedStartResidueInProtein(),
                                                                      peptide->getOneBasedEndResidueInProtein(),
                                                                      peptide->getCleavageSpecificityForFdrCategory(),
                                                                      peptide->getPeptideDescription(),
                                                                      peptide->getMissedCleavages(),
                                                                      modDict, peptide->NumFixedMods);
#ifdef ORIG
                    loopFragments.AddRange(peptideWithLoop->Fragment(dissociationType, FragmentationTerminus::Both).Where([&] (std::any p) {
                                delete peptideWithLoop;
                                return p::TerminusFragment->Terminus == FragmentationTerminus::C &&
                                    p::TerminusFragment::AminoAcidPosition <= position1;
                            }));
#endif
                    auto tmp2 = peptideWithLoop->Fragment(dissociationType, FragmentationTerminus::Both);
                    for ( auto p : tmp2 ) {
                        if ( p->TerminusFragment->Terminus == FragmentationTerminus::C &&
                             p->TerminusFragment->AminoAcidPosition <= position1 ) {
                            loopFragments.push_back(p);
                        }
                    }
                    
                    
                    //AllTheoreticalFragmentsLists.emplace(loopPositions, loopFragments);
                    XLTuple thistuple = std::make_tuple(loopPositions[0],loopPositions[1]);
                    AllTheoreticalFragmentsLists.emplace(thistuple, loopFragments);
                    
                    delete peptideWithLoop;
                }
            }
            
            return AllTheoreticalFragmentsLists;
        }
    }
}
