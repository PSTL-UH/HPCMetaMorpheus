#include "IndexingEngine.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "../MetaMorpheusException.h"
#include "../EventArgs/ProgressEventArgs.h"
#include "../GlobalVariables.h"

#include "IndexingResults.h"
#include "bankersrounding.h"

using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;
namespace EngineLayer
{
    namespace Indexing
    {
        
        IndexingEngine::IndexingEngine(std::vector<Protein*> &proteinList, std::vector<Modification*> &variableModifications,
                                       std::vector<Modification*> &fixedModifications, int currentPartition, DecoyType decoyType,
                                       CommonParameters *commonParams, double maxFragmentSize,
                                       bool generatePrecursorIndex, std::vector<std::string> &proteinDatabases,
                                       std::vector<std::string> &nestedIds) :
            MetaMorpheusEngine(commonParams, nestedIds), ProteinList(proteinList), FixedModifications(fixedModifications),
            VariableModifications(variableModifications), CurrentPartition(currentPartition + 1),
            decoyType(decoyType), MaxFragmentSize(maxFragmentSize),
            GeneratePrecursorIndex(generatePrecursorIndex), ProteinDatabases(proteinDatabases)
        {
        }

        std::string IndexingEngine::ToString()
        {
            auto sb = new StringBuilder();
#ifdef ORIG
            sb->appendLine("Databases: " + std::string::Join(",", ProteinDatabases.OrderBy([&] (std::any p)   {
                            p->Name;
			})->Select([&] (std::any p){
                                delete sb;
				return p->Name + ":" + p->CreationTime;
                            }
                            )));
#endif
            std::sort(ProteinDatabases.begin(), ProteinDatabases.end() );

#ifdef NOT_NOW
            // not sorting for now by creation time
            std::vector<std::string> vs={"Databases: "};
            for ( auto p : ProteinDatabases ) {
                vs.push_back((p->Name + ":" + std::to_string(p->CreationTime)));
            }
            sb->appendLine(StringHelper::join(vs, ','));
#endif
            sb->appendLine(StringHelper::join(ProteinDatabases, ','));
            
            sb->appendLine("Partitions: " + std::to_string(CurrentPartition) + "/" + std::to_string(commonParameters->getTotalPartitions()));
            sb->appendLine("Precursor Index: " + StringHelper::toString(GeneratePrecursorIndex));
            sb->appendLine("Search Decoys: " + std::to_string(static_cast<int>(decoyType)));
            sb->appendLine("Number of proteins: " + std::to_string(ProteinList.size()));
            sb->appendLine("Number of fixed mods: " + std::to_string(FixedModifications.size()));
            sb->appendLine("Number of variable mods: " + std::to_string(VariableModifications.size()));
            sb->appendLine("Dissociation Type: " + std::to_string(static_cast<int>(commonParameters->getDissociationType())));
            sb->appendLine("protease: " + commonParameters->getDigestionParams()->getProtease()->ToString());
            sb->appendLine("initiatorMethionineBehavior: " + std::to_string(static_cast<int>(commonParameters->getDigestionParams()->getInitiatorMethionineBehavior())));
            sb->appendLine("maximumMissedCleavages: " + commonParameters->getDigestionParams()->getMaxMissedCleavages());
            sb->appendLine("minPeptideLength: " + commonParameters->getDigestionParams()->getMinPeptideLength());
            sb->appendLine("maxPeptideLength: " + commonParameters->getDigestionParams()->getMaxPeptideLength());
            sb->appendLine("maximumVariableModificationIsoforms: " + commonParameters->getDigestionParams()->getMaxModificationIsoforms());
            sb->appendLine("digestionTerminus: " + std::to_string(static_cast<int>(commonParameters->getDigestionParams()->getFragmentationTerminus())));
            sb->appendLine("maxModsForEachPeptide: " + commonParameters->getDigestionParams()->getMaxModsForPeptide());
            sb->appendLine("cleavageSpecificity: " + std::to_string(static_cast<int>(commonParameters->getDigestionParams()->getSearchModeType())));
            sb->appendLine("specificProtease: " + commonParameters->getDigestionParams()->getSpecificProtease()->ToString());
            
#ifdef ORIG
            sb->append("Localizeable mods: " + ProteinList.Select([&] (std::any b)   {
                        b::OneBasedPossibleLocalizedModifications->Count;
                    }).Sum());
#endif
            int sum=0;
            for ( auto b: ProteinList ) {
                sum += b->getOneBasedPossibleLocalizedModifications().size();
            }
            sb->append("Localizeable mods: " + std::to_string(sum));

            std::string s = sb->toString();
            delete sb;
            return s;
        }

        MetaMorpheusEngineResults *IndexingEngine::RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            
            // digest database
            std::vector<PeptideWithSetModifications*> globalPeptides;
            
#ifdef ORIG
            //ParallelOptions *tempVar = new ParallelOptions();
            //tempVar->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
            //Parallel::ForEach(Partitioner::Create(0, ProteinList.size()), tempVar, [&] (range, loopState)  {
            //        std::vector<PeptideWithSetModifications*> localPeptides;
            //        
            //        for (int i = range::Item1; i < range::Item2; i++) {
            //        //Stop loop if canceled
            //        if (GlobalVariables::getStopLoops()) {
            //               loopState::Stop();
            //               return;
            //        }
                        
            // localPeptides.AddRange(ProteinList[i]->Digest(commonParameters->getDigestionParams(),
            //                                                          FixedModifications, VariableModifications));
                                                
#endif
            for ( int i = 0; i < (int)ProteinList.size(); i++ ) {
                progress++;
                std::vector<Modification*> *fixedmodis = const_cast<std::vector<Modification*>*> (&FixedModifications);
                std::vector<Modification*> *varmodis = const_cast<std::vector<Modification*>*>(&VariableModifications);
                std::vector<PeptideWithSetModifications*> t = ProteinList[i]->Digest(commonParameters->getDigestionParams(),
                                                *fixedmodis, *varmodis);
                for ( auto p : t ) {
                    globalPeptides.push_back(p);
                }
                auto percentProgress = static_cast<int>((progress / ProteinList.size()) * 100);
                        
                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    std::vector<std::string> vs(nestedIds.begin(), nestedIds.end() );
                    ProgressEventArgs tempVar2(percentProgress, "Digesting proteins...", vs);
                    ReportProgress(&tempVar2);
                }
            }
            
#ifdef ORIG
            {
                std::lock_guard<std::mutex> lock(globalPeptides);
                globalPeptides.AddRange(localPeptides);
            }
#endif
            //});
            
            // sort peptides by mass
#ifdef ORIG
            auto peptidesSortedByMass = globalPeptides.OrderBy([&] (std::any p)  {
                    p::MonoisotopicMass;
                }).ToList();
#endif
            std::vector<PeptideWithSetModifications*> peptidesSortedByMass(globalPeptides.begin(), globalPeptides.end());
            std::sort(peptidesSortedByMass.begin(), peptidesSortedByMass.end(), [&] (PeptideWithSetModifications* l, PeptideWithSetModifications* r) {
                    return l->getMonoisotopicMass() < r->getMonoisotopicMass();
                });
            globalPeptides.clear();
            
            // create fragment index
            std::vector<std::vector<int>> fragmentIndex;
            
            try
            {
                fragmentIndex = std::vector<std::vector<int>>(static_cast<int>(std::ceil(MaxFragmentSize)) * FragmentBinsPerDalton + 1);
            }
            catch (int e1)
            {
                throw MetaMorpheusException("Max fragment mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
            }
            
            // populate fragment index
            progress = 0;
            oldPercentProgress = 0;
            for (int peptideId = 0; peptideId < (int) peptidesSortedByMass.size(); peptideId++)
            {
 #ifdef ORIG
                //auto fragmentMasses = peptidesSortedByMass[peptideId].Fragment(commonParameters->getDissociationType(),
                //              commonParameters->getDigestionParams()->FragmentationTerminus)->Select([&] (std::any m)   {
                //                      m::NeutralMass;
                //                  }).ToList();
#endif
                auto t = peptidesSortedByMass[peptideId]->Fragment(commonParameters->getDissociationType(),
                                                                   commonParameters->getDigestionParams()->getFragmentationTerminus());
                std::vector<double> fragmentMasses;
                for ( auto m: t ) {
                    fragmentMasses.push_back(m->NeutralMass);
                }
                for (auto theoreticalFragmentMass : fragmentMasses)
                {
                    if (theoreticalFragmentMass < MaxFragmentSize && theoreticalFragmentMass > 0)
                    {
                        int fragmentBin = static_cast<int>(BankersRounding::round(theoreticalFragmentMass * FragmentBinsPerDalton));
                        
                        if (fragmentIndex[fragmentBin].empty())
                        {
                            fragmentIndex[fragmentBin] = {peptideId};
                        }
                        else
                        {
                            fragmentIndex[fragmentBin].push_back(peptideId);
                        }
                    }
                }
                
                progress++;
                auto percentProgress = static_cast<int>((progress / peptidesSortedByMass.size()) * 100);
                
                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    std::vector<std::string> vs(nestedIds.begin(), nestedIds.end() );
                    ProgressEventArgs tempVar3(percentProgress, "Fragmenting peptides...", vs);
                    ReportProgress(&tempVar3);
                }
            }
            
            std::vector<std::vector<int>> precursorIndex;
            
            if (GeneratePrecursorIndex)
            {
                // create precursor index
                try
                {
                    precursorIndex = std::vector<std::vector<int>>(static_cast<int>(std::ceil(MaxFragmentSize)) * FragmentBinsPerDalton + 1);
                }
                catch (int e2)
                {
                    throw MetaMorpheusException("Max precursor mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
                }
                progress = 0;
                oldPercentProgress = 0;
                std::vector<std::string> vs(nestedIds.begin(), nestedIds.end() );
                ProgressEventArgs tempVar4(0, "Creating precursor index...", vs);
                ReportProgress(&tempVar4);
                
                for (int i = 0; i < (int)peptidesSortedByMass.size(); i++)
                {
                    double mass = peptidesSortedByMass[i]->getMonoisotopicMass();
                    if (!std::isnan(mass))
                    {
                        if (mass > MaxFragmentSize) //if the precursor is larger than the index allows, then stop adding precursors
                        {
                            break;
                        }
                        
                        int precursorBin = static_cast<int>(BankersRounding::round(mass * FragmentBinsPerDalton));
                        
                        if (precursorIndex[precursorBin].empty())
                        {
                            precursorIndex[precursorBin] = {i};
                        }
                        else
                        {
                            precursorIndex[precursorBin].push_back(i);
                        }
                    }
                    progress++;
                    auto percentProgress = static_cast<int>((progress / peptidesSortedByMass.size()) * 100);
                    
                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        std::vector<std::string> vs(nestedIds.begin(), nestedIds.end() );
                        ProgressEventArgs tempVar5(percentProgress, "Creating precursor index...", vs);
                        ReportProgress(&tempVar5);
                    }
                }
            }
            
            return new IndexingResults(peptidesSortedByMass, fragmentIndex, precursorIndex, this);
        }
    }
}
