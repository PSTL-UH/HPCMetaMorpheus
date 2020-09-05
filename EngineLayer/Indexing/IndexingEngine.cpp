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
                                       std::vector<std::string> nestedIds, int verbosityLevel) :
            MetaMorpheusEngine(commonParams, nestedIds, verbosityLevel),
            ProteinList(proteinList), FixedModifications(fixedModifications),
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
            auto var1 = decoyType;
            sb->appendLine("Search Decoys: " + DecoyTypeToString(var1));
            sb->appendLine("Number of proteins: " + std::to_string(ProteinList.size()));
            sb->appendLine("Number of fixed mods: " + std::to_string(FixedModifications.size()));
            sb->appendLine("Number of variable mods: " + std::to_string(VariableModifications.size()));
            auto var2 = commonParameters->getDissociationType();
            sb->appendLine("Dissociation Type: " + GetDissocationType::GetDissocationTypeAsString(var2));
            sb->appendLine("protease: " + commonParameters->getDigestionParams()->getProtease()->ToString());
            auto var3 = commonParameters->getDigestionParams()->getInitiatorMethionineBehavior();
            sb->appendLine("initiatorMethionineBehavior: " + InitiatorMethionineBehaviorToString(var3));
            sb->appendLine("maximumMissedCleavages: " + std::to_string(commonParameters->getDigestionParams()->getMaxMissedCleavages()));
            sb->appendLine("minPeptideLength: " + std::to_string(commonParameters->getDigestionParams()->getMinPeptideLength()));
            sb->appendLine("maxPeptideLength: " + std::to_string(commonParameters->getDigestionParams()->getMaxPeptideLength()));
            sb->appendLine("maximumVariableModificationIsoforms: " + std::to_string(commonParameters->getDigestionParams()->getMaxModificationIsoforms()));
            auto var4 = commonParameters->getDigestionParams()->getFragmentationTerminus();
            sb->appendLine("digestionTerminus: " + FragmentationTerminusToString(var4));
            sb->appendLine("maxModsForEachPeptide: " + std::to_string(commonParameters->getDigestionParams()->getMaxModsForPeptide()));
            auto var5 = commonParameters->getDigestionParams()->getSearchModeType();
            sb->appendLine("cleavageSpecificity: " + CleavageSpecificityExtension::GetCleavageSpecificityAsString(var5));
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
            
            auto result = new IndexingResults(this);
            std::vector<PeptideWithSetModifications*>& peptidesSortedByMass = result->getPeptideIndex();            
#ifdef ORIG
            //ParallelOptions *tempVar = new ParallelOptions();
            //tempVar->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
            //Parallel::ForEach(Partitioner::Create(0, ProteinList.size()), tempVar, [&] (range, loopState)  {
            //        std::vector<PeptideWithSetModifications*> localPeptides;
#endif
            // digest database
            for ( int i = 0; i < (int)ProteinList.size(); i++ ) {
                progress++;
                std::vector<Modification*> *fixedmodis = const_cast<std::vector<Modification*>*> (&FixedModifications);
                std::vector<Modification*> *varmodis = const_cast<std::vector<Modification*>*>(&VariableModifications);
                auto localPeptides = ProteinList[i]->Digest(commonParameters->getDigestionParams(),
                                                            *fixedmodis, *varmodis);

                peptidesSortedByMass.insert(peptidesSortedByMass.end(), localPeptides.begin(), localPeptides.end() );

                auto percentProgress = static_cast<int>((progress / ProteinList.size()) * 100);                        
                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    std::vector<std::string> vs(nestedIds.begin(), nestedIds.end() );
                    ProgressEventArgs tempVar2(percentProgress, "Digesting proteins...", vs);
                    ReportProgress(&tempVar2);
                }
            }
            //});
            
            // sort peptides by mass
            std::sort(peptidesSortedByMass.begin(), peptidesSortedByMass.end(), [&]
                      (PeptideWithSetModifications* l, PeptideWithSetModifications* r) {
                          return l->getMonoisotopicMass() < r->getMonoisotopicMass();
                      });
            
            // create fragment index
            std::vector<std::vector<int>>& fragmentIndex = result->getFragmentIndex();            
            try
            {
                fragmentIndex = std::vector<std::vector<int>>(static_cast<int>(std::ceil(MaxFragmentSize)) *
                                                              FragmentBinsPerDalton + 1);
            }
            catch (int e1)
            {
                std::string s = "Max fragment mass too large for indexing engine; try \"Classic Search\" mode, "
                    "or make the maximum fragment mass smaller";
                throw MetaMorpheusException(s);
            }
            
            // populate fragment index
            progress = 0;
            oldPercentProgress = 0;
            for (int peptideId = 0; peptideId < (int) peptidesSortedByMass.size(); peptideId++)
            {
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
            
            std::vector<std::vector<int>>& precursorIndex = result->getPrecursorIndex();
            
            if (GeneratePrecursorIndex)
            {
                // create precursor index
                try
                {
                    precursorIndex = std::vector<std::vector<int>>(static_cast<int>(std::ceil(MaxFragmentSize)) *
                                                                   FragmentBinsPerDalton + 1);
                }
                catch (int e2)
                {
                    std::string s = "Max precursor mass too large for indexing engine; try \"Classic Search\" mode,"
                                                " or make the maximum fragment mass smaller";
                    throw MetaMorpheusException(s);
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
            
            return result;
        }
    }
}
