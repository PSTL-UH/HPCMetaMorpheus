#include "IndexingEngine.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "../MetaMorpheusException.h"
#include "../EventArgs/ProgressEventArgs.h"
#include "../GlobalVariables.h"

#include "IndexingResults.h"
#include "bankersrounding.h"

#include "omp.h"

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
            int oldPercentProgress = 0;
            ReportEngineProgress("Digesting proteins...", oldPercentProgress);            

            auto result = new IndexingResults(this);
            std::vector<PeptideWithSetModifications*>& peptidesSortedByMass = result->getPeptideIndex();            
            int ProteinListsize  = (int) ProteinList.size();
            
            // digest database
#pragma omp parallel
            {
                int tid = omp_get_thread_num();
                int num_threads = omp_get_num_threads();
                double progress = 0;
            
#pragma omp for schedule(guided)
                for ( int i = 0; i < ProteinListsize; i++ ) {
                    progress++;
                    std::vector<Modification*> *fixedmodis = const_cast<std::vector<Modification*>*> (&FixedModifications);
                    std::vector<Modification*> *varmodis = const_cast<std::vector<Modification*>*>(&VariableModifications);
                    auto localPeptides = ProteinList[i]->Digest(commonParameters->getDigestionParams(),
                                                                *fixedmodis, *varmodis);
                    
#pragma omp critical
                    {
                        peptidesSortedByMass.insert(peptidesSortedByMass.end(), localPeptides.begin(), localPeptides.end() );
                    }
                    
                    if ( tid == 0 )
                    {
                        auto percentProgress = static_cast<int>(((progress* num_threads) / ProteinListsize) * 100);                        
                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportEngineProgress("Digesting proteins...", percentProgress);
                        }
                    }
                }
            }
                
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
            oldPercentProgress = 0;
            ReportEngineProgress("Fragmenting peptides...", oldPercentProgress);

            int peptidesSortedByMasssize = (int) peptidesSortedByMass.size();

#pragma omp parallel
            {
                int tid = omp_get_thread_num();
                int num_threads = omp_get_num_threads();
                std::map<int, std::vector<int>> local_fragmentIndex;
                double progress = 0;

#pragma omp for schedule(guided)
                for (int peptideId = 0; peptideId < peptidesSortedByMasssize; peptideId++)
                {
                    auto t = peptidesSortedByMass[peptideId]->Fragment(commonParameters->getDissociationType(),
                                           commonParameters->getDigestionParams()->getFragmentationTerminus());
                    std::vector<double> fragmentMasses;
                    for ( auto m: t ) {
                        fragmentMasses.push_back(m->NeutralMass);
                        delete m;
                    }
                    for (auto theoreticalFragmentMass : fragmentMasses)
                    {
                        if (theoreticalFragmentMass < MaxFragmentSize && theoreticalFragmentMass > 0)
                        {
                            int fragmentBin = static_cast<int>(BankersRounding::round(theoreticalFragmentMass * FragmentBinsPerDalton));
                            
                            if ( local_fragmentIndex.find(fragmentBin) == local_fragmentIndex.end())
                            {
                                local_fragmentIndex[fragmentBin] = {peptideId};
                            }
                            else
                            {
                                local_fragmentIndex[fragmentBin].push_back(peptideId);
                            }
                        }
                    }
                    
                    if ( tid == 0 )
                    {
                        progress++;
                        auto percentProgress = static_cast<int>(((progress* num_threads) / peptidesSortedByMasssize) * 100);
                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportEngineProgress("Fragmenting peptides...", oldPercentProgress);
                        }
                    }
                }

#pragma omp critical
                {
                    for ( auto p = local_fragmentIndex.begin(); p != local_fragmentIndex.end() ; p++ ) {
                        int key = p->first;
                        fragmentIndex[key].insert(fragmentIndex[key].end(), p->second.begin(), p->second.end() );
                    }
                }
            }

            // fragmentIndex needs to be sorted to match the result of the sequential case. Not sure whether it is
            // required for correctness though.
#pragma omp parallel for schedule(guided)
            for ( int i = 0; i < fragmentIndex.size(); i++  ) {
                if (!fragmentIndex[i].empty() ) {
                    std::sort(fragmentIndex[i].begin(), fragmentIndex[i].end() );
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

                double progress = 0;
                oldPercentProgress = 0;
                ReportEngineProgress ("Creating precursor index...", oldPercentProgress);

                // no OpenMP directives for this loop because of the break statement below.
                for (int i = 0; i < peptidesSortedByMasssize; i++)
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
                    auto percentProgress = static_cast<int>(((progress) / peptidesSortedByMasssize) * 100);
                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportEngineProgress ("Creating precursor index...", oldPercentProgress);
                    }
                }
            }

            return result;
        }
    }
}
