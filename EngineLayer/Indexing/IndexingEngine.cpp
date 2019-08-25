#include "IndexingEngine.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "../MetaMorpheusException.h"
#include "../EventArgs/ProgressEventArgs.h"
#include "IndexingResults.h"

using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;
namespace EngineLayer
{
	namespace Indexing
	{

		IndexingEngine::IndexingEngine(std::vector<Protein*> &proteinList, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, int currentPartition, DecoyType *decoyType, CommonParameters *commonParams, double maxFragmentSize, bool generatePrecursorIndex, std::vector<FileInfo*> &proteinDatabases, std::vector<std::string> &nestedIds) : MetaMorpheusEngine(commonParams, nestedIds), ProteinList(proteinList), FixedModifications(fixedModifications), VariableModifications(variableModifications), CurrentPartition(currentPartition + 1), DecoyType(decoyType), MaxFragmentSize(maxFragmentSize), GeneratePrecursorIndex(generatePrecursorIndex), ProteinDatabases(proteinDatabases)
		{
		}

		std::string IndexingEngine::ToString()
		{
			auto sb = new StringBuilder();
			sb->appendLine("Databases: " + std::string::Join(",", ProteinDatabases.OrderBy([&] (std::any p)
			{
				p->Name;
			})->Select([&] (std::any p)
			{
			delete sb;
				return p->Name + ":" + p::CreationTime;
			})));
			sb->appendLine("Partitions: " + std::to_string(CurrentPartition) + "/" + std::to_string(commonParameters->getTotalPartitions()));
			sb->appendLine("Precursor Index: " + StringHelper::toString(GeneratePrecursorIndex));
			sb->appendLine("Search Decoys: " + DecoyType);
			sb->appendLine("Number of proteins: " + std::to_string(ProteinList.size()));
			sb->appendLine("Number of fixed mods: " + std::to_string(FixedModifications.size()));
			sb->appendLine("Number of variable mods: " + std::to_string(VariableModifications.size()));
			sb->appendLine("Dissociation Type: " + commonParameters->getDissociationType());

			sb->appendLine("protease: " + commonParameters->getDigestionParams()->Protease);
			sb->appendLine("initiatorMethionineBehavior: " + commonParameters->getDigestionParams()->InitiatorMethionineBehavior);
			sb->appendLine("maximumMissedCleavages: " + commonParameters->getDigestionParams()->MaxMissedCleavages);
			sb->appendLine("minPeptideLength: " + commonParameters->getDigestionParams()->MinPeptideLength);
			sb->appendLine("maxPeptideLength: " + commonParameters->getDigestionParams()->MaxPeptideLength);
			sb->appendLine("maximumVariableModificationIsoforms: " + commonParameters->getDigestionParams()->MaxModificationIsoforms);
			sb->appendLine("digestionTerminus: " + commonParameters->getDigestionParams()->FragmentationTerminus);
			sb->appendLine("maxModsForEachPeptide: " + commonParameters->getDigestionParams()->MaxModsForPeptide);
			sb->appendLine("cleavageSpecificity: " + commonParameters->getDigestionParams()->SearchModeType);
			sb->appendLine("specificProtease: " + commonParameters->getDigestionParams()->SpecificProtease);

			sb->append("Localizeable mods: " + ProteinList.Select([&] (std::any b)
			{
				b::OneBasedPossibleLocalizedModifications->Count;
			}).Sum());

			delete sb;
			return sb->toString();
		}

		MetaMorpheusEngineResults *IndexingEngine::RunSpecific()
		{
			double progress = 0;
			int oldPercentProgress = 0;

			// digest database
			std::vector<PeptideWithSetModifications*> globalPeptides;

			ParallelOptions *tempVar = new ParallelOptions();
			tempVar->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
			Parallel::ForEach(Partitioner::Create(0, ProteinList.size()), tempVar, [&] (range, loopState)
			{
				std::vector<PeptideWithSetModifications*> localPeptides;
            
				for (int i = range::Item1; i < range::Item2; i++)
				{
					// Stop loop if canceled
					if (GlobalVariables::getStopLoops())
					{
						loopState::Stop();
						return;
					}
            
            
					localPeptides.AddRange(ProteinList[i]->Digest(commonParameters->getDigestionParams(), FixedModifications, VariableModifications));
            
            
					progress++;
					auto percentProgress = static_cast<int>((progress / ProteinList.size()) * 100);
            
					if (percentProgress > oldPercentProgress)
					{
						oldPercentProgress = percentProgress;
						ProgressEventArgs tempVar2(percentProgress, "Digesting proteins...", nestedIds);
						ReportProgress(&tempVar2);
					}
				}
            
				{
						std::lock_guard<std::mutex> lock(globalPeptides);
					globalPeptides.AddRange(localPeptides);
				}
			});

			// sort peptides by mass
			auto peptidesSortedByMass = globalPeptides.OrderBy([&] (std::any p)
			{
				p::MonoisotopicMass;
			}).ToList();
			globalPeptides.clear();

			// create fragment index
			std::vector<std::vector<int>> fragmentIndex;

			try
			{
				fragmentIndex = std::vector<std::vector<int>>(static_cast<int>(std::ceil(MaxFragmentSize)) * FragmentBinsPerDalton + 1);
			}
			catch (const OutOfMemoryException &e1)
			{
				throw MetaMorpheusException("Max fragment mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
			}

			// populate fragment index
			progress = 0;
			oldPercentProgress = 0;
			for (int peptideId = 0; peptideId < peptidesSortedByMass.size(); peptideId++)
			{
				auto fragmentMasses = peptidesSortedByMass[peptideId].Fragment(commonParameters->getDissociationType(), commonParameters->getDigestionParams()->FragmentationTerminus)->Select([&] (std::any m)
				{
					m::NeutralMass;
				}).ToList();

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
					ProgressEventArgs tempVar3(percentProgress, "Fragmenting peptides...", nestedIds);
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
				catch (const OutOfMemoryException &e2)
				{
					throw MetaMorpheusException("Max precursor mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
				}
				progress = 0;
				oldPercentProgress = 0;
				ProgressEventArgs tempVar4(0, "Creating precursor index...", nestedIds);
				ReportProgress(&tempVar4);

				for (int i = 0; i < peptidesSortedByMass.size(); i++)
				{
					double mass = peptidesSortedByMass[i].MonoisotopicMass;
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
						ProgressEventArgs tempVar5(percentProgress, "Creating precursor index...", nestedIds);
						ReportProgress(&tempVar5);
					}
				}
			}

			return new IndexingResults(peptidesSortedByMass, fragmentIndex, precursorIndex, this);
		}
	}
}
