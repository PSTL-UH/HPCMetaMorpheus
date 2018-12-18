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

		IndexingEngine::IndexingEngine(std::vector<Protein*> &proteinList, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, int currentPartition, DecoyType *decoyType, CommonParameters *commonParams, double maxFragmentSize, bool generatePrecursorIndex, std::vector<FileInfo*> &proteinDatabases, std::vector<std::wstring> &nestedIds) : MetaMorpheusEngine(commonParams, nestedIds), ProteinList(proteinList), FixedModifications(fixedModifications), VariableModifications(variableModifications), CurrentPartition(currentPartition + 1), DecoyType(decoyType), MaxFragmentSize(maxFragmentSize), GeneratePrecursorIndex(generatePrecursorIndex), ProteinDatabases(proteinDatabases)
		{
		}

		std::wstring IndexingEngine::ToString()
		{
			auto sb = new StringBuilder();
			sb->appendLine(L"Databases: " + std::wstring::Join(L",", ProteinDatabases.OrderBy([&] (std::any p)
			{
				p->Name;
			})->Select([&] (std::any p)
			{
			delete sb;
				return p->Name + L":" + p::CreationTime;
			})));
			sb->appendLine(L"Partitions: " + std::to_wstring(CurrentPartition) + L"/" + std::to_wstring(commonParameters->getTotalPartitions()));
			sb->appendLine(L"Precursor Index: " + StringHelper::toString(GeneratePrecursorIndex));
			sb->appendLine(L"Search Decoys: " + DecoyType);
			sb->appendLine(L"Number of proteins: " + std::to_wstring(ProteinList.size()));
			sb->appendLine(L"Number of fixed mods: " + std::to_wstring(FixedModifications.size()));
			sb->appendLine(L"Number of variable mods: " + std::to_wstring(VariableModifications.size()));
			sb->appendLine(L"Dissociation Type: " + commonParameters->getDissociationType());

			sb->appendLine(L"protease: " + commonParameters->getDigestionParams()->Protease);
			sb->appendLine(L"initiatorMethionineBehavior: " + commonParameters->getDigestionParams()->InitiatorMethionineBehavior);
			sb->appendLine(L"maximumMissedCleavages: " + commonParameters->getDigestionParams()->MaxMissedCleavages);
			sb->appendLine(L"minPeptideLength: " + commonParameters->getDigestionParams()->MinPeptideLength);
			sb->appendLine(L"maxPeptideLength: " + commonParameters->getDigestionParams()->MaxPeptideLength);
			sb->appendLine(L"maximumVariableModificationIsoforms: " + commonParameters->getDigestionParams()->MaxModificationIsoforms);
			sb->appendLine(L"digestionTerminus: " + commonParameters->getDigestionParams()->FragmentationTerminus);
			sb->appendLine(L"maxModsForEachPeptide: " + commonParameters->getDigestionParams()->MaxModsForPeptide);
			sb->appendLine(L"cleavageSpecificity: " + commonParameters->getDigestionParams()->SearchModeType);
			sb->appendLine(L"specificProtease: " + commonParameters->getDigestionParams()->SpecificProtease);

			sb->append(L"Localizeable mods: " + ProteinList.Select([&] (std::any b)
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
						ProgressEventArgs tempVar2(percentProgress, L"Digesting proteins...", nestedIds);
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
				throw MetaMorpheusException(L"Max fragment mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
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
					ProgressEventArgs tempVar3(percentProgress, L"Fragmenting peptides...", nestedIds);
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
					throw MetaMorpheusException(L"Max precursor mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
				}
				progress = 0;
				oldPercentProgress = 0;
				ProgressEventArgs tempVar4(0, L"Creating precursor index...", nestedIds);
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
						ProgressEventArgs tempVar5(percentProgress, L"Creating precursor index...", nestedIds);
						ReportProgress(&tempVar5);
					}
				}
			}

			return new IndexingResults(peptidesSortedByMass, fragmentIndex, precursorIndex, this);
		}
	}
}
