#include "ClassicSearchEngine.h"
#include "../PrecursorSearchModes/MassDiffAcceptor.h"
#include "../PeptideSpectralMatch.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "../ScanWithIndexAndNotchInfo.h"

using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
	namespace ClassicSearch
	{

		ClassicSearchEngine::ClassicSearchEngine(std::vector<PeptideSpectralMatch*> &globalPsms, std::vector<Ms2ScanWithSpecificMass*> &arrayOfSortedMS2Scans, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<Protein*> &proteinList, MassDiffAcceptor *searchMode, CommonParameters *commonParameters, std::vector<std::wstring> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds), SearchMode(searchMode), Proteins(proteinList), FixedModifications(fixedModifications), VariableModifications(variableModifications), PeptideSpectralMatches(globalPsms), ArrayOfSortedMS2Scans(arrayOfSortedMS2Scans), MyScanPrecursorMasses(arrayOfSortedMS2Scans.Select([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray())
		{
		}

		MetaMorpheusEngineResults *ClassicSearchEngine::RunSpecific()
		{
			Status(L"Getting ms2 scans...");

			double proteinsSearched = 0;
			int oldPercentProgress = 0;

			// one lock for each MS2 scan; a scan can only be accessed by one thread at a time
			auto myLocks = std::vector<std::any>(PeptideSpectralMatches.size());
			for (int i = 0; i < myLocks.size(); i++)
			{
				myLocks[i] = std::any();
			}

			Status(L"Performing classic search...");

			if (Proteins.Any())
			{
				ParallelOptions *tempVar = new ParallelOptions();
				tempVar->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
				Parallel::ForEach(Partitioner::Create(0, Proteins.size()), tempVar, [&] (partitionRange, loopState)
				{
					for (int i = partitionRange::Item1; i < partitionRange::Item2; i++)
					{
						// Stop loop if canceled
						if (GlobalVariables::getStopLoops())
						{
							loopState::Stop();
							return;
						}
            
						// digest each protein into peptides and search for each peptide in all spectra within precursor mass tolerance
						for (PeptideWithSetModifications *peptide : Proteins[i]->Digest(commonParameters->getDigestionParams(), FixedModifications, VariableModifications))
						{
							std::vector<Product*> peptideTheorProducts = peptide->Fragment(commonParameters->getDissociationType(), commonParameters->getDigestionParams()->FragmentationTerminus).ToList();
            
							for (auto scan : GetAcceptableScans(peptide->MonoisotopicMass, SearchMode))
							{
								std::vector<MatchedFragmentIon*> matchedIons = MatchFragmentIons(scan->TheScan, peptideTheorProducts, commonParameters);
            
								double thisScore = CalculatePeptideScore(scan->TheScan->getTheScan(), matchedIons, 0);
            
								bool meetsScoreCutoff = thisScore >= commonParameters->getScoreCutoff();
            
								// this is thread-safe because even if the score improves from another thread writing to this PSM,
								// the lock combined with AddOrReplace method will ensure thread safety
								if (meetsScoreCutoff || commonParameters->getCalculateEValue())
								{
									{
									// valid hit (met the cutoff score); lock the scan to prevent other threads from accessing it
											std::lock_guard<std::mutex> lock(myLocks[scan->ScanIndex]);
										bool scoreImprovement = PeptideSpectralMatches[scan->ScanIndex] == nullptr || (thisScore - PeptideSpectralMatches[scan->ScanIndex]->getRunnerUpScore()) > -PeptideSpectralMatch::ToleranceForScoreDifferentiation;
            
										if (scoreImprovement)
										{
											if (PeptideSpectralMatches[scan->ScanIndex] == nullptr)
											{
												PeptideSpectralMatches[scan->ScanIndex] = new PeptideSpectralMatch(peptide, scan->Notch, thisScore, scan->ScanIndex, scan->TheScan, commonParameters->getDigestionParams(), matchedIons);
											}
											else
											{
												PeptideSpectralMatches[scan->ScanIndex]->AddOrReplace(peptide, thisScore, scan->Notch, commonParameters->getReportAllAmbiguity(), matchedIons);
											}
										}
            
										if (commonParameters->getCalculateEValue())
										{
											PeptideSpectralMatches[scan->ScanIndex]->getAllScores().push_back(thisScore);
										}
									}
								}
							}
						}
            
						// report search progress (proteins searched so far out of total proteins in database)
						proteinsSearched++;
						auto percentProgress = static_cast<int>((proteinsSearched / Proteins.size()) * 100);
            
						if (percentProgress > oldPercentProgress)
						{
							oldPercentProgress = percentProgress;
							ProgressEventArgs tempVar2(percentProgress, L"Performing classic search... ", nestedIds);
							ReportProgress(&tempVar2);
						}
					}
				});
			}

			// remove peptides below the score cutoff that were stored to calculate expectation values
			if (commonParameters->getCalculateEValue())
			{
				for (int i = 0; i < PeptideSpectralMatches.size(); i++)
				{
					if (PeptideSpectralMatches[i] != nullptr && PeptideSpectralMatches[i]->getScore() < commonParameters->getScoreCutoff())
					{
						PeptideSpectralMatches[i] = nullptr;
					}
				}
			}

			for (PeptideSpectralMatch *psm : PeptideSpectralMatches.Where([&] (std::any p)
			{
				return p != nullptr;
			}))
			{
				psm::ResolveAllAmbiguities();
			}

			return new MetaMorpheusEngineResults(this);
		}

		std::vector<ScanWithIndexAndNotchInfo*> ClassicSearchEngine::GetAcceptableScans(double peptideMonoisotopicMass, MassDiffAcceptor *searchMode)
		{
			for (auto allowedIntervalWithNotch : searchMode->GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMonoisotopicMass).ToList())
			{
				DoubleRange *allowedInterval = allowedIntervalWithNotch->AllowedInterval;
				int scanIndex = GetFirstScanWithMassOverOrEqual(allowedInterval->Minimum);
				if (scanIndex < ArrayOfSortedMS2Scans.size())
				{
					auto scanMass = MyScanPrecursorMasses[scanIndex];
					while (scanMass <= allowedInterval->Maximum)
					{
						auto scan = ArrayOfSortedMS2Scans[scanIndex];
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
						yield return new ScanWithIndexAndNotchInfo(scan, allowedIntervalWithNotch->getNotch(), scanIndex);
						scanIndex++;
						if (scanIndex == ArrayOfSortedMS2Scans.size())
						{
							break;
						}

						scanMass = MyScanPrecursorMasses[scanIndex];
					}
				}
			}
		}

		int ClassicSearchEngine::GetFirstScanWithMassOverOrEqual(double minimum)
		{
			int index = Array::BinarySearch(MyScanPrecursorMasses, minimum);
			if (index < 0)
			{
				index = ~index;
			}

			// index of the first element that is larger than value
			return index;
		}
	}
}
