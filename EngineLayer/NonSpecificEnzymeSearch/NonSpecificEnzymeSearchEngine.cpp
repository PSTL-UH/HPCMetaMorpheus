#include "NonSpecificEnzymeSearchEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "../PrecursorSearchModes/MassDiffAcceptor.h"
#include "../MetaMorpheusEngineResults.h"
#include "../EventArgs/ProgressEventArgs.h"

using namespace Chemistry;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::ModernSearch;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
	namespace NonSpecificEnzymeSearch
	{

const double NonSpecificEnzymeSearchEngine::WaterMonoisotopicMass = PeriodicTable::GetElement(L"H").PrincipalIsotope::AtomicMass * 2 + PeriodicTable::GetElement(L"O").PrincipalIsotope::AtomicMass;

		NonSpecificEnzymeSearchEngine::NonSpecificEnzymeSearchEngine(std::vector<std::vector<PeptideSpectralMatch*>> &globalPsms, std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans, std::vector<PeptideWithSetModifications*> &peptideIndex, std::vector<std::vector<int>&> &fragmentIndex, std::vector<std::vector<int>&> &precursorIndex, int currentPartition, CommonParameters *CommonParameters, MassDiffAcceptor *massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled, std::vector<std::wstring> &nestedIds) : ModernSearchEngine(nullptr, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, CommonParameters, massDiffAcceptor, maximumMassThatFragmentIonScoreIsDoubled, nestedIds), PrecursorIndex(precursorIndex), MinimumPeptideLength(commonParameters->getDigestionParams()->MinPeptideLength)
		{
			GlobalCategorySpecificPsms = globalPsms;
			ModifiedParametersNoComp = commonParameters->CloneWithNewTerminus(, std::make_optional(false));
			ProductTypesToSearch = DissociationTypeCollection::ProductsFromDissociationType[commonParameters->getDissociationType()].Intersect(TerminusSpecificProductTypes::ProductIonTypesFromSpecifiedTerminus[commonParameters->getDigestionParams()->FragmentationTerminus])->ToList();
		}

		MetaMorpheusEngineResults *NonSpecificEnzymeSearchEngine::RunSpecific()
		{
			double progress = 0;
			int oldPercentProgress = 0;
			ProgressEventArgs tempVar(oldPercentProgress, L"Performing nonspecific search... " + std::to_wstring(CurrentPartition) + L"/" + std::to_wstring(commonParameters->getTotalPartitions()), nestedIds);
			ReportProgress(&tempVar);

			unsigned char byteScoreCutoff = static_cast<unsigned char>(commonParameters->getScoreCutoff());

			ParallelOptions *tempVar2 = new ParallelOptions();
			tempVar2->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
			Parallel::ForEach(Partitioner::Create(0, ListOfSortedMs2Scans.size()), tempVar2, [&] (std::any range)
			{
				std::vector<unsigned char> scoringTable(PeptideIndex.size());
				std::unordered_set<int> idsOfPeptidesPossiblyObserved;
            
				for (int i = range::Item1; i < range::Item2; i++)
				{
					// empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
					Array::Clear(scoringTable, 0, scoringTable.Length);
					idsOfPeptidesPossiblyObserved.Clear();
					Ms2ScanWithSpecificMass *scan = ListOfSortedMs2Scans[i];
            
					//get bins to add points to
					std::vector<int> allBinsToSearch = GetBinsToSearch(scan);
            
					//the entire indexed scoring is done here
					for (int j = 0; j < allBinsToSearch.Count; j++)
					{
						std::for_each(FragmentIndex[allBinsToSearch[j]].begin(), FragmentIndex[allBinsToSearch[j]].end(), [&] (std::any id)
						{
					scoringTable[id]++;
						});
					}
            
					//populate ids of possibly observed with those containing allowed precursor masses
					std::vector<AllowedIntervalWithNotch*> validIntervals = MassDiffAcceptor->GetAllowedPrecursorMassIntervalsFromObservedMass(scan->getPrecursorMass()).ToList(); //get all valid notches
					for (auto interval : validIntervals)
					{
						int obsPrecursorFloorMz = static_cast<int>(std::floor(interval->AllowedInterval->Minimum * FragmentBinsPerDalton));
						int obsPrecursorCeilingMz = static_cast<int>(std::ceil(interval->AllowedInterval->Maximum * FragmentBinsPerDalton));
            
						for (auto pt : ProductTypesToSearch)
						{
							int dissociationBinShift = static_cast<int>(BankersRounding::round((WaterMonoisotopicMass - DissociationTypeCollection::GetMassShiftFromProductType(pt)) * FragmentBinsPerDalton));
							int lowestBin = obsPrecursorFloorMz - dissociationBinShift;
							int highestBin = obsPrecursorCeilingMz - dissociationBinShift;
							for (int bin = lowestBin; bin <= highestBin; bin++)
							{
								if (bin < FragmentIndex.size() && FragmentIndex[bin].size() > 0)
								{
									std::for_each(FragmentIndex[bin].begin(), FragmentIndex[bin].end(), [&] (std::any id)
									{
					idsOfPeptidesPossiblyObserved.Add(id);
									});
								}
							}
						}
            
						for (int bin = obsPrecursorFloorMz; bin <= obsPrecursorCeilingMz; bin++) //no bin shift, since they're precursor masses
						{
							if (bin < PrecursorIndex.size() && PrecursorIndex[bin].size() > 0)
							{
								std::for_each(PrecursorIndex[bin].begin(), PrecursorIndex[bin].end(), [&] (std::any id)
								{
					idsOfPeptidesPossiblyObserved.Add(id);
								});
							}
						}
					}
            
					// done with initial scoring; refine scores and create PSMs
					if (idsOfPeptidesPossiblyObserved.Any())
					{
						int maxInitialScore = idsOfPeptidesPossiblyObserved.Max([&] (std::any id)
						{
					scoringTable[id];
						}) + 1;
						while (maxInitialScore > commonParameters->getScoreCutoff()) //go through all until we hit the end
						{
							maxInitialScore--;
							for (int id : idsOfPeptidesPossiblyObserved.Where([&] (std::any id)
							{
					return scoringTable[id] == maxInitialScore;
							}))
							{
								PeptideWithSetModifications *peptide = PeptideIndex[id];
								std::vector<Product*> peptideTheorProducts = peptide->Fragment(commonParameters->getDissociationType(), commonParameters->getDigestionParams()->FragmentationTerminus).ToList();
            
								std::tuple<int, PeptideWithSetModifications*> notchAndUpdatedPeptide = Accepts(peptideTheorProducts, scan->getPrecursorMass(), peptide, commonParameters->getDigestionParams()->FragmentationTerminus, MassDiffAcceptor);
								int notch = notchAndUpdatedPeptide.Item1;
								if (notch >= 0)
								{
									peptide = notchAndUpdatedPeptide.Item2;
									peptideTheorProducts = peptide->Fragment(commonParameters->getDissociationType(), FragmentationTerminus::Both).ToList();
									std::vector<MatchedFragmentIon*> matchedIons = MatchFragmentIons(scan, peptideTheorProducts, ModifiedParametersNoComp);
            
									double thisScore = CalculatePeptideScore(scan->getTheScan(), matchedIons, MaxMassThatFragmentIonScoreIsDoubled);
									if (thisScore > commonParameters->getScoreCutoff())
									{
										std::vector<PeptideSpectralMatch*> localPeptideSpectralMatches = GlobalCategorySpecificPsms[static_cast<int>(FdrClassifier::GetCleavageSpecificityCategory(peptide->CleavageSpecificityForFdrCategory))];
										if (localPeptideSpectralMatches[i] == nullptr)
										{
											localPeptideSpectralMatches[i] = new PeptideSpectralMatch(peptide, notch, thisScore, i, scan, commonParameters->getDigestionParams(), matchedIons);
										}
										else
										{
											localPeptideSpectralMatches[i]->AddOrReplace(peptide, thisScore, notch, commonParameters->getReportAllAmbiguity(), matchedIons);
										}
									}
								}
							}
						}
					}
					// report search progress
					progress++;
					int percentProgress = static_cast<int>((progress / ListOfSortedMs2Scans.size()) * 100);
            
					if (percentProgress > oldPercentProgress)
					{
						oldPercentProgress = percentProgress;
						ProgressEventArgs tempVar3(percentProgress, L"Performing nonspecific search... " + std::to_wstring(CurrentPartition) + L"/" + std::to_wstring(commonParameters->getTotalPartitions()), nestedIds);
						ReportProgress(&tempVar3);
					}
				}
			});
			return new MetaMorpheusEngineResults(this);
		}

		std::tuple<int, PeptideWithSetModifications*> NonSpecificEnzymeSearchEngine::Accepts(std::vector<Product*> &fragments, double scanPrecursorMass, PeptideWithSetModifications *peptide, FragmentationTerminus *fragmentationTerminus, MassDiffAcceptor *searchMode)
		{
			//all masses in N and CTerminalMasses are b-ion masses, which are one water away from a full peptide
			int localminPeptideLength = commonParameters->getDigestionParams()->MinPeptideLength;

			for (int i = localminPeptideLength - 1; i < fragments.size(); i++) //minus one start, because fragment 1 is at index 0
			{
				Product *fragment = fragments[i];
				double theoMass = fragment->NeutralMass - DissociationTypeCollection::GetMassShiftFromProductType(fragment->ProductType) + WaterMonoisotopicMass;
				int notch = searchMode->Accepts(scanPrecursorMass, theoMass);
				if (notch >= 0)
				{
					PeptideWithSetModifications *updatedPwsm = nullptr;
					if (fragmentationTerminus == FragmentationTerminus::N)
					{
						int endResidue = peptide->OneBasedStartResidueInProtein + fragment->TerminusFragment.FragmentNumber - 1; //-1 for one based index
						std::unordered_map<int, Modification*> updatedMods;
						for (auto mod : peptide->AllModsOneIsNterminus)
						{
							if (mod.first < endResidue - peptide->OneBasedStartResidueInProtein + 3) //check if we cleaved it off, +1 for N-terminus being mod 1 and first residue being mod 2, +1 again for the -1 on end residue for one based index, +1 (again) for the one-based start residue
							{
								updatedMods.emplace(mod.first, mod.second);
							}
						}
						updatedPwsm = new PeptideWithSetModifications(peptide->Protein, peptide->DigestionParams, peptide->OneBasedStartResidueInProtein, endResidue, CleavageSpecificity::Unknown, L"", 0, updatedMods, 0);
					}
					else
					{
						int startResidue = peptide->OneBasedEndResidueInProtein - fragment->TerminusFragment.FragmentNumber + 1; //plus one for one based index
						std::unordered_map<int, Modification*> updatedMods; //updateMods
						int indexShift = startResidue - peptide->OneBasedStartResidueInProtein;
						for (auto mod : peptide->AllModsOneIsNterminus)
						{
							if (mod.first > indexShift + 1) //check if we cleaved it off, +1 for N-terminus being mod 1 and first residue being 2
							{
								int key = mod.first - indexShift;
								updatedMods.emplace(key, mod.second);
							}
						}
						updatedPwsm = new PeptideWithSetModifications(peptide->Protein, peptide->DigestionParams, startResidue, peptide->OneBasedEndResidueInProtein, CleavageSpecificity::Unknown, L"", 0, updatedMods, 0);
					}

					delete updatedPwsm;
					return std::tuple<int, PeptideWithSetModifications*>(notch, updatedPwsm);
				}
				else if (theoMass > scanPrecursorMass)
				{
					break;
				}
			}
			//if the theoretical and experimental have the same mass
			if (fragments.size() > localminPeptideLength)
			{
				double totalMass = peptide->MonoisotopicMass; // + Constants.ProtonMass;
				int notch = searchMode->Accepts(scanPrecursorMass, totalMass);
				if (notch >= 0)
				{
					//need to update so that the cleavage specificity is recorded
					PeptideWithSetModifications *updatedPwsm = new PeptideWithSetModifications(peptide->Protein, peptide->DigestionParams, peptide->OneBasedStartResidueInProtein, peptide->OneBasedEndResidueInProtein, CleavageSpecificity::Unknown, L"", 0, peptide->AllModsOneIsNterminus, peptide->NumFixedMods);

					delete updatedPwsm;
					return std::tuple<int, PeptideWithSetModifications*>(notch, updatedPwsm);
				}
			}
			return std::tuple<int, PeptideWithSetModifications*>(-1, nullptr);
		}

		std::vector<PeptideSpectralMatch*> NonSpecificEnzymeSearchEngine::ResolveFdrCategorySpecificPsms(std::vector<std::vector<PeptideSpectralMatch*>&> &AllPsms, int numNotches, const std::wstring &taskId, CommonParameters *commonParameters)
		{
			//update all psms with peptide info
			AllPsms.ToList()->Where([&] (std::any psmArray)
			{
				return psmArray != nullptr;
			}).ToList()->ForEach([&] (std::any psmArray)
			{
				psmArray::Where([&] (std::any psm)
				{
					return psm != nullptr;
				}).ToList()->ForEach([&] (std::any psm)
				{
					psm::ResolveAllAmbiguities();
				});
			});

			for (auto psmsArray : AllPsms)
			{
				if (psmsArray.size() > 0)
				{
					std::vector<PeptideSpectralMatch*> cleanedPsmsArray = psmsArray.Where([&] (std::any b)
					{
						return b != nullptr;
					}).OrderByDescending([&] (std::any b)
					{
						b::Score;
					}).ThenBy([&] (std::any b)
					{
						b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
					}).GroupBy([&] (std::any b)
					{
						(b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass);
					})->Select([&] (std::any b)
					{
						b::First();
					}).ToList();

					FdrAnalysisEngine tempVar(cleanedPsmsArray, numNotches, commonParameters, new std::vector<std::wstring> {taskId});
					(&tempVar)->Run();

					for (int i = 0; i < psmsArray.size(); i++)
					{
						if (psmsArray[i] != nullptr)
						{
							if (psmsArray[i]->getFdrInfo() == nullptr) //if it was grouped in the cleanedPsmsArray
							{
								psmsArray[i] = nullptr;
							}
						}
					}
				}
			}

			std::vector<int> ranking(AllPsms.size()); //high int is good ranking
			std::vector<int> indexesOfInterest;
			for (int i = 0; i < ranking.size(); i++)
			{
				if (AllPsms[i].size() > 0)
				{
					ranking[i] = AllPsms[i].Where([&] (std::any x)
					{
						return x != nullptr;
					})->Count([&] (std::any x)
					{
						return x::FdrInfo::QValue <= 0.01;
					}); //set ranking as number of psms above 1% FDR
					indexesOfInterest.push_back(i);
				}
			}

			//get the index of the category with the highest ranking
			int majorCategoryIndex = indexesOfInterest[0];
			for (int i = 1; i < indexesOfInterest.size(); i++)
			{
				int currentCategoryIndex = indexesOfInterest[i];
				if (ranking[currentCategoryIndex] > ranking[majorCategoryIndex])
				{
					majorCategoryIndex = currentCategoryIndex;
				}
			}

			//update other category q-values
			//There's a chance of weird categories getting a random decoy before a random target, but we don't want to give that target a q value of zero.
			//We can't just take the q of the first decoy, because if the target wasn't random (score = 40), but there are no other targets before the decoy (score = 5), then we're incorrectly dinging the target
			//The current solution is such that if a minor category has a lower q value than it's corresponding score in the major category, then its q-value is changed to what it would be in the major category
			std::vector<PeptideSpectralMatch*> majorCategoryPsms = AllPsms[majorCategoryIndex].Where([&] (std::any x)
			{
				return x != nullptr;
			}).OrderByDescending([&] (std::any x)
			{
				x::Score;
			}).ToList(); //get sorted major category
			for (int i = 0; i < indexesOfInterest.size(); i++)
			{
				int minorCategoryIndex = indexesOfInterest[i];
				if (minorCategoryIndex != majorCategoryIndex)
				{
					std::vector<PeptideSpectralMatch*> minorCategoryPsms = AllPsms[minorCategoryIndex].Where([&] (std::any x)
					{
						return x != nullptr;
					}).OrderByDescending([&] (std::any x)
					{
						x::Score;
					}).ToList(); //get sorted minor category
					int minorPsmIndex = 0;
					int majorPsmIndex = 0;
					while (minorPsmIndex < minorCategoryPsms.size() && majorPsmIndex < majorCategoryPsms.size()) //while in the lists
					{
						PeptideSpectralMatch *majorPsm = majorCategoryPsms[majorPsmIndex];
						PeptideSpectralMatch *minorPsm = minorCategoryPsms[minorPsmIndex];
						//major needs to be a lower score than the minor
						if (majorPsm->getScore() > minorPsm->getScore())
						{
							majorPsmIndex++;
						}
						else
						{
							if (majorPsm->getFdrInfo()->getQValue() > minorPsm->getFdrInfo()->getQValue())
							{
								minorPsm->getFdrInfo()->setQValue(majorPsm->getFdrInfo()->getQValue());
							}
							minorPsmIndex++;
						}
					}
					//wrap up if we hit the end of the major category
					while (minorPsmIndex < minorCategoryPsms.size())
					{
						PeptideSpectralMatch *majorPsm = majorCategoryPsms[majorPsmIndex - 1]; //-1 because it's out of index right now
						PeptideSpectralMatch *minorPsm = minorCategoryPsms[minorPsmIndex];
						if (majorPsm->getFdrInfo()->getQValue() > minorPsm->getFdrInfo()->getQValue())
						{
							minorPsm->getFdrInfo()->setQValue(majorPsm->getFdrInfo()->getQValue());
						}
						minorPsmIndex++;
					}
				}
			}

			int numTotalSpectraWithPrecursors = AllPsms[indexesOfInterest[0]].size();
			std::vector<PeptideSpectralMatch*> bestPsmsList;
			for (int i = 0; i < numTotalSpectraWithPrecursors; i++)
			{
				PeptideSpectralMatch *bestPsm = nullptr;
				double lowestQ = std::numeric_limits<double>::max();
				int bestIndex = -1;
				for (auto index : indexesOfInterest) //foreach category
				{
					PeptideSpectralMatch *currentPsm = AllPsms[index][i];
					if (currentPsm != nullptr)
					{
						double currentQValue = currentPsm->getFdrInfo()->getQValue();
						if (currentQValue < lowestQ || (currentQValue == lowestQ && currentPsm->getScore() > bestPsm->getScore()))
						{
							if (bestIndex != -1)
							{
								//remove the old one so we don't use it for fdr later
								AllPsms[bestIndex][i] = nullptr;
							}
							bestPsm = currentPsm;
							lowestQ = currentQValue;
							bestIndex = index;
						}
						else //remove the old one so we don't use it for fdr later
						{
							AllPsms[index][i] = nullptr;
						}
					}
				}
				if (bestPsm != nullptr)
				{
					bestPsmsList.push_back(bestPsm);
				}
			}

			//It's probable that psms from some categories were removed by psms from other categories.
			//however, the fdr is still affected by their presence, since it was calculated before their removal.
			for (auto psmsArray : AllPsms)
			{
				if (psmsArray.size() > 0)
				{
					std::vector<PeptideSpectralMatch*> cleanedPsmsArray = psmsArray.Where([&] (std::any b)
					{
						return b != nullptr;
					}).OrderByDescending([&] (std::any b)
					{
						b::Score;
					}).ThenBy([&] (std::any b)
					{
						b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
					}).ToList();

					FdrAnalysisEngine tempVar2(cleanedPsmsArray, numNotches, commonParameters, new std::vector<std::wstring> {taskId});
					(&tempVar2)->Run();
				}
			}

			return bestPsmsList.OrderBy([&] (std::any b)
			{
				b::FdrInfo::QValue;
			}).ThenByDescending([&] (std::any b)
			{
				b::Score;
			}).ToList();
		}
	}
}
