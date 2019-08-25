#include "FdrAnalysisEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "FdrAnalysisResults.h"
#include "../GlobalVariables.h"

using namespace MathNet::Numerics;
namespace EngineLayer
{
	namespace FdrAnalysis
	{

		FdrAnalysisEngine::FdrAnalysisEngine(std::vector<PeptideSpectralMatch*> &psms, int massDiffAcceptorNumNotches, CommonParameters *commonParameters, std::vector<std::string> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds), MassDiffAcceptorNumNotches(massDiffAcceptorNumNotches), UseDeltaScore(commonParameters->getUseDeltaScore()), CalculateEValue(commonParameters->getCalculateEValue()), ScoreCutoff(commonParameters->getScoreCutoff())
		{
			AllPsms = psms;
		}

		MetaMorpheusEngineResults *FdrAnalysisEngine::RunSpecific()
		{
			FdrAnalysisResults *myAnalysisResults = new FdrAnalysisResults(this);

			Status("Running FDR analysis...");
			DoFalseDiscoveryRateAnalysis(myAnalysisResults);

			myAnalysisResults->setPsmsWithin1PercentFdr(AllPsms.size()([&] (std::any b)
			{
//C# TO C++ CONVERTER TODO TASK: A 'delete myAnalysisResults' statement was not added since myAnalysisResults was passed to a method or constructor. Handle memory management manually.
				return b::FdrInfo::QValue < 0.01;
			}));

//C# TO C++ CONVERTER TODO TASK: A 'delete myAnalysisResults' statement was not added since myAnalysisResults was passed to a method or constructor. Handle memory management manually.
			return myAnalysisResults;
		}

		void FdrAnalysisEngine::DoFalseDiscoveryRateAnalysis(FdrAnalysisResults *myAnalysisResults)
		{
			// Stop if canceled
			if (GlobalVariables::getStopLoops())
			{
				return;
			}

			// calculate FDR on a per-protease basis (targets and decoys for a specific protease)
			auto psmsGroupedByProtease = AllPsms.GroupBy([&] (std::any p)
			{
				p::DigestionParams::Protease;
			});

			for (auto proteasePsms : psmsGroupedByProtease)
			{
				auto psms = proteasePsms->ToList();

				// generate the null distribution for e-value calculations
				double globalMeanScore = 0;
				int globalMeanCount = 0;

				if (CalculateEValue && psms.Any())
				{
					std::vector<double> combinedScores;

					for (auto psm : psms)
					{
						std::sort(psm->AllScores.begin(), psm->AllScores.end());
						combinedScores.insert(combinedScores.end(), psm->getAllScores().begin(), psm->getAllScores().end());

						//remove top scoring peptide
						if (combinedScores.Any())
						{
							combinedScores.pop_back();
						}
					}

					if (combinedScores.Any())
					{
						globalMeanScore = combinedScores.Average();
						globalMeanCount = static_cast<int>(static_cast<double>(combinedScores.size()) / psms.size());
					}
					else
					{
						// should be a very rare case... if there are PSMs but each PSM only has one hit
						globalMeanScore = 0;
						globalMeanCount = 0;
					}
				}

				//Calculate delta scores for the psms (regardless of if we are using them)
				for (auto psm : psms)
				{
					if (psm != nullptr)
					{
						psm->CalculateDeltaScore(ScoreCutoff);
					}
				}

				//determine if Score or DeltaScore performs better
				if (UseDeltaScore)
				{
					constexpr double qValueCutoff = 0.01; //optimize to get the most PSMs at a 1% FDR

					std::vector<PeptideSpectralMatch*> scoreSorted = psms.OrderByDescending([&] (std::any b)
					{
						b::Score;
					}).ThenBy([&] (std::any b)
					{
						b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
					}).GroupBy([&] (std::any b)
					{
						return std::tuple < std::string;
					}, int, std::optional<double>>(b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass))->Select([&] (std::any b)
					{
						b::First();
					}).ToList();
					int ScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);
					scoreSorted = psms.OrderByDescending([&] (std::any b)
					{
						b::DeltaScore;
					}).ThenBy([&] (std::any b)
					{
						b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
					}).GroupBy([&] (std::any b)
					{
						return std::tuple < std::string;
					}, int, std::optional<double>>(b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass))->Select([&] (std::any b)
					{
						b::First();
					}).ToList();
					int DeltaScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);

					//sort by best method
					myAnalysisResults->setDeltaScoreImprovement(DeltaScorePSMs > ScorePSMs);
					psms = myAnalysisResults->getDeltaScoreImprovement() ? psms.OrderByDescending([&] (std::any b)
					{
						b::DeltaScore;
					}).ThenBy([&] (std::any b)
					{
						b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
					}).ToList() : psms.OrderByDescending([&] (std::any b)
					{
						b::Score;
					}).ThenBy([&] (std::any b)
					{
						b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
					}).ToList();
				}
				else //sort by score
				{
					psms = psms.OrderByDescending([&] (std::any b)
					{
						b::Score;
					}).ThenBy([&] (std::any b)
					{
						b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
					}).ToList();
				}

				double cumulativeTarget = 0;
				double cumulativeDecoy = 0;

				//set up arrays for local FDRs
				std::vector<double> cumulativeTargetPerNotch(MassDiffAcceptorNumNotches + 1);
				std::vector<double> cumulativeDecoyPerNotch(MassDiffAcceptorNumNotches + 1);

				//Assign FDR values to PSMs
				for (int i = 0; i < psms.size(); i++)
				{
					// Stop if canceled
					if (GlobalVariables::getStopLoops())
					{
						break;
					}

					PeptideSpectralMatch *psm = psms[i];
					std::optional<int> tempVar = psm.getNotch();
					int notch = tempVar ? tempVar : MassDiffAcceptorNumNotches;
					if (psm->getIsDecoy())
					{
						// the PSM can be ambiguous between a target and a decoy sequence
						// in that case, count it as the fraction of decoy hits
						// e.g. if the PSM matched to 1 target and 2 decoys, it counts as 2/3 decoy
						double decoyHits = 0;
						double totalHits = 0;
						auto hits = psm->BestMatchingPeptides.GroupBy([&] (std::any p)
						{
							p::Peptide::FullSequence;
						});
						for (auto hit : hits)
						{
							if (hit->First().Peptide.Protein.IsDecoy)
							{
								decoyHits++;
							}
							totalHits++;
						}

						cumulativeDecoy += decoyHits / totalHits;
						cumulativeDecoyPerNotch[notch] += decoyHits / totalHits;
					}
					else
					{
						cumulativeTarget++;
						cumulativeTargetPerNotch[notch]++;
					}

					double qValue = std::min(1, cumulativeDecoy / cumulativeTarget);
					double qValueNotch = std::min(1, cumulativeDecoyPerNotch[notch] / cumulativeTargetPerNotch[notch]);

					double maximumLikelihood = 0;
					double eValue = 0;
					double eScore = 0;
					if (CalculateEValue)
					{
						eValue = GetEValue(psm, globalMeanCount, globalMeanScore, maximumLikelihood);
						eScore = -Math::Log(eValue, 10);
					}

					psm->SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetPerNotch[notch], cumulativeDecoyPerNotch[notch], qValueNotch, maximumLikelihood, eValue, eScore, CalculateEValue);
				}

				// set q-value thresholds such that a lower scoring PSM can't have 
				// a higher confidence than a higher scoring PSM
				//Populate min qValues
				double qValueThreshold = 1.0;
				std::vector<double> qValueNotchThreshold(MassDiffAcceptorNumNotches + 1);
				for (int i = 0; i < qValueNotchThreshold.size(); i++)
				{
					qValueNotchThreshold[i] = 1.0;
				}

				for (int i = psms.size() - 1; i >= 0; i--)
				{
					PeptideSpectralMatch *psm = psms[i];

					// threshold q-values
					if (psm->getFdrInfo()->getQValue() > qValueThreshold)
					{
						psm->getFdrInfo()->setQValue(qValueThreshold);
					}
					else if (psm->getFdrInfo()->getQValue() < qValueThreshold)
					{
						qValueThreshold = psm->getFdrInfo()->getQValue();
					}

					// threshold notch q-values
					std::optional<int> tempVar2 = psm.getNotch();
					int notch = tempVar2 ? tempVar2 : MassDiffAcceptorNumNotches;
					if (psm->getFdrInfo()->getQValueNotch() > qValueNotchThreshold[notch])
					{
						psm->getFdrInfo()->setQValueNotch(qValueNotchThreshold[notch]);
					}
					else if (psm->getFdrInfo()->getQValueNotch() < qValueNotchThreshold[notch])
					{
						qValueNotchThreshold[notch] = psm->getFdrInfo()->getQValueNotch();
					}
				}
			}
		}

		double FdrAnalysisEngine::GetEValue(PeptideSpectralMatch *psm, int globalMeanCount, double globalMeanScore, double &maximumLikelihood)
		{
			// get all of the PSM's scores for all hits, sort them, then remove the last value (the best score)
			std::vector<double> scoresWithoutBestHit;
			scoresWithoutBestHit.insert(scoresWithoutBestHit.end(), psm->getAllScores().begin(), psm->getAllScores().end());
			std::sort(scoresWithoutBestHit.begin(), scoresWithoutBestHit.end());

			if (scoresWithoutBestHit.Any())
			{
				scoresWithoutBestHit.pop_back();
			}

			// this is the "default" case for when there are no scores except the best hit
			// it uses a global mean score (all scores across all PSMs) to generate the null Poisson distribution
			// this will be overriden by the next few lines if there are enough scores in this PSM to estimate a null distribution
			double preValue = SpecialFunctions::GammaLowerRegularized(globalMeanScore, psm->getScore());
			maximumLikelihood = globalMeanScore;

			// calculate single-spectrum evalue if there are enough hits besides the best scoring peptide
			if (psm->getScore() == 0)
			{
				preValue = 1;
				maximumLikelihood = 0;
			}
			else if (scoresWithoutBestHit.Any())
			{
				maximumLikelihood = scoresWithoutBestHit.Average();

				// this is the cumulative distribution for the poisson at each score up to but not including the score of the winner.
				// This is the probability that the winner has of getting that score at random by matching against a SINGLE spectrum
				if (maximumLikelihood > 0)
				{
					preValue = SpecialFunctions::GammaLowerRegularized(maximumLikelihood, psm->getScore());
				}
			}

			// Now the probability of getting the winner's score goes up for each spectrum match.
			// We multiply the preValue by the number of theoretical spectrum within the tolerance to get this new probability.
			int count = scoresWithoutBestHit.size();
			if (count == 0)
			{
				count = globalMeanCount;
			}

			double probabilityOfScore = 1 - std::pow(preValue, count);
			return count * probabilityOfScore;
		}

		int FdrAnalysisEngine::GetNumPSMsAtqValueCutoff(std::vector<PeptideSpectralMatch*> &psms, double qValueCutoff)
		{
			int cumulative_target = 0;
			int cumulative_decoy = 0;
			for (auto psm : psms)
			{
				if (psm->getIsDecoy())
				{
					cumulative_decoy++;
					if (static_cast<double>(cumulative_decoy) / cumulative_target >= qValueCutoff)
					{
						return cumulative_target;
					}
				}
				else
				{
					cumulative_target++;
				}
			}
			return cumulative_target;
		}
	}
}
