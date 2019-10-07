#include "DataPointAcquisitionEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "LabeledDataPoint.h"
#include "DataPointAquisitionResults.h"

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics::AminoAcidPolymer;
namespace EngineLayer
{
	namespace Calibration
	{

		DataPointAcquisitionEngine::DataPointAcquisitionEngine(std::vector<PeptideSpectralMatch*> &goodIdentifications, MsDataFile *myMsDataFile, Tolerance *mzToleranceForMs1Search, int minMS1isotopicPeaksNeededForConfirmedIdentification, CommonParameters *commonParameters, std::vector<std::string> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds), GoodIdentifications(goodIdentifications), MyMsDataFile(myMsDataFile), MzToleranceForMs1Search(mzToleranceForMs1Search), MinMS1isotopicPeaksNeededForConfirmedIdentification(minMS1isotopicPeaksNeededForConfirmedIdentification)
		{
		}

		MetaMorpheusEngineResults *DataPointAcquisitionEngine::RunSpecific()
		{
			Status("Extracting data points:");
			// The final training point list

			int numMs1MassChargeCombinationsConsidered = 0;
			int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
			int numMs2MassChargeCombinationsConsidered = 0;
			int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
			std::vector<LabeledDataPoint*> Ms1List;
			std::vector<LabeledDataPoint*> Ms2List;

			int numIdentifications = GoodIdentifications.size();

			std::any lockObj = std::any();
			std::any lockObj2 = std::any();
			ParallelOptions *tempVar = new ParallelOptions();
			tempVar->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
			Parallel::ForEach(Partitioner::Create(0, numIdentifications), tempVar, [&] (fff, loopState)
			{
				for (int matchIndex = fff::Item1; matchIndex < fff::Item2; matchIndex++)
				{
					// Stop loop if canceled
					if (GlobalVariables::getStopLoops())
					{
						loopState::Stop();
						return;
					}
            
					PeptideSpectralMatch *identification = GoodIdentifications[matchIndex];
            
					// Each identification has an MS2 spectrum attached to it.
					int ms2scanNumber = identification->getScanNumber();
					int peptideCharge = identification->getScanPrecursorCharge();
					if (identification->getFullSequence() == "")
					{
						continue;
					}
            
					auto representativeSinglePeptide = identification->BestMatchingPeptides.First().Peptide;
            
					// Get the peptide, don't forget to add the modifications!!!!
					auto SequenceWithChemicalFormulas = representativeSinglePeptide->SequenceWithChemicalFormulas;
					if (SequenceWithChemicalFormulas == nullptr || representativeSinglePeptide->AllModsOneIsNterminus.Any([&] (std::any b)
					{
					return b->Value->NeutralLosses != nullptr;
					}))
					{
						continue;
					}
            
					Peptide *coolPeptide = new Peptide(SequenceWithChemicalFormulas);
            
					auto ms2tuple = SearchMS2Spectrum(MyMsDataFile->GetOneBasedScan(ms2scanNumber), identification);
            
					{
							std::lock_guard<std::mutex> lock(lockObj2);
						Ms2List.AddRange(ms2tuple);
					}
            
					// Calculate isotopic distribution of the full peptide
					auto dist = IsotopicDistribution::GetDistribution(coolPeptide->GetChemicalFormula(), FineResolutionForIsotopeDistCalculation, 0.001);
            
					std::vector<double> theoreticalMasses = dist->Masses->ToArray();
					std::vector<double> theoreticalIntensities = dist->Intensities->ToArray();
            
					Array::Sort(theoreticalIntensities, theoreticalMasses, Comparer<double>::Create([&] (x, y)
					{
					y->CompareTo(x);
					}));
            
					auto ms1tupleBack = SearchMS1Spectra(theoreticalMasses, theoreticalIntensities, ms2scanNumber, -1, peptideCharge, identification);
            
					auto ms1tupleForward = SearchMS1Spectra(theoreticalMasses, theoreticalIntensities, ms2scanNumber, 1, peptideCharge, identification);
            
					{
							std::lock_guard<std::mutex> lock(lockObj);
						Ms1List.AddRange(ms1tupleBack->Item1);
						numMs1MassChargeCombinationsConsidered += ms1tupleBack->Item2;
						numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += ms1tupleBack->Item3;
						Ms1List.AddRange(ms1tupleForward->Item1);
						numMs1MassChargeCombinationsConsidered += ms1tupleForward->Item2;
						numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += ms1tupleForward->Item3;
					}
				}
			});

			// datapoints are ordered because they were acquired in a parallized search and we want repeatable results
			return new DataPointAquisitionResults(this, GoodIdentifications, Ms1List.OrderBy([&] (std::any p)
			{
				p::ScanNumber;
			}).ThenBy([&] (std::any p)
			{
				p::ExperimentalMz;
			}).ToList(), Ms2List.OrderBy([&] (std::any p)
			{
				p::ScanNumber;
			}).ThenBy([&] (std::any p)
			{
				p::ExperimentalMz;
			}).ToList(), numMs1MassChargeCombinationsConsidered, numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks, numMs2MassChargeCombinationsConsidered, numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
		}

		std::vector<LabeledDataPoint*> DataPointAcquisitionEngine::SearchMS2Spectrum(MsDataScan *ms2DataScan, PeptideSpectralMatch *identification)
		{
			std::vector<LabeledDataPoint*> result;

			if (ms2DataScan->MassSpectrum->Size == 0)
			{
				return result;
			}

			for (auto matchedIon : identification->MatchedFragmentIons)
			{
				double exptPeakMz = matchedIon->Mz;
				double exptPeakIntensity = matchedIon->Intensity;
				double injTime = (ms2DataScan->InjectionTime != nullptr) ? ms2DataScan->InjectionTime : NAN;

				LabeledDataPoint tempVar(exptPeakMz, ms2DataScan->OneBasedScanNumber, std::log(ms2DataScan->TotalIonCurrent), std::log(injTime), std::log(exptPeakIntensity), matchedIon->NeutralTheoreticalProduct.NeutralMass.ToMz(matchedIon->Charge), identification);
				result.push_back(&tempVar);
			}
			return result;
		}
	}

    //private (List<LabeledDataPoint>, int, int) SearchMS1Spectra(double[] theoreticalMasses, double[] theoreticalIntensities,
    //         int ms2spectrumIndex, int direction, int peptideCharge, PeptideSpectralMatch identification)
    //		{
    //			List<LabeledDataPoint> result = new List<LabeledDataPoint>();
    //			int numMs1MassChargeCombinationsConsidered = 0;
    //			int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
    //
    //			int theIndex;
    //			theIndex = direction == 1 ? ms2spectrumIndex : ms2spectrumIndex - 1;
    //
    //			bool addedAscan = true;
    //
    //			int highestKnownChargeForThisPeptide = peptideCharge;
    //			while (theIndex >= 1 && theIndex <= MyMsDataFile.NumSpectra && addedAscan)
    //			{
    //				int countForThisScan = 0;
    //				if (MyMsDataFile.GetOneBasedScan(theIndex).MsnOrder > 1)
    //				{
    //					theIndex += direction;
    //					continue;
    //				}
    //				addedAscan = false;
    //				var fullMS1scan = MyMsDataFile.GetOneBasedScan(theIndex);
    //				var scanWindowRange = fullMS1scan.ScanWindowRange;
    //				var fullMS1spectrum = fullMS1scan.MassSpectrum;
    //				if (fullMS1spectrum.Size == 0)
    //					break;
    //
    //				bool startingToAddCharges = false;
    //				int chargeToLookAt = 1;
    //				do
    //				{
    //					if (theoreticalMasses[0].ToMz(chargeToLookAt) > scanWindowRange.Maximum)
    //					{
    //						chargeToLookAt++;
    //						continue;
    //					}
    //					if (theoreticalMasses[0].ToMz(chargeToLookAt) < scanWindowRange.Minimum)
    //						break;
    //					var trainingPointsToAverage = new List<LabeledDataPoint>();
    //					foreach (double a in theoreticalMasses)
    //					{
    //						double theMZ = a.ToMz(chargeToLookAt);
    //
    //						var npwr = fullMS1spectrum.NumPeaksWithinRange(MzToleranceForMs1Search.GetMinimumValue(theMZ), MzToleranceForMs1Search.GetMaximumValue(theMZ));
    //						if (npwr == 0)
    //						{
    //							break;
    //						}
    //						numMs1MassChargeCombinationsConsidered++;
    //						if (npwr > 1)
    //						{
    //							numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks++;
    //							continue;
    //						}
    //
    //						var closestPeakIndex = fullMS1spectrum.GetClosestPeakIndex(theMZ);
    //						var closestPeakMZ = fullMS1spectrum.XArray[closestPeakIndex.Value];
    //
    //						highestKnownChargeForThisPeptide = Math.Max(highestKnownChargeForThisPeptide, chargeToLookAt);
    //						trainingPointsToAverage.Add(new LabeledDataPoint(closestPeakMZ, -1, double.NaN, double.NaN, Math.Log(fullMS1spectrum.YArray[closestPeakIndex.Value]), theMZ, nullptr));
    //					}
    //					// If started adding and suddnely stopped, go to next one, no need to look at higher charges
    //					if (trainingPointsToAverage.Count == 0 && startingToAddCharges)
    //					{
    //						break;
    //					}
    //					if ((trainingPointsToAverage.Count == 0 || (trainingPointsToAverage.Count == 1 && theoreticalIntensities[0] < 0.65)) && (peptideCharge <= chargeToLookAt))
    //					{
    //						break;
    //					}
    //					if ((trainingPointsToAverage.Count == 1 && theoreticalIntensities[0] < 0.65) || trainingPointsToAverage.Count < Math.Min(MinMS1isotopicPeaksNeededForConfirmedIdentification, theoreticalIntensities.Count()))
    //					{
    //					}
    //					else
    //					{
    //						addedAscan = true;
    //						startingToAddCharges = true;
    //						countForThisScan++;
    //						result.Add(new LabeledDataPoint(trainingPointsToAverage.Select(b => b.ExperimentalMz).Average(), fullMS1scan.OneBasedScanNumber, Math.Log(fullMS1scan.TotalIonCurrent), fullMS1scan.InjectionTime.HasValue ? Math.Log(fullMS1scan.InjectionTime.Value) : double.NaN, trainingPointsToAverage.Select(b => b.LogIntensity).Average(), trainingPointsToAverage.Select(b => b.TheoreticalMz).Average(), identification));
    //					}
    //					chargeToLookAt++;
    //				} while (chargeToLookAt <= highestKnownChargeForThisPeptide + 1);
    //				theIndex += direction;
    //			}
    //			return (result, numMs1MassChargeCombinationsConsidered, numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
    //		}
    
}
