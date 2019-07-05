#include "CalibrationEngine.h"
#include "DataPointAquisitionResults.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "LabeledDataPoint.h"

using namespace MassSpectrometry;
namespace EngineLayer
{
	namespace Calibration
	{

		MsDataFile *CalibrationEngine::getCalibratedDataFile() const
		{
			return privateCalibratedDataFile;
		}

		void CalibrationEngine::setCalibratedDataFile(MsDataFile *value)
		{
			privateCalibratedDataFile = value;
		}

		CalibrationEngine::CalibrationEngine(MsDataFile *myMSDataFile, DataPointAquisitionResults *datapoints, CommonParameters *commonParameters, std::vector<std::string> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds), MyMsDataFile(myMSDataFile), Datapoints(datapoints)
		{
		}

		MetaMorpheusEngineResults *CalibrationEngine::RunSpecific()
		{
			Status("Calibrating spectra");
			std::vector<LabeledDataPoint*> &ms1Points = Datapoints->getMs1List();
			std::vector<LabeledDataPoint*> &ms2Points = Datapoints->getMs2List();
			std::vector<MsDataScan*> originalScans = MyMsDataFile->GetAllScansList();

			std::vector<MsDataScan*> ms1Scans;
			std::vector<MsDataScan*> ms2Scans;
			//separate scans by msnOrder, because the mass accuracy varies between the two
			for (auto scan : originalScans)
			{
				if (scan->MsnOrder == 1)
				{
					ms1Scans.push_back(scan);
				}
				else
				{
					ms2Scans.push_back(scan);
				}
			}

			//create a way to go from the scan number to the scan order. This can be a single array, since there shouldn't be any overlap of scan numbers between ms1 and ms2 scans
			std::vector<int> scanNumberToScanPlacement = std::vector<int>(originalScans.Max([&] (std::any x)
			{
				x::OneBasedScanNumber;
			}) + 1);
			for (int i = 0; i < ms1Scans.size(); i++)
			{
				scanNumberToScanPlacement[ms1Scans[i]->OneBasedScanNumber] = i;
			}
			for (int i = 0; i < ms2Scans.size(); i++)
			{
				scanNumberToScanPlacement[ms2Scans[i]->OneBasedScanNumber] = i;
			}

			//Populate the weighted average relative error for each scan, where index of returned array is the placement
			std::vector<double> ms1RelativeErrors = PopulateErrors(ms1Points, scanNumberToScanPlacement, ms1Scans.size());
			std::vector<double> ms2RelativeErrors = PopulateErrors(ms2Points, scanNumberToScanPlacement, ms2Scans.size());

			//generate new scans
			std::vector<MsDataScan*> calibratedScans(originalScans.size());

			//hard copy original scans
			for (int i = 0; i < originalScans.size(); i++)
			{
				calibratedScans[i] = originalScans[i];
			}

			//apply a smoothing function, so that outlier scans aren't wildly shifted
			std::vector<double> ms1SmoothedErrors = SmoothErrors(ms1RelativeErrors);
			std::vector<double> ms2SmoothedErrors = SmoothErrors(ms2RelativeErrors);

			//calibrate the data
			int ms1Index = 0;
			int ms2Index = 0;
			double mostRecentMS1SmoothedError = ms1SmoothedErrors.FirstOrDefault(); //this is needed to update the precursor mass error of MS2 scans
			for (int scanIndex = 0; scanIndex < calibratedScans.size(); scanIndex++) //foreach scan
			{
				MsDataScan *originalScan = originalScans[scanIndex]; //get original scan
				if (originalScan->MsnOrder == 1) //if ms1
				{
					mostRecentMS1SmoothedError = ms1SmoothedErrors[ms1Index]; //update the mass error
					calibratedScans[scanIndex] = CalibrateScan(originalScan, mostRecentMS1SmoothedError);
					ms1Index++;
				}
				else //if ms2
				{
					calibratedScans[scanIndex] = CalibrateScan(originalScan, ms2SmoothedErrors[ms2Index], std::make_optional(mostRecentMS1SmoothedError));
					ms2Index++;
				}
			}

			MsDataFile tempVar(calibratedScans, MyMsDataFile->SourceFile);
			setCalibratedDataFile(&tempVar);
			return new MetaMorpheusEngineResults(this);
		}

		std::vector<double> CalibrationEngine::PopulateErrors(std::vector<LabeledDataPoint*> &datapoints, std::vector<int> &scanNumberToScanPlacement, int arrayLength)
		{
			//return an array of weighted average relative errors for each scan
			std::vector<double> averageRelativeErrors(arrayLength);

			int currentScanNumber = datapoints.front().ScanNumber; //get the first scan number. This should be an ordered list
			std::vector<(double massError, double logIntensity)*> localRelativeErrors; //create a list to store the mass error

			for (auto datapoint : datapoints)
			{
				if (datapoint->ScanNumber == currentScanNumber)
				{
					localRelativeErrors.push_back((datapoint->RelativeMzError, datapoint->LogIntensity));
				}
				else
				{
					//wrap up the old set
					averageRelativeErrors[scanNumberToScanPlacement[currentScanNumber]] = CalculateAverageRelativeErrors(localRelativeErrors);

					//update
					currentScanNumber = datapoint->ScanNumber;
					localRelativeErrors = {(datapoint->RelativeMzError, datapoint->LogIntensity)};
				}
			}
			//finish loop
			averageRelativeErrors[scanNumberToScanPlacement[currentScanNumber]] = CalculateAverageRelativeErrors(localRelativeErrors);

			return averageRelativeErrors;
		}

		double CalibrationEngine::CalculateAverageRelativeErrors(std::vector<(double massError, double logIntensity)*> &localRelativeErrors)
		{
			//double logIntensityToSubtract = localRelativeErrors.Min(x => x.logIntensity) - 1; // normalize each log intensity so that the minimum log intensity is 1. 
			//Convert from log to actual intensity to more heavily weight intensities.
			double weightedSumOfErrors = localRelativeErrors.Sum([&] (std::any x)
			{
				return x::massError * std::pow(10, x::logIntensity);
			});
			double sumOfIntensities = localRelativeErrors.Sum([&] (std::any x)
			{
				std::pow(10, x::logIntensity);
			});
			return weightedSumOfErrors / sumOfIntensities;
		}

		std::vector<double> CalibrationEngine::SmoothErrors(std::vector<double> &relativeErrors)
		{
			//impute missing values
			//not all scans are guarenteed to contain data points. We can infer these data point with nearby points.
			for (int index = 0; index < relativeErrors.size(); index++)
			{
				if (relativeErrors[index] == 0) //if there were no points, then the value should be perfectly zero (the double default)
				{
					int startingBlankIndex = index;
					//increase the index until we find the next scan containing a data point
					while (index < relativeErrors.size() && relativeErrors[index] == 0)
					{
						index++;
					}
					double nextError = index == relativeErrors.size() ? relativeErrors[startingBlankIndex - 1] : relativeErrors[index]; //can't go all the way through without any data points, the original function checks for enough data points (where enough is more than zero)
					double previousError = startingBlankIndex > 0 ? relativeErrors[startingBlankIndex - 1] : nextError;
					int numberOfConsecutiveScansWithoutDataPoints = index - startingBlankIndex;
					for (int tempIndex = 0; tempIndex < numberOfConsecutiveScansWithoutDataPoints; tempIndex++)
					{
						relativeErrors[startingBlankIndex + tempIndex] = ((tempIndex + 1) * nextError + (numberOfConsecutiveScansWithoutDataPoints - tempIndex - 1) * previousError) / numberOfConsecutiveScansWithoutDataPoints;
					}
				}
			}

			//get correction factor for each scan, where half of the smooth comes from both directions
			std::vector<double> smoothedErrors(relativeErrors.size());
			int leftIndex = 0; //starting scan index used for numbers less than the current scan
			int rightIndex = 1;
			double smoothedCorrectionFactor = 0; //this variable is the sum of all nearby errors, to be later divided by the number of summed errors for a smoothed average error
			//for scan #1 (index 0)
			//no left scans, because we're at the beginning of the file. Just populate the first "numberOfScansUsedForSmoothingOnEachSide" scan errors on the right side

			while (rightIndex < NumberOfScansUsedForSmoothingOnEachSide && rightIndex < relativeErrors.size())
			{
				smoothedCorrectionFactor += relativeErrors[rightIndex];
				rightIndex++;
			}

			//for each scan
			for (int i = 0; i < relativeErrors.size(); i++)
			{
				//for left index, remove an error if 
				if (i > NumberOfScansUsedForSmoothingOnEachSide)
				{
					smoothedCorrectionFactor -= relativeErrors[leftIndex];
					leftIndex++;
				}

				if (rightIndex < relativeErrors.size()) //need to check so we don't run off the end of the file
				{
					smoothedCorrectionFactor += relativeErrors[rightIndex];
					rightIndex++;
				}

				smoothedErrors[i] = smoothedCorrectionFactor / (rightIndex - leftIndex);
			}
			return smoothedErrors;
		}

		MsDataScan *CalibrationEngine::CalibrateScan(MsDataScan *oldScan, double smoothedRelativeError, std::optional<double> &precursorSmoothedRelativeError)
		{
			double correctionFactor = 1 - smoothedRelativeError; //create the multiplier. Positive mass errors mean that the experimental mass was greater than the theoretical, so we want to shift the experimental DOWN
			std::vector<double> originalMzs = oldScan->MassSpectrum.XArray;
			std::vector<double> calibratedMzs(originalMzs.size());
			//calibrate the mzs. Because peaks are in mz and we are making a mz shift, we don't need to deconvolute
			for (int i = 0; i < originalMzs.size(); i++)
			{
				calibratedMzs[i] = originalMzs[i] * correctionFactor;
			}

			//update precursor values (if applicable)
			std::optional<double> selectedIonMz = std::nullopt;
			std::optional<double> isolationMz = std::nullopt;
			std::optional<double> selectedIonMonoisotopicGuessMz = std::nullopt;
			if (precursorSmoothedRelativeError)
			{
				correctionFactor = 1 - precursorSmoothedRelativeError.value();
				if (oldScan->SelectedIonMZ.HasValue)
				{
					selectedIonMz = oldScan->SelectedIonMZ * correctionFactor;
				}
				if (oldScan->IsolationMz.HasValue)
				{
					isolationMz = oldScan->IsolationMz * correctionFactor;
				}
				if (oldScan->SelectedIonMZ.HasValue)
				{
					selectedIonMonoisotopicGuessMz = oldScan->SelectedIonMonoisotopicGuessMz * correctionFactor;
				}
			}

			//create new calibrated spectrum
			MzSpectrum *calibratedSpectrum = new MzSpectrum(calibratedMzs, oldScan->MassSpectrum.YArray, false);

			delete calibratedSpectrum;
			return new MsDataScan(calibratedSpectrum, oldScan->OneBasedScanNumber, oldScan->MsnOrder, oldScan->IsCentroid, oldScan->Polarity, oldScan->RetentionTime, oldScan->ScanWindowRange, oldScan->ScanFilter, oldScan->MzAnalyzer, oldScan->TotalIonCurrent, oldScan->InjectionTime, oldScan->NoiseData, oldScan->NativeId, selectedIonMz, oldScan->SelectedIonChargeStateGuess, oldScan->SelectedIonIntensity, isolationMz, oldScan->IsolationWidth, oldScan->DissociationType, oldScan->OneBasedPrecursorScanNumber, selectedIonMonoisotopicGuessMz);
		}
	}
}
