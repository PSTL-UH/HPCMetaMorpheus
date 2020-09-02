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
        
        CalibrationEngine::CalibrationEngine(MsDataFile *myMSDataFile, DataPointAquisitionResults *datapoints,
                                             CommonParameters *commonParameters,
                                             std::vector<std::string> nestedIds, int verbosityLevel) :
            MetaMorpheusEngine(commonParameters, nestedIds, verbosityLevel),
            MyMsDataFile(myMSDataFile), Datapoints(datapoints)
        {
        }
        
        MetaMorpheusEngineResults *CalibrationEngine::RunSpecific()
        {
            Status("Calibrating spectra");
            std::vector<LabeledDataPoint*> ms1Points = Datapoints->getMs1List();
            std::vector<LabeledDataPoint*> ms2Points = Datapoints->getMs2List();
            std::vector<MsDataScan*> originalScans = MyMsDataFile->GetAllScansList();
            
            std::vector<MsDataScan*> ms1Scans;
            std::vector<MsDataScan*> ms2Scans;
            //separate scans by msnOrder, because the mass accuracy varies between the two
            for (auto scan : originalScans)
            {
                if (scan->getMsnOrder() == 1)
                {
                    ms1Scans.push_back(scan);
                }
                else
                {
                    ms2Scans.push_back(scan);
                }
            }
            
            //create a way to go from the scan number to the scan order. This can be a single array, since there shouldn't be any overlap of scan numbers between ms1 and ms2 scans
#ifdef ORIG
            std::vector<int> scanNumberToScanPlacement = std::vector<int>(originalScans.Max([&] (std::any x)
            {
                x::OneBasedScanNumber;
            }) + 1);
#endif
            int num = 0;
            for ( auto x : originalScans ) {
                if ( x->getOneBasedScanNumber() > num ) {
                    num = x->getOneBasedScanNumber();
                }
            }
            std::vector<int> scanNumberToScanPlacement = std::vector<int>(num+1);
            
            for (int i = 0; i < (int)ms1Scans.size(); i++)
            {
                scanNumberToScanPlacement[ms1Scans[i]->getOneBasedScanNumber()] = i;
            }
            for (int i = 0; i < (int)ms2Scans.size(); i++)
            {
                scanNumberToScanPlacement[ms2Scans[i]->getOneBasedScanNumber()] = i;
            }
            
            //Populate the weighted average relative error for each scan, where index of returned array is the placement
            std::vector<double> ms1RelativeErrors = PopulateErrors(ms1Points, scanNumberToScanPlacement, ms1Scans.size());
            std::vector<double> ms2RelativeErrors = PopulateErrors(ms2Points, scanNumberToScanPlacement, ms2Scans.size());
            
            //generate new scans
            std::vector<MsDataScan*> calibratedScans(originalScans.size());
            
            //hard copy original scans
            for (int i = 0; i < (int)originalScans.size(); i++)
            {
                calibratedScans[i] = originalScans[i];
            }
            
            //apply a smoothing function, so that outlier scans aren't wildly shifted
            std::vector<double> ms1SmoothedErrors = SmoothErrors(ms1RelativeErrors);
            std::vector<double> ms2SmoothedErrors = SmoothErrors(ms2RelativeErrors);
            
            //calibrate the data
            int ms1Index = 0;
            int ms2Index = 0;
            
            //this is needed to update the precursor mass error of MS2 scans
            //double mostRecentMS1SmoothedError = ms1SmoothedErrors.FirstOrDefault(); 
            double mostRecentMS1SmoothedError = ms1SmoothedErrors.front(); 
            for (int scanIndex = 0; scanIndex < (int)calibratedScans.size(); scanIndex++) //foreach scan
            {
                MsDataScan *originalScan = originalScans[scanIndex]; //get original scan
                if (originalScan->getMsnOrder() == 1) //if ms1
                {
                    mostRecentMS1SmoothedError = ms1SmoothedErrors[ms1Index]; //update the mass error
                    std::optional<double> l = std::nullopt;
                    calibratedScans[scanIndex] = CalibrateScan(originalScan, mostRecentMS1SmoothedError, l);
                    ms1Index++;
                }
                else //if ms2
                {
                    std::optional<double> m =  std::make_optional(mostRecentMS1SmoothedError);
                    calibratedScans[scanIndex] = CalibrateScan(originalScan, ms2SmoothedErrors[ms2Index], m);
                    ms2Index++;
                }
            }
            
            MsDataFile tempVar(calibratedScans, MyMsDataFile->getSourceFile());
            setCalibratedDataFile(&tempVar);
            return new MetaMorpheusEngineResults(this);
        }
        
        std::vector<double> CalibrationEngine::PopulateErrors(std::vector<LabeledDataPoint*> &datapoints,
                                                              std::vector<int> &scanNumberToScanPlacement,
                                                              int arrayLength)
        {
            //return an array of weighted average relative errors for each scan
            std::vector<double> averageRelativeErrors(arrayLength);
            //get the first scan number. This should be an ordered list
            int currentScanNumber = datapoints.front()->ScanNumber;
            //create a list to store the mass error
            std::vector<std::tuple<double, double>> localRelativeErrors; 
            
            for (auto datapoint : datapoints)
            {
                if (datapoint->ScanNumber == currentScanNumber)
                {
                    localRelativeErrors.push_back(std::make_tuple(datapoint->RelativeMzError, datapoint->LogIntensity));
                }
                else
                {
                    //wrap up the old set
                    averageRelativeErrors[scanNumberToScanPlacement[currentScanNumber]] = CalculateAverageRelativeErrors(localRelativeErrors);
                    
                    //update
                    currentScanNumber = datapoint->ScanNumber;
                    localRelativeErrors.push_back(std::make_tuple(datapoint->RelativeMzError, datapoint->LogIntensity));
                }
            }
            //finish loop
            averageRelativeErrors[scanNumberToScanPlacement[currentScanNumber]] = CalculateAverageRelativeErrors(localRelativeErrors);
            
            return averageRelativeErrors;
        }
        
        double CalibrationEngine::CalculateAverageRelativeErrors(std::vector<std::tuple<double, double>> &localRelativeErrors)
        {
            // normalize each log intensity so that the minimum log intensity is 1. 
            //double logIntensityToSubtract = localRelativeErrors.Min(x => x.logIntensity) - 1; 
            //Convert from log to actual intensity to more heavily weight intensities.
#ifdef ORIG
            double weightedSumOfErrors = localRelativeErrors.Sum([&] (std::any x) {
                    return x::massError * std::pow(10, x::logIntensity);
                });
#endif
            double weightedSumOfErrors = 0.0;
            for ( auto x: localRelativeErrors ) {
                weightedSumOfErrors += std::get<0>(x)  * std::pow(10, std::get<1>(x));
            }
            
#ifdef ORIG
            double sumOfIntensities = localRelativeErrors.Sum([&] (std::any x) {
                    std::pow(10, x::logIntensity);
                });
#endif
            double sumOfIntensities = 0.0;
            for ( auto x: localRelativeErrors ) {
                sumOfIntensities +=  std::pow(10, std::get<1>(x));
            }
            
            return weightedSumOfErrors / sumOfIntensities;
        }
        
        std::vector<double> CalibrationEngine::SmoothErrors(std::vector<double> &relativeErrors)
        {
            //impute missing values
            //not all scans are guarenteed to contain data points. We can infer these data point with nearby points.
            for (int index = 0; index < (int)relativeErrors.size(); index++)
            {
                if (relativeErrors[index] == 0) //if there were no points, then the value should be perfectly zero (the double default)
                {
                    int startingBlankIndex = index;
                    //increase the index until we find the next scan containing a data point
                    while (index < (int)relativeErrors.size() && relativeErrors[index] == 0)
                    {
                        index++;
                    }
                    //can't go all the way through without any data points, the original function checks for
                    //enough data points (where enough is more than zero)
                    double nextError = index == (int)relativeErrors.size() ? relativeErrors[startingBlankIndex - 1] : relativeErrors[index]; 
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
            //this variable is the sum of all nearby errors, to be later divided by the number of
            //summed errors for a smoothed average error
            double smoothedCorrectionFactor = 0; 
            //for scan #1 (index 0)
            //no left scans, because we're at the beginning of the file. Just populate the
            //first "numberOfScansUsedForSmoothingOnEachSide" scan errors on the right side
                
            while (rightIndex < NumberOfScansUsedForSmoothingOnEachSide && rightIndex < (int)relativeErrors.size())
            {
                smoothedCorrectionFactor += relativeErrors[rightIndex];
                rightIndex++;
            }
            
            //for each scan
            for (int i = 0; i < (int)relativeErrors.size(); i++)
            {
                //for left index, remove an error if 
                if (i > NumberOfScansUsedForSmoothingOnEachSide)
                {
                    smoothedCorrectionFactor -= relativeErrors[leftIndex];
                    leftIndex++;
                }
                
                if (rightIndex < (int)relativeErrors.size()) //need to check so we don't run off the end of the file
                {
                    smoothedCorrectionFactor += relativeErrors[rightIndex];
                    rightIndex++;
                }
                
                smoothedErrors[i] = smoothedCorrectionFactor / (rightIndex - leftIndex);
            }
            return smoothedErrors;
        }
        
        MsDataScan *CalibrationEngine::CalibrateScan(MsDataScan *oldScan, double smoothedRelativeError,
                                                     std::optional<double> &precursorSmoothedRelativeError)
        {
            //create the multiplier. Positive mass errors mean that the experimental mass was
            //greater than the theoretical, so we want to shift the experimental DOWN
            double correctionFactor = 1 - smoothedRelativeError; 
            std::vector<double> originalMzs = oldScan->getMassSpectrum()->getXArray();
            std::vector<double> calibratedMzs(originalMzs.size());
            //calibrate the mzs. Because peaks are in mz and we are making a mz shift, we don't need to deconvolute
            for (int i = 0; i < (int)originalMzs.size(); i++)
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
                if (oldScan->getSelectedIonMZ().has_value())
                {
                    selectedIonMz = oldScan->getSelectedIonMZ().value() * correctionFactor;
                }
                if (oldScan->getIsolationMz().has_value())
                {
                    isolationMz = oldScan->getIsolationMz().value() * correctionFactor;
                }
                if (oldScan->getSelectedIonMZ().has_value())
                {
                    selectedIonMonoisotopicGuessMz = oldScan->getSelectedIonMonoisotopicGuessMz().value() * correctionFactor;
                }
            }
            
            //create new calibrated spectrum
            std::vector<double> yArray = oldScan->getMassSpectrum()->getYArray();
            MzSpectrum *calibratedSpectrum = new MzSpectrum(calibratedMzs,
                                                            yArray,
                                                            false);
            
            return new MsDataScan(calibratedSpectrum,
                                  oldScan->getOneBasedScanNumber(),
                                  oldScan->getMsnOrder(),
                                  oldScan->getIsCentroid(),
                                  oldScan->getPolarity(),
                                  oldScan->getRetentionTime(),
                                  oldScan->getScanWindowRange(),
                                  oldScan->getScanFilter(),
                                  oldScan->getMzAnalyzer(),
                                  oldScan->getTotalIonCurrent(),
                                  oldScan->getInjectionTime(),
                                  oldScan->getNoiseData(),
                                  oldScan->getNativeId(),
                                  selectedIonMz,
                                  oldScan->getSelectedIonChargeStateGuess(),
                                  oldScan->getSelectedIonIntensity(),
                                  isolationMz,
                                  oldScan->getIsolationWidth(),
                                  oldScan->getDissociationType(),
                                  oldScan->getOneBasedPrecursorScanNumber(),
                                  selectedIonMonoisotopicGuessMz);
            //delete calibratedSpectrum;
        }
    }
}
