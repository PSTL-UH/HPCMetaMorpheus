#include "DataPointAcquisitionEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "LabeledDataPoint.h"
#include "DataPointAquisitionResults.h"
#include "../GlobalVariables.h"

#include "Chemistry/ClassExtensions.h"
#include "Sort.h"

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics::AminoAcidPolymer;

namespace EngineLayer
{
    namespace Calibration
    {
        
        DataPointAcquisitionEngine::DataPointAcquisitionEngine(std::vector<PeptideSpectralMatch*> &goodIdentifications,
                                                               MsDataFile *myMsDataFile, Tolerance *mzToleranceForMs1Search,
                                                               int minMS1isotopicPeaksNeededForConfirmedIdentification,
                                                               CommonParameters *commonParameters,
                                         std::vector<std::string> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds),
                                                           GoodIdentifications(goodIdentifications),
                                                           MyMsDataFile(myMsDataFile),
                                                           MzToleranceForMs1Search(mzToleranceForMs1Search),
                                                           MinMS1isotopicPeaksNeededForConfirmedIdentification(minMS1isotopicPeaksNeededForConfirmedIdentification)
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
            
#ifdef ORIG
            std::any lockObj = std::any();
            std::any lockObj2 = std::any();
            ParallelOptions *tempVar = new ParallelOptions();
            tempVar->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
            Parallel::ForEach(Partitioner::Create(0, numIdentifications), tempVar, [&] (fff, loopState)   );
#endif
            // for (int matchIndex = fff::Item1; matchIndex < fff::Item2; matchIndex++)
            for (int matchIndex = 0; matchIndex < numIdentifications; matchIndex++)
            {
                // Stop loop if canceled
                if (GlobalVariables::getStopLoops())
                {
                    //loopState::Stop();
                    //return;
                    break;
                }
                
                PeptideSpectralMatch *identification = GoodIdentifications[matchIndex];
                
                // Each identification has an MS2 spectrum attached to it.
                int ms2scanNumber = identification->getScanNumber();
                int peptideCharge = identification->getScanPrecursorCharge();
                if (identification->getFullSequence() == "")
                {
                    continue;
                }
                
                auto representativeSinglePeptide = std::get<1>(identification->getBestMatchingPeptides().front());
                
                // Get the peptide, don't forget to add the modifications!!!!
                auto SequenceWithChemicalFormulas = representativeSinglePeptide->getSequenceWithChemicalFormulas();
#ifdef ORIG
                if (SequenceWithChemicalFormulas == nullptr ||
                    representativeSinglePeptide->AllModsOneIsNterminus.Any([&] (std::any b)  {
                            return b->Value->NeutralLosses != nullptr;
                        }));
#endif
                bool l = false;
                for ( auto b: representativeSinglePeptide->getAllModsOneIsNterminus() ) {
                    if ( !b.second->getNeutralLosses().empty()  ) {
                        l = true;
                        break;
                    }
                }
                if ( SequenceWithChemicalFormulas.empty() || l)
                {
                    continue;
                }
                
                Peptide *coolPeptide = new Peptide(SequenceWithChemicalFormulas);
                
                auto ms2tuple = SearchMS2Spectrum(MyMsDataFile->GetOneBasedScan(ms2scanNumber), identification);
                
                {
                    //std::lock_guard<std::mutex> lock(lockObj2);
                    //Ms2List.AddRange(ms2tuple);
                    for ( auto m: ms2tuple ) {
                        Ms2List.push_back(m);
                    }
                }
                
                // Calculate isotopic distribution of the full peptide
                auto dist = IsotopicDistribution::GetDistribution(coolPeptide->GetChemicalFormula(),
                                                                  FineResolutionForIsotopeDistCalculation, 0.001);
                
                std::vector<double> theoreticalMasses = dist->getMasses(); //->ToArray();
                std::vector<double> theoreticalIntensities = dist->getIntensities(); // ->ToArray();
                
#ifdef ORIG
                Array::Sort(theoreticalIntensities, theoreticalMasses, Comparer<double>::Create([&] (x, y ) {
                            y->CompareTo(x);
                        }));
#endif
                Sort::SortPairs(theoreticalIntensities, theoreticalMasses, theoreticalMasses.size());
                
                std::tuple<std::vector<LabeledDataPoint>, int, int> ms1tupleBack = SearchMS1Spectra(theoreticalMasses,
                                                                                                    theoreticalIntensities,
                                                                                                    ms2scanNumber, -1,
                                                                                                    peptideCharge,
                                                                                                    *identification);
                
                std::tuple<std::vector<LabeledDataPoint>, int, int> ms1tupleForward = SearchMS1Spectra(theoreticalMasses,
                                                                                                       theoreticalIntensities,
                                                                                                       ms2scanNumber, 1,
                                                                                                       peptideCharge,
                                                                                                       *identification);
                
                {
                    //std::lock_guard<std::mutex> lock(lockObj);
                    //Ms1List.AddRange(ms1tupleBack->Item1);
                    for ( auto m : std::get<0>(ms1tupleBack) ) {
                        Ms1List.push_back(&m);
                    }
                    
                    //numMs1MassChargeCombinationsConsidered += ms1tupleBack->Item2;
                    numMs1MassChargeCombinationsConsidered += std::get<1>(ms1tupleBack);
                    numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += std::get<2>(ms1tupleBack);
                    
                    //Ms1List.AddRange(ms1tupleForward->Item1);
                    for ( auto m : std::get<0>(ms1tupleForward) ) {
                        Ms1List.push_back(&m);
                    }
                    
                    numMs1MassChargeCombinationsConsidered += std::get<1>(ms1tupleForward);
                    numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += std::get<2>(ms1tupleForward);
                }
            }
            
            
            // datapoints are ordered because they were acquired in a parallized search and we want repeatable results
#ifdef ORIG
            return new DataPointAquisitionResults(this, GoodIdentifications, Ms1List.OrderBy([&] (std::any p) {
                        p::ScanNumber;
                    }).ThenBy([&] (std::any p) {
                            p::ExperimentalMz;
                        }).ToList(), Ms2List.OrderBy([&] (std::any p)	{
                                p::ScanNumber;
                            }).ThenBy([&] (std::any p)  {
                                    p::ExperimentalMz;
                                }).ToList(), numMs1MassChargeCombinationsConsidered,
                numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
                numMs2MassChargeCombinationsConsidered,
                numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
#endif
            // OrderBy and ThenBy
            std::stable_sort (Ms1List.begin(), Ms1List.end(), [&] (LabeledDataPoint* p1, LabeledDataPoint* p2) {
                    auto p1p = p1->ScanNumber;
                    auto p2p = p2->ScanNumber;
                    if ( p1p < p2p ) return true;
                    if ( p1p > p2p ) return true;

                    return p1->ExperimentalMz<  p2->ExperimentalMz;
                });
            std::stable_sort (Ms2List.begin(), Ms2List.end(), [&] (LabeledDataPoint* p1, LabeledDataPoint* p2) {
                    auto p1p = p1->ScanNumber;
                    auto p2p = p2->ScanNumber;
                    if ( p1p < p2p ) return true;
                    if ( p1p > p2p ) return true;

                    return p1->ExperimentalMz<  p2->ExperimentalMz;
                });
            return new DataPointAquisitionResults(this, GoodIdentifications, Ms1List, Ms2List,
                                                  numMs1MassChargeCombinationsConsidered,
                                                  numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
                                                  numMs2MassChargeCombinationsConsidered,
                                                  numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            
        }
    
    
        std::vector<LabeledDataPoint*> DataPointAcquisitionEngine::SearchMS2Spectrum(MsDataScan *ms2DataScan,
                                                                                     PeptideSpectralMatch *identification)
        {
            std::vector<LabeledDataPoint*> result;
            
            if (ms2DataScan->getMassSpectrum()->getSize() == 0)
            {
                return result;
            }
            
            for (auto matchedIon : identification->getMatchedFragmentIons())
            {
                double exptPeakMz = matchedIon->Mz;
                double exptPeakIntensity = matchedIon->Intensity;
                double injTime = ms2DataScan->getInjectionTime().has_value() ? ms2DataScan->getInjectionTime().value() : NAN;
                
                LabeledDataPoint tempVar(exptPeakMz, ms2DataScan->getOneBasedScanNumber(),
                                         std::log(ms2DataScan->getTotalIonCurrent()), std::log(injTime),
                                         std::log(exptPeakIntensity),
                                         Chemistry::ClassExtensions::ToMz(matchedIon->NeutralTheoreticalProduct->NeutralMass, matchedIon->Charge),
                                         identification);
                result.push_back(&tempVar);
            }
            return result;
        }
        
        
        //private (List<LabeledDataPoint>, int, int) SearchMS1Spectra(double[] theoreticalMasses, double[] theoreticalIntensities,
        //         int ms2spectrumIndex, int direction, int peptideCharge, PeptideSpectralMatch identification)
        std::tuple<std::vector<LabeledDataPoint>, int, int> DataPointAcquisitionEngine::SearchMS1Spectra(std::vector<double> theoreticalMasses,
                                                                                                         std::vector<double> theoreticalIntensities,
                                                                                                         int ms2spectrumIndex, int direction,
                                                                                                         int peptideCharge,
                                                                                                         PeptideSpectralMatch identification)
        {
            std::vector<LabeledDataPoint> result;
            int numMs1MassChargeCombinationsConsidered = 0;
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
            
            int theIndex;
            theIndex = direction == 1 ? ms2spectrumIndex : ms2spectrumIndex - 1;
            bool addedAscan = true;
            
            int highestKnownChargeForThisPeptide = peptideCharge;
            while (theIndex >= 1 && theIndex <= MyMsDataFile->getNumSpectra() && addedAscan)
            {
                int countForThisScan = 0;
                if (MyMsDataFile->GetOneBasedScan(theIndex)->getMsnOrder() > 1)
                {
                    theIndex += direction;
                    continue;
                }
                addedAscan = false;
                auto fullMS1scan = MyMsDataFile->GetOneBasedScan(theIndex);
                auto scanWindowRange = fullMS1scan->getScanWindowRange();
                auto fullMS1spectrum = fullMS1scan->getMassSpectrum();
                if (fullMS1spectrum->getSize() == 0)
                    break;
                
                bool startingToAddCharges = false;
                int chargeToLookAt = 1;
                do
                {
                    if (Chemistry::ClassExtensions::ToMz(theoreticalMasses[0], chargeToLookAt) > scanWindowRange->getMaximum())
                    {
                        chargeToLookAt++;
                        continue;
                    }
                    if (Chemistry::ClassExtensions::ToMz(theoreticalMasses[0], chargeToLookAt) < scanWindowRange->getMinimum())
                        break;
                    auto trainingPointsToAverage = new std::vector<LabeledDataPoint>();
                    for (double a : theoreticalMasses)
                    {
                        double theMZ = Chemistry::ClassExtensions::ToMz(a,chargeToLookAt);
                        
                        auto npwr = fullMS1spectrum->NumPeaksWithinRange(MzToleranceForMs1Search->GetMinimumValue(theMZ),
                                                                        MzToleranceForMs1Search->GetMaximumValue(theMZ));
                        if (npwr == 0)
                        {
                            break;
                        }
                        numMs1MassChargeCombinationsConsidered++;
                        if (npwr > 1)
                        {
                            numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks++;
                            continue;
                        }
                        
                        auto closestPeakIndex = fullMS1spectrum->GetClosestPeakIndex(theMZ);
                        auto closestPeakMZ = fullMS1spectrum->getXArray()[closestPeakIndex.value()];
                        
                        highestKnownChargeForThisPeptide = std::max(highestKnownChargeForThisPeptide, chargeToLookAt);
                        LabeledDataPoint *ldp = new LabeledDataPoint(closestPeakMZ, -1, NAN, NAN,
                                                                     std::log(fullMS1spectrum->getYArray()[closestPeakIndex.value()]),
                                                                     theMZ, nullptr);
                        trainingPointsToAverage->push_back(*ldp);
                    }
                    // If started adding and suddnely stopped, go to next one, no need to look at higher charges
                    if (trainingPointsToAverage->size() == 0 && startingToAddCharges)
                    {
                        break;
                    }
                    if ((trainingPointsToAverage->size() == 0 ||
                         (trainingPointsToAverage->size() == 1 && theoreticalIntensities[0] < 0.65)) &&
                        (peptideCharge <= chargeToLookAt))
                    {
                        break;
                    }
                    if ((trainingPointsToAverage->size() == 1 && theoreticalIntensities[0] < 0.65) ||
                        (int)trainingPointsToAverage->size() < std::min(MinMS1isotopicPeaksNeededForConfirmedIdentification, (int)theoreticalIntensities.size()))
                    {
                    }
                    else
                    {
                        addedAscan = true;
                        startingToAddCharges = true;
                        countForThisScan++;
#ifdef ORIG
                        result.Add(new LabeledDataPoint(trainingPointsToAverage.Select(b => b.ExperimentalMz).Average(),
                                                        fullMS1scan.OneBasedScanNumber,
                                                        Math.Log(fullMS1scan.TotalIonCurrent),
                                                        fullMS1scan.InjectionTime.HasValue ? Math.Log(fullMS1scan.InjectionTime.Value) : double.NaN,
                                                        trainingPointsToAverage.Select(b => b.LogIntensity).Average(),
                                                        trainingPointsToAverage.Select(b => b.TheoreticalMz).Average(),
                                                        identification));
#endif
                        
                        double average1 = 0.0, average2 = 0.0, average3 = 0.0;
                        for ( auto b = trainingPointsToAverage->begin(); b != trainingPointsToAverage->end(); b++  ) {
                            average1 += b->ExperimentalMz;
                            average2 += b->LogIntensity;
                            average3 += b->TheoreticalMz;
                        }
                        average1 = average1 / trainingPointsToAverage->size();
                        average2 = average2 / trainingPointsToAverage->size();
                        average3 = average3 / trainingPointsToAverage->size();
                        
                        LabeledDataPoint *ldp = new LabeledDataPoint(average1,
                                                          fullMS1scan->getOneBasedScanNumber(),
                                                          std::log(fullMS1scan->getTotalIonCurrent()),
                                                          fullMS1scan->getInjectionTime().has_value() ? std::log(fullMS1scan->getInjectionTime().value()) : NAN,
                                                          average2, average3, &identification);
                         result.push_back (*ldp);
                    }
                    chargeToLookAt++;
                } while (chargeToLookAt <= highestKnownChargeForThisPeptide + 1);
                theIndex += direction;
            }
            return (std::make_tuple(result, numMs1MassChargeCombinationsConsidered, numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks));
        }    
    }
}
