#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <cmath>
#include <any>
#include <mutex>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }
namespace EngineLayer { namespace Calibration { class LabeledDataPoint; } }

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics::AminoAcidPolymer;

namespace EngineLayer
{
	namespace Calibration
	{
		class DataPointAcquisitionEngine : public MetaMorpheusEngine
		{
		private:
			static constexpr double FineResolutionForIsotopeDistCalculation = 0.1;

			const std::vector<PeptideSpectralMatch*> GoodIdentifications;
			MsDataFile *const MyMsDataFile;
			Tolerance *const MzToleranceForMs1Search;
			const int MinMS1isotopicPeaksNeededForConfirmedIdentification;

		public:
			virtual ~DataPointAcquisitionEngine()
			{
				delete MyMsDataFile;
				delete MzToleranceForMs1Search;
			}

			DataPointAcquisitionEngine(std::vector<PeptideSpectralMatch*> &goodIdentifications, MsDataFile *myMsDataFile, Tolerance *mzToleranceForMs1Search, int minMS1isotopicPeaksNeededForConfirmedIdentification, CommonParameters *commonParameters, std::vector<std::string> &nestedIds);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;

//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//			private (List<LabeledDataPoint>, int, int) SearchMS1Spectra(double[] theoreticalMasses, double[] theoreticalIntensities, int ms2spectrumIndex, int direction, int peptideCharge, PeptideSpectralMatch identification)
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

		private:
			static std::vector<LabeledDataPoint*> SearchMS2Spectrum(MsDataScan *ms2DataScan, PeptideSpectralMatch *identification);
		};
	}
}
