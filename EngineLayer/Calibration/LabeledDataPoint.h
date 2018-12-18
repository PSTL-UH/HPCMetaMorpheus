#pragma once

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }

namespace EngineLayer
{
	namespace Calibration
	{
		class LabeledDataPoint
		{
		public:
			const double ExperimentalMz;
			const int ScanNumber;
			const double LogTotalIonCurrent;
			const double LogInjectionTime;
			const double LogIntensity;
			const double TheoreticalMz;
			const double RelativeMzError;
			PeptideSpectralMatch *const Identification;

			virtual ~LabeledDataPoint()
			{
				delete Identification;
			}

			LabeledDataPoint(double experimentalMz, int scanNumber, double logTotalIonCurrent, double logInjectionTime, double logIntensity, double theoreticalMz, PeptideSpectralMatch *identification);
		};
	}
}
