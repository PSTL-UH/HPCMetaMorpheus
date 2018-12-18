#include "LabeledDataPoint.h"
#include "../PeptideSpectralMatch.h"

namespace EngineLayer
{
	namespace Calibration
	{

		LabeledDataPoint::LabeledDataPoint(double experimentalMz, int scanNumber, double logTotalIonCurrent, double logInjectionTime, double logIntensity, double theoreticalMz, PeptideSpectralMatch *identification) : ExperimentalMz(experimentalMz), ScanNumber(scanNumber), LogTotalIonCurrent(logTotalIonCurrent), LogInjectionTime(logInjectionTime), LogIntensity(logIntensity), TheoreticalMz(theoreticalMz), RelativeMzError((experimentalMz - theoreticalMz) / theoreticalMz), Identification(identification)
		{
		}
	}
}
