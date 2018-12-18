#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <cmath>
#include <optional>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { namespace Calibration { class DataPointAquisitionResults; } }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }
namespace EngineLayer { namespace Calibration { class LabeledDataPoint; } }

using namespace MassSpectrometry;

namespace EngineLayer
{
	namespace Calibration
	{
		class CalibrationEngine : public MetaMorpheusEngine
		{
		private:
			MsDataFile *privateCalibratedDataFile;

			static constexpr int NumberOfScansUsedForSmoothingOnEachSide = 100;
			MsDataFile *const MyMsDataFile;
			DataPointAquisitionResults *const Datapoints;
			public:
				virtual ~CalibrationEngine()
				{
					delete MyMsDataFile;
					delete Datapoints;
				}

				MsDataFile *getCalibratedDataFile() const;
				void setCalibratedDataFile(MsDataFile *value);

			CalibrationEngine(MsDataFile *myMSDataFile, DataPointAquisitionResults *datapoints, CommonParameters *commonParameters, std::vector<std::wstring> &nestedIds);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;

		private:
			static std::vector<double> PopulateErrors(std::vector<LabeledDataPoint*> &datapoints, std::vector<int> &scanNumberToScanPlacement, int arrayLength);

			static double CalculateAverageRelativeErrors(std::vector<(double massError, double logIntensity)*> &localRelativeErrors);

			static std::vector<double> SmoothErrors(std::vector<double> &relativeErrors);

			static MsDataScan *CalibrateScan(MsDataScan *oldScan, double smoothedRelativeError, std::optional<double> &precursorSmoothedRelativeError = std::nullopt);
		};
	}
}
