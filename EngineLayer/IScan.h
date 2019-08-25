#pragma once

#include <string>
#include <optional>

namespace EngineLayer
{
	class IScan
	{
	public:
		virtual std::string getFullFilePath() const = 0;
		virtual int getOneBasedScanNumber() const = 0;
		virtual std::optional<int> getOneBasedPrecursorScanNumber() const = 0;
		virtual double getRetentionTime() const = 0;
		virtual int getNumPeaks() const = 0;
		virtual double getTotalIonCurrent() const = 0;
		virtual int getPrecursorCharge() const = 0;
		virtual double getPrecursorMonoisotopicPeakMz() const = 0;
		virtual double getPrecursorMass() const = 0;
	};
}
