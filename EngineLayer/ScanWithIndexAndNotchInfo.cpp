#include "ScanWithIndexAndNotchInfo.h"
#include "Ms2ScanWithSpecificMass.h"

namespace EngineLayer
{

	ScanWithIndexAndNotchInfo::ScanWithIndexAndNotchInfo(Ms2ScanWithSpecificMass *theScan, int notch, int scanIndex)
	{
		TheScan = theScan;
		Notch = notch;
		ScanIndex = scanIndex;
	}
}
