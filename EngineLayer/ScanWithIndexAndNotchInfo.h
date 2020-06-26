#pragma once

#include "Ms2ScanWithSpecificMass.h"

namespace EngineLayer
{
    class ScanWithIndexAndNotchInfo
    {
    public:
        Ms2ScanWithSpecificMass *TheScan;
        int Notch = 0;
        int ScanIndex = 0;
        
        virtual ~ScanWithIndexAndNotchInfo()
        {
            //delete TheScan;
        }
        
        ScanWithIndexAndNotchInfo(Ms2ScanWithSpecificMass *theScan, int notch, int scanIndex);
    };
}
