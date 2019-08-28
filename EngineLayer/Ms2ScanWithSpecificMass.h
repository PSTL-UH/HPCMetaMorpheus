#pragma once

#include "IScan.h"
#include <string>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <optional>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
//namespace EngineLayer { class CommonParameters; }
#include "CommonParameters.h"

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

namespace EngineLayer
{
    class Ms2ScanWithSpecificMass : public IScan
    {
    private:
        MsDataScan *privateTheScan;
        double privatePrecursorMonoisotopicPeakMz = 0;
        double privatePrecursorMass = 0;
        int privatePrecursorCharge = 0;
        std::string privateFullFilePath;
        std::vector<IsotopicEnvelope*> privateExperimentalFragments;
        
    public:
        Ms2ScanWithSpecificMass(MsDataScan *mzLibScan, double precursorMonoisotopicPeakMz, int precursorCharge, const std::string &fullFilePath, CommonParameters *commonParam, std::vector<IsotopicEnvelope*> &neutralExperimentalFragments);
        
        MsDataScan *getTheScan() const;
        double getPrecursorMonoisotopicPeakMz() const override;
        double getPrecursorMass() const override;
        int getPrecursorCharge() const override;
        std::string getFullFilePath() const override;
        std::vector<IsotopicEnvelope*> getExperimentalFragments() const;
        void setExperimentalFragments(const std::vector<IsotopicEnvelope*> &value);
    private:
        std::vector<double> DeconvolutedMonoisotopicMasses;
        
    public:
        int getOneBasedScanNumber() const override;
        
        std::optional<int> getOneBasedPrecursorScanNumber() const override;
        
        double getRetentionTime() const override;
        
        int getNumPeaks() const override;
        
        double getTotalIonCurrent() const override;
        
        static std::vector<IsotopicEnvelope*> GetNeutralExperimentalFragments(MsDataScan *scan, CommonParameters *commonParam);
        
        IsotopicEnvelope *GetClosestExperimentalFragmentMass(double theoreticalNeutralMass);
        
    private:
        std::optional<int> GetClosestFragmentMass(double mass);
    };
}
