#include "LowTheoreticalDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

    OpenLowTheoSearchMode::OpenLowTheoSearchMode() : MassDiffAcceptor("OpenLow")
    {
    }
    
    int OpenLowTheoSearchMode::Accepts(double scanPrecursorMass, double peptideMass)
    {
        return scanPrecursorMass > peptideMass - 1 ? 0 : -1;
    }
    
    std::vector<AllowedIntervalWithNotch*> OpenLowTheoSearchMode::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;
        auto  tempVar = new DoubleRange(-std::numeric_limits<double>::infinity(), peptideMonoisotopicMass + 1);
        result.push_back(new AllowedIntervalWithNotch(tempVar, 0));
        return result;
    }

    std::vector<AllowedIntervalWithNotch*> OpenLowTheoSearchMode::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;
        auto tempVar = new DoubleRange (peptideMonoisotopicMass - 1, std::numeric_limits<double>::infinity());
        result.push_back( new AllowedIntervalWithNotch(tempVar, 0));
        return result;
    }

    std::string OpenLowTheoSearchMode::ToProseString()
    {
        return ("unboundedHigh");
    }
    
    std::string OpenLowTheoSearchMode::ToString()
    {
        return getFileNameAddition() + " OpenHighSearch";
    }
}
