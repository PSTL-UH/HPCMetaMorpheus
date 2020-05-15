#include "OpenMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

    OpenSearchMode::OpenSearchMode() : MassDiffAcceptor("OpenSearch")
    {
    }
    
    int OpenSearchMode::Accepts(double scanPrecursorMass, double peptideMass)
    {
        return 0;
    }
    
    std::vector<AllowedIntervalWithNotch*> OpenSearchMode::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;
        auto  tempVar = new DoubleRange(-std::numeric_limits<double>::infinity(),
                                        std::numeric_limits<double>::infinity());
        result.push_back( new AllowedIntervalWithNotch(tempVar, 0));
        return result;
    }
    
    std::vector<AllowedIntervalWithNotch*> OpenSearchMode::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;
        auto tempVar = new DoubleRange (-std::numeric_limits<double>::infinity(),
                                        std::numeric_limits<double>::infinity());
        result.push_back( new AllowedIntervalWithNotch(tempVar, 0));
        return result;
    }
    
    std::string OpenSearchMode::ToProseString()
    {
        return ("unbounded");
    }
    
    std::string OpenSearchMode::ToString()
    {
        return getFileNameAddition() + " OpenSearch";
    }
}
