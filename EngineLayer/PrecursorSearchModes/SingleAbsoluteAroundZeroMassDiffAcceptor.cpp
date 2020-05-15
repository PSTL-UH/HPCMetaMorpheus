#include "SingleAbsoluteAroundZeroMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

    SingleAbsoluteAroundZeroSearchMode::SingleAbsoluteAroundZeroSearchMode(double value) : MassDiffAcceptor(std::to_string(value) + "daltonsAroundZero"), Value(value)
    {
    }
    
    int SingleAbsoluteAroundZeroSearchMode::Accepts(double scanPrecursorMass, double peptideMass)
	{
            return std::abs(scanPrecursorMass - peptideMass) < Value ? 0 : -1;
	}
    
    std::vector<AllowedIntervalWithNotch*> SingleAbsoluteAroundZeroSearchMode::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;            
        auto tempVar = new DoubleRange (peptideMonoisotopicMass - Value, peptideMonoisotopicMass + Value);
        result.push_back(new AllowedIntervalWithNotch(tempVar, 0));
        return result;
    }
    
    std::vector<AllowedIntervalWithNotch*> SingleAbsoluteAroundZeroSearchMode::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;
        auto  tempVar = new DoubleRange (peptideMonoisotopicMass - Value, peptideMonoisotopicMass + Value);
        result.push_back(new AllowedIntervalWithNotch(tempVar, 0));
        return result;
    }
    
    std::string SingleAbsoluteAroundZeroSearchMode::ToString()
    {
        return getFileNameAddition() + " daltonsAroundZero " + std::to_string(Value);
    }
    
    std::string SingleAbsoluteAroundZeroSearchMode::ToProseString()
    {
        return (std::to_string(Value) + " Da around zero");
    }
}
