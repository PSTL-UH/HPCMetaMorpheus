#include "SinglePpmAroundZeroMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

    SinglePpmAroundZeroSearchMode::SinglePpmAroundZeroSearchMode(double ppmTolerance) : MassDiffAcceptor(std::to_string(ppmTolerance) + "ppmAroundZero"), PpmTolerance(ppmTolerance)
    {
    }
    
    int SinglePpmAroundZeroSearchMode::Accepts(double scanPrecursorMass, double peptideMass)
    {
        return std::abs(((scanPrecursorMass - peptideMass) / (peptideMass)) * 1e6) < PpmTolerance ? 0 : -1;
    }
    
    std::vector<AllowedIntervalWithNotch*> SinglePpmAroundZeroSearchMode::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;
        auto diff = PpmTolerance / 1e6 * peptideMonoisotopicMass;
        auto  tempVar = new DoubleRange(peptideMonoisotopicMass - diff, peptideMonoisotopicMass + diff);
        result.push_back( new AllowedIntervalWithNotch(tempVar, 0));
        return result;
    }
    
    std::vector<AllowedIntervalWithNotch*> SinglePpmAroundZeroSearchMode::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;
        auto diff = PpmTolerance / 1e6 * peptideMonoisotopicMass;
        auto  tempVar = new DoubleRange(peptideMonoisotopicMass - diff, peptideMonoisotopicMass + diff);
        result.push_back(new AllowedIntervalWithNotch(tempVar, 0));
        return result;
    }
    
    std::string SinglePpmAroundZeroSearchMode::ToProseString()
    {
        return (std::to_string(PpmTolerance) + " ppm around zero");
    }
    
    std::string SinglePpmAroundZeroSearchMode::ToString()
    {
        return getFileNameAddition() + " ppmAroundZero " + std::to_string(PpmTolerance);
    }
}
