#pragma once

#include "MassDiffAcceptor.h"
#include <string>
#include <vector>
#include <cmath>

#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{
    class SinglePpmAroundZeroSearchMode : public MassDiffAcceptor
    {
    private:
        const double PpmTolerance;
        
    public:
        SinglePpmAroundZeroSearchMode(double ppmTolerance);
        
        int Accepts(double scanPrecursorMass, double peptideMass) override;
        
        std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) override;
        std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) override;
        
        std::string ToProseString() override;
        
        std::string ToString();
    };
}
