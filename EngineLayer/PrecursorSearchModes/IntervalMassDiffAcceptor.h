#pragma once

#include "MassDiffAcceptor.h"
#include <string>
#include <vector>

#include "AllowedIntervalWithNotch.h"
#include "MzLibUtil.h"
using namespace MzLibUtil;

namespace EngineLayer
{
    class IntervalMassDiffAcceptor : public MassDiffAcceptor
    {
    private:
        std::vector<DoubleRange*> Intervals;
        std::vector<double> Means;
        
    public:
        IntervalMassDiffAcceptor(const std::string &fileNameAddition, std::vector<DoubleRange*> &doubleRanges);
        
        int Accepts(double scanPrecursorMass, double peptideMass) override;
        
        std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) override;
        
        std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) override;
        
        std::string ToString();
        
        std::string ToProseString() override;
    };
}
