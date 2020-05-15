#pragma once

#include <string>
#include <vector>

#include "AllowedIntervalWithNotch.h"

namespace EngineLayer
{
    class MassDiffAcceptor
    {
    private:
        int privateNumNotches = 0;
        std::string privateFileNameAddition;
        
    protected:
        MassDiffAcceptor(const std::string &fileNameAddition);
        
    public:
        int getNumNotches() const;
        void setNumNotches(int value);
        std::string getFileNameAddition() const;
        
        /// <summary>
        /// If acceptable, returns 0 or greater, negative means does not accept
        /// </summary>
        /// <param name="scanPrecursorMass"></param>
        /// <param name="peptideMass"></param>
        /// <returns></returns>
        virtual int Accepts(double scanPrecursorMass, double peptideMass) = 0;
        
        virtual std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) = 0;
        virtual std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) = 0;
        
        virtual std::string ToProseString() = 0;
    };
}
