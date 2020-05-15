#include "DotMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"
#include "stringhelper.h"
#include "Search.h"

using namespace MzLibUtil;

namespace EngineLayer
{

    DotMassDiffAcceptor::DotMassDiffAcceptor(const std::string &FileNameAddition, std::vector<double> &acceptableMassShifts, Tolerance *tol) : MassDiffAcceptor(FileNameAddition), tolerance(tol)
    {
        AcceptableSortedMassShifts = acceptableMassShifts;
        std::sort (AcceptableSortedMassShifts.begin(), AcceptableSortedMassShifts.end() );
        setNumNotches(AcceptableSortedMassShifts.size());
    }

    int DotMassDiffAcceptor::Accepts(double scanPrecursorMass, double peptideMass)
    {
        // index of the first element that is larger than or equal to value
        int index = BinarySearch(AcceptableSortedMassShifts, scanPrecursorMass - peptideMass);
        if (index < 0)
        {
            index = ~index;
        }
        
        if (index < (int)AcceptableSortedMassShifts.size() &&
            tolerance->Within(scanPrecursorMass, AcceptableSortedMassShifts[index] + peptideMass))
        {
            return index;
        }

        if (index > 0 && tolerance->Within(scanPrecursorMass, AcceptableSortedMassShifts[index - 1] + peptideMass))
        {
            return index - 1;
        }
        return -1;
    }
    
    std::vector<AllowedIntervalWithNotch*> DotMassDiffAcceptor::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;
        for (int j = 0; j < (int)AcceptableSortedMassShifts.size(); j++)
        {
            //yield return new AllowedIntervalWithNotch(Tolerance->GetRange(peptideMonoisotopicMass +
            //                                          AcceptableSortedMassShifts[j]), j);
            auto elem = new AllowedIntervalWithNotch(tolerance->GetRange(peptideMonoisotopicMass +
                                                                         AcceptableSortedMassShifts[j]), j);
            result.push_back(elem);
        }
        return result;
    }

    std::vector<AllowedIntervalWithNotch*> DotMassDiffAcceptor::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
    {        
        std::vector<AllowedIntervalWithNotch*> result;
        for (int j = 0; j < (int)AcceptableSortedMassShifts.size(); j++)
        {                
            //yield return new AllowedIntervalWithNotch(Tolerance->GetRange(peptideMonoisotopicMass -
            //                                          AcceptableSortedMassShifts[j]), j);
            auto elem = new AllowedIntervalWithNotch(tolerance->GetRange(peptideMonoisotopicMass -
                                                                         AcceptableSortedMassShifts[j]), j);
            result.push_back(elem);
        }
        return result;
    }
    
    std::string DotMassDiffAcceptor::ToString()
    {
        std::vector<std::string> vec;
        for ( auto p:  AcceptableSortedMassShifts ) {
            vec.push_back(std::to_string(p));
        }
        return getFileNameAddition() + " dot " + std::to_string(tolerance->getValue()) + " " +
            StringHelper::join( vec, ',');
    }
    
    std::string DotMassDiffAcceptor::ToProseString()
    {
        std::vector<std::string> vec;
        for ( auto p:  AcceptableSortedMassShifts ) {
            vec.push_back(std::to_string(p));
        }
        return (std::to_string(tolerance->getValue()) + " around " +
                StringHelper::join( vec, ',') + " Da");
    }
}
