#include "IntervalMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"
#include "Search.h"

using namespace MzLibUtil;

namespace EngineLayer
{

    IntervalMassDiffAcceptor::IntervalMassDiffAcceptor(const std::string &fileNameAddition,
                                                       std::vector<DoubleRange*> &doubleRanges) :
        MassDiffAcceptor(fileNameAddition)
    {
        Intervals = doubleRanges;
        std::sort (Intervals.begin(), Intervals.end(), [&] (DoubleRange *l, DoubleRange *r) {
                return l->getMean() < r->getMean();
                    });
           
        for ( auto b : Intervals ) {
            Means.push_back(b->getMean() );
        }
    }
    
    int IntervalMassDiffAcceptor::Accepts(double scanPrecursorMass, double peptideMass)
    {
        double diff = scanPrecursorMass - peptideMass;
        
        int index = BinarySearch(Means, diff);
        if (index >= 0)
        {
            return 0;
        }
        index = ~index;
        // Two options: either it's the index of the first element greater than diff, or len if diff greater than all
		if (index < (int)Means.size() && Intervals[index]->Contains(diff))
		{
                    return 0;
		}
		if (index > 0 && Intervals[index - 1]->Contains(diff))
		{
                    return 0;
		}
		return -1;
    }
    
    std::vector<AllowedIntervalWithNotch*> IntervalMassDiffAcceptor::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;
        for ( auto b: Intervals )  {
            auto tempVar = new DoubleRange (peptideMonoisotopicMass + b->getMinimum(),
                                            peptideMonoisotopicMass + b->getMaximum());
            result.push_back(new AllowedIntervalWithNotch(tempVar, 0));
        }
        return result;
    }
    
    std::vector<AllowedIntervalWithNotch*> IntervalMassDiffAcceptor::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
    {
        std::vector<AllowedIntervalWithNotch*> result;

        for ( auto b: Intervals )  {
            auto tempVar = new DoubleRange (peptideMonoisotopicMass - b->getMaximum(),
                                            peptideMonoisotopicMass - b->getMinimum());
            result.push_back(new AllowedIntervalWithNotch(tempVar, 0));
        }
        return result;
    }
    
    std::string IntervalMassDiffAcceptor::ToString()
    {
        std::vector<std::string> vec;
        for ( auto p : Intervals) {
            vec.push_back(p->ToString());
        }
        return getFileNameAddition() + " interval " + StringHelper::join(vec, ',');
    }

    std::string IntervalMassDiffAcceptor::ToProseString()
    {
        std::vector<std::string> vec;
        for ( auto p : Intervals) {
            vec.push_back(p->ToString());
        }
        return "the mass (Da) interval(s) " + StringHelper::join(vec, ',');
    }
}
