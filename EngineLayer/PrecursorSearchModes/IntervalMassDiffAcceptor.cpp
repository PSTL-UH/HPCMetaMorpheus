#include "IntervalMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

	IntervalMassDiffAcceptor::IntervalMassDiffAcceptor(const std::string &fileNameAddition, std::vector<DoubleRange*> &doubleRanges) : MassDiffAcceptor(fileNameAddition), Intervals(doubleRanges.OrderBy([&] (std::any b)
	{
			b::Mean;
	}).ToList()), Means(Intervals.Select([&] (std::any b)
	{
			b::Mean;
		})->ToArray())
		{
	}

	int IntervalMassDiffAcceptor::Accepts(double scanPrecursorMass, double peptideMass)
	{
		double diff = scanPrecursorMass - peptideMass;

		int index = Array::BinarySearch(Means, diff);
		if (index >= 0)
		{
			return 0;
		}
		index = ~index;
		// Two options: either it's the index of the first element greater than diff, or len if diff greater than all
		if (index < Means.size() && Intervals[index]->Contains(diff))
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
		return Intervals.Select([&] (std::any b)
		{
			DoubleRange tempVar(peptideMonoisotopicMass + b::Minimum, peptideMonoisotopicMass + b::Maximum);
			new AllowedIntervalWithNotch(&tempVar, 0);
		});
	}

	std::vector<AllowedIntervalWithNotch*> IntervalMassDiffAcceptor::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
	{
		return Intervals.Select([&] (std::any b)
		{
			DoubleRange tempVar(peptideMonoisotopicMass - b::Maximum, peptideMonoisotopicMass - b::Minimum);
			new AllowedIntervalWithNotch(&tempVar, 0);
		});
	}

	std::string IntervalMassDiffAcceptor::ToString()
	{
		return getFileNameAddition() + " interval " + std::string::Join(",", Intervals);
	}

	std::string IntervalMassDiffAcceptor::ToProseString()
	{
		return ("the mass (Da) interval(s) " + std::string::Join(", ", Intervals));
	}
}
