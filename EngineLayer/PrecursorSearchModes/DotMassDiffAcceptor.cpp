#include "DotMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

	DotMassDiffAcceptor::DotMassDiffAcceptor(const std::string &FileNameAddition, std::vector<double> &acceptableMassShifts, Tolerance *tol) : MassDiffAcceptor(FileNameAddition), AcceptableSortedMassShifts(acceptableMassShifts.OrderBy([&] (std::any b)
	{
			return b;
	})->ToArray()), Tolerance(tol)
	{
		setNumNotches(AcceptableSortedMassShifts.size());
	}

	int DotMassDiffAcceptor::Accepts(double scanPrecursorMass, double peptideMass)
	{
		// index of the first element that is larger than or equal to value
		int index = Array::BinarySearch(AcceptableSortedMassShifts, scanPrecursorMass - peptideMass);
		if (index < 0)
		{
			index = ~index;
		}

		if (index < AcceptableSortedMassShifts.size() && Tolerance->Within(scanPrecursorMass, AcceptableSortedMassShifts[index] + peptideMass))
		{
			return index;
		}

		if (index > 0 && Tolerance->Within(scanPrecursorMass, AcceptableSortedMassShifts[index - 1] + peptideMass))
		{
			return index - 1;
		}
		return -1;
	}

	std::vector<AllowedIntervalWithNotch*> DotMassDiffAcceptor::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
	{
		for (int j = 0; j < AcceptableSortedMassShifts.size(); j++)
		{
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
			yield return new AllowedIntervalWithNotch(Tolerance->GetRange(peptideMonoisotopicMass + AcceptableSortedMassShifts[j]), j);
		}
	}

	std::vector<AllowedIntervalWithNotch*> DotMassDiffAcceptor::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
	{
		for (int j = 0; j < AcceptableSortedMassShifts.size(); j++)
		{
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
			yield return new AllowedIntervalWithNotch(Tolerance->GetRange(peptideMonoisotopicMass - AcceptableSortedMassShifts[j]), j);
		}
	}

	std::string DotMassDiffAcceptor::ToString()
	{
		return getFileNameAddition() + " dot " + Tolerance->ToString() + " " + std::string::Join(",", AcceptableSortedMassShifts);
	}

	std::string DotMassDiffAcceptor::ToProseString()
	{
		return (Tolerance->ToString() + " around " + std::string::Join(",", AcceptableSortedMassShifts) + " Da");
	}
}
