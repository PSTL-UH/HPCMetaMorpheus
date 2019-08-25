#include "SingleAbsoluteAroundZeroMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

	SingleAbsoluteAroundZeroSearchMode::SingleAbsoluteAroundZeroSearchMode(double value) : MassDiffAcceptor(std::to_string(value) + L"daltonsAroundZero"), Value(value)
	{
	}

	int SingleAbsoluteAroundZeroSearchMode::Accepts(double scanPrecursorMass, double peptideMass)
	{
		return std::abs(scanPrecursorMass - peptideMass) < Value ? 0 : -1;
	}

	std::vector<AllowedIntervalWithNotch*> SingleAbsoluteAroundZeroSearchMode::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
	{
		DoubleRange tempVar(peptideMonoisotopicMass - Value, peptideMonoisotopicMass + Value);
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 0);
	}

	std::vector<AllowedIntervalWithNotch*> SingleAbsoluteAroundZeroSearchMode::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
	{
		DoubleRange tempVar(peptideMonoisotopicMass - Value, peptideMonoisotopicMass + Value);
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 0);
	}

	std::string SingleAbsoluteAroundZeroSearchMode::ToString()
	{
		return getFileNameAddition() + " daltonsAroundZero " + std::to_string(Value);
	}

	std::string SingleAbsoluteAroundZeroSearchMode::ToProseString()
	{
		return (std::string::Format("{0:0.000}", Value) + " Da around zero");
	}
}
