#include "SingleAbsoluteAroundZeroMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

	SingleAbsoluteAroundZeroSearchMode::SingleAbsoluteAroundZeroSearchMode(double value) : MassDiffAcceptor(std::to_wstring(value) + L"daltonsAroundZero"), Value(value)
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

	std::wstring SingleAbsoluteAroundZeroSearchMode::ToString()
	{
		return getFileNameAddition() + L" daltonsAroundZero " + std::to_wstring(Value);
	}

	std::wstring SingleAbsoluteAroundZeroSearchMode::ToProseString()
	{
		return (std::wstring::Format(L"{0:0.000}", Value) + L" Da around zero");
	}
}
