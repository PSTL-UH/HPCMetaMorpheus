#include "LowTheoreticalDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

	OpenLowTheoSearchMode::OpenLowTheoSearchMode() : MassDiffAcceptor(L"OpenLow")
	{
	}

	int OpenLowTheoSearchMode::Accepts(double scanPrecursorMass, double peptideMass)
	{
		return scanPrecursorMass > peptideMass - 1 ? 0 : -1;
	}

	std::vector<AllowedIntervalWithNotch*> OpenLowTheoSearchMode::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
	{
		DoubleRange tempVar(-std::numeric_limits<double>::infinity(), peptideMonoisotopicMass + 1);
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 0);
	}

	std::vector<AllowedIntervalWithNotch*> OpenLowTheoSearchMode::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
	{
		DoubleRange tempVar(peptideMonoisotopicMass - 1, std::numeric_limits<double>::infinity());
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 0);
	}

	std::wstring OpenLowTheoSearchMode::ToProseString()
	{
		return (L"unboundedHigh");
	}

	std::wstring OpenLowTheoSearchMode::ToString()
	{
		return getFileNameAddition() + L" OpenHighSearch";
	}
}
