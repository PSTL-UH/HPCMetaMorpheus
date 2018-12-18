#include "OpenMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

	OpenSearchMode::OpenSearchMode() : MassDiffAcceptor(L"OpenSearch")
	{
	}

	int OpenSearchMode::Accepts(double scanPrecursorMass, double peptideMass)
	{
		return 0;
	}

	std::vector<AllowedIntervalWithNotch*> OpenSearchMode::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
	{
		DoubleRange tempVar(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 0);
	}

	std::vector<AllowedIntervalWithNotch*> OpenSearchMode::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
	{
		DoubleRange tempVar(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 0);
	}

	std::wstring OpenSearchMode::ToProseString()
	{
		return (L"unbounded");
	}

	std::wstring OpenSearchMode::ToString()
	{
		return getFileNameAddition() + L" OpenSearch";
	}
}
