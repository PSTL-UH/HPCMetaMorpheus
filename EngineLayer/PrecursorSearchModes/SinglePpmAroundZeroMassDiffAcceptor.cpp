#include "SinglePpmAroundZeroMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

	SinglePpmAroundZeroSearchMode::SinglePpmAroundZeroSearchMode(double ppmTolerance) : MassDiffAcceptor(std::to_wstring(ppmTolerance) + L"ppmAroundZero"), PpmTolerance(ppmTolerance)
	{
	}

	int SinglePpmAroundZeroSearchMode::Accepts(double scanPrecursorMass, double peptideMass)
	{
		return std::abs(((scanPrecursorMass - peptideMass) / (peptideMass)) * 1e6) < PpmTolerance ? 0 : -1;
	}

	std::vector<AllowedIntervalWithNotch*> SinglePpmAroundZeroSearchMode::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
	{
		auto diff = PpmTolerance / 1e6 * peptideMonoisotopicMass;
		DoubleRange tempVar(peptideMonoisotopicMass - diff, peptideMonoisotopicMass + diff);
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 0);
	}

	std::vector<AllowedIntervalWithNotch*> SinglePpmAroundZeroSearchMode::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
	{
		auto diff = PpmTolerance / 1e6 * peptideMonoisotopicMass;
		DoubleRange tempVar(peptideMonoisotopicMass - diff, peptideMonoisotopicMass + diff);
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 0);
	}

	std::wstring SinglePpmAroundZeroSearchMode::ToProseString()
	{
		return (std::wstring::Format(L"{0:0.0}", PpmTolerance) + L" ppm around zero");
	}

	std::wstring SinglePpmAroundZeroSearchMode::ToString()
	{
		return getFileNameAddition() + L" ppmAroundZero " + std::to_wstring(PpmTolerance);
	}
}
