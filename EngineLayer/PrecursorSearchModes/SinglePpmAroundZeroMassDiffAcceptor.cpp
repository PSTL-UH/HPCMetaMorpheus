#include "SinglePpmAroundZeroMassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

	SinglePpmAroundZeroSearchMode::SinglePpmAroundZeroSearchMode(double ppmTolerance) : MassDiffAcceptor(std::to_string(ppmTolerance) + L"ppmAroundZero"), PpmTolerance(ppmTolerance)
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

	std::string SinglePpmAroundZeroSearchMode::ToProseString()
	{
		return (std::string::Format("{0:0.0}", PpmTolerance) + " ppm around zero");
	}

	std::string SinglePpmAroundZeroSearchMode::ToString()
	{
		return getFileNameAddition() + " ppmAroundZero " + std::to_string(PpmTolerance);
	}
}
