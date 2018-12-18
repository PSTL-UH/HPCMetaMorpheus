#pragma once

#include "MassDiffAcceptor.h"
#include <string>
#include <vector>
#include <cmath>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class AllowedIntervalWithNotch; }

using namespace MzLibUtil;

namespace EngineLayer
{
	class SinglePpmAroundZeroSearchMode : public MassDiffAcceptor
	{
	private:
		const double PpmTolerance;

	public:
		SinglePpmAroundZeroSearchMode(double ppmTolerance);

		int Accepts(double scanPrecursorMass, double peptideMass) override;

		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) override;
		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) override;

		std::wstring ToProseString() override;

		std::wstring ToString() override;
	};
}
