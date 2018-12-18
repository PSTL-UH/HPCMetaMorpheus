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
	class SingleAbsoluteAroundZeroSearchMode : public MassDiffAcceptor
	{
	private:
		const double Value;

	public:
		SingleAbsoluteAroundZeroSearchMode(double value);

		int Accepts(double scanPrecursorMass, double peptideMass) override;

		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) override;
		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) override;

		std::wstring ToString() override;

		std::wstring ToProseString() override;
	};
}
