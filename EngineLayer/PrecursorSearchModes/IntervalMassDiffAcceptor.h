#pragma once

#include "MassDiffAcceptor.h"
#include <string>
#include <vector>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class AllowedIntervalWithNotch; }

using namespace MzLibUtil;

namespace EngineLayer
{
	class IntervalMassDiffAcceptor : public MassDiffAcceptor
	{
	private:
		const std::vector<DoubleRange*> Intervals;
		std::vector<double> const Means;

	public:
		IntervalMassDiffAcceptor(const std::wstring &fileNameAddition, std::vector<DoubleRange*> &doubleRanges);

		int Accepts(double scanPrecursorMass, double peptideMass) override;

		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) override;

		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) override;

		std::wstring ToString() override;

		std::wstring ToProseString() override;
	};
}
