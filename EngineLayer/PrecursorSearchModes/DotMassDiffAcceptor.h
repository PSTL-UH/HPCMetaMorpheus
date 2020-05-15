#pragma once

#include "MassDiffAcceptor.h"
#include <string>
#include <vector>

#include "MzLibUtil.h"
#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{
	class DotMassDiffAcceptor : public MassDiffAcceptor
	{
	private:
		std::vector<double> AcceptableSortedMassShifts;
		Tolerance * const tolerance;

	public:
		virtual ~DotMassDiffAcceptor()
		{
			delete tolerance;
		}

		DotMassDiffAcceptor(const std::string &FileNameAddition, std::vector<double> &acceptableMassShifts, Tolerance *tol);

		int Accepts(double scanPrecursorMass, double peptideMass) override;

		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) override;

		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) override;

		std::string ToString();

		std::string ToProseString() override;
	};
}
