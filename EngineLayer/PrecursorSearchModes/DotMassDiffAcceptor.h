#pragma once

#include "MassDiffAcceptor.h"
#include <string>
#include <vector>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class AllowedIntervalWithNotch; }

using namespace MzLibUtil;

namespace EngineLayer
{
	class DotMassDiffAcceptor : public MassDiffAcceptor
	{
	private:
		std::vector<double> const AcceptableSortedMassShifts;
		Tolerance *const Tolerance;

	public:
		virtual ~DotMassDiffAcceptor()
		{
			delete Tolerance;
		}

		DotMassDiffAcceptor(const std::string &FileNameAddition, std::vector<double> &acceptableMassShifts, Tolerance *tol);

		int Accepts(double scanPrecursorMass, double peptideMass) override;

		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) override;

		std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) override;

		std::string ToString();

		std::string ToProseString() override;
	};
}
