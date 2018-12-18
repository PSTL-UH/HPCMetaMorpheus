#pragma once

#include <string>

namespace EngineLayer
{
	namespace CrosslinkSearch
	{
		enum class CrosslinkerType
		{
			DSSO,
			DSS,
			DisulfideBond,
			DSBU,
			UserDefined
		};

		class Crosslinker
		{
		private:
			std::wstring privateCrosslinkerModSites;
			std::wstring privateCrosslinkerModSites2;
			std::wstring privateCrosslinkerName;
			bool privateCleavable = false;
			double privateTotalMass = 0;
			double privateCleaveMassShort = 0;
			double privateCleaveMassLong = 0;
			double privateLoopMass = 0;
			double privateDeadendMassH2O = 0;
			double privateDeadendMassNH2 = 0;
			double privateDeadendMassTris = 0;

		public:
			Crosslinker(const std::wstring &crosslinkerModSites, const std::wstring &crosslinkerModSites2, const std::wstring &crosslinkerName, bool cleavable, double totalMass, double cleaveMassShort, double cleaveMassLong, double loopMass, double deadendMassH2O, double deadendMassNH2, double deadendMassTris);

			Crosslinker();

				std::wstring getCrosslinkerModSites() const;
				void setCrosslinkerModSites(const std::wstring &value);
				std::wstring getCrosslinkerModSites2() const;
				void setCrosslinkerModSites2(const std::wstring &value);
				std::wstring getCrosslinkerName() const;
				void setCrosslinkerName(const std::wstring &value);
				bool getCleavable() const;
				void setCleavable(bool value);
				double getTotalMass() const;
				void setTotalMass(double value);
				double getCleaveMassShort() const;
				void setCleaveMassShort(double value);
				double getCleaveMassLong() const;
				void setCleaveMassLong(double value);
				double getLoopMass() const;
				void setLoopMass(double value);
				double getDeadendMassH2O() const;
				void setDeadendMassH2O(double value);
				double getDeadendMassNH2() const;
				void setDeadendMassNH2(double value);
				double getDeadendMassTris() const;
				void setDeadendMassTris(double value);

			Crosslinker *SelectCrosslinker(CrosslinkerType type);
		};
	}
}
