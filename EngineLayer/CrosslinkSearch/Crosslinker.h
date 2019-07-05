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
			std::string privateCrosslinkerModSites;
			std::string privateCrosslinkerModSites2;
			std::string privateCrosslinkerName;
			bool privateCleavable = false;
			double privateTotalMass = 0;
			double privateCleaveMassShort = 0;
			double privateCleaveMassLong = 0;
			double privateLoopMass = 0;
			double privateDeadendMassH2O = 0;
			double privateDeadendMassNH2 = 0;
			double privateDeadendMassTris = 0;

		public:
			Crosslinker(const std::string &crosslinkerModSites, const std::string &crosslinkerModSites2, const std::string &crosslinkerName, bool cleavable, double totalMass, double cleaveMassShort, double cleaveMassLong, double loopMass, double deadendMassH2O, double deadendMassNH2, double deadendMassTris);

			Crosslinker();

				std::string getCrosslinkerModSites() const;
				void setCrosslinkerModSites(const std::string &value);
				std::string getCrosslinkerModSites2() const;
				void setCrosslinkerModSites2(const std::string &value);
				std::string getCrosslinkerName() const;
				void setCrosslinkerName(const std::string &value);
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
