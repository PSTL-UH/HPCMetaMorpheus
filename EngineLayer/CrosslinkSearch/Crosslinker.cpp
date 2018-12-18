#include "Crosslinker.h"

namespace EngineLayer
{
	namespace CrosslinkSearch
	{

		Crosslinker::Crosslinker(const std::wstring &crosslinkerModSites, const std::wstring &crosslinkerModSites2, const std::wstring &crosslinkerName, bool cleavable, double totalMass, double cleaveMassShort, double cleaveMassLong, double loopMass, double deadendMassH2O, double deadendMassNH2, double deadendMassTris)
		{
			setCrosslinkerModSites(crosslinkerModSites);
			setCrosslinkerModSites2(crosslinkerModSites2);
			setCrosslinkerName(crosslinkerName);
			setCleavable(cleavable);
			setTotalMass(totalMass);
			setCleaveMassShort(cleaveMassShort);
			setCleaveMassLong(cleaveMassLong);
			setLoopMass(loopMass);
			setDeadendMassH2O(deadendMassH2O);
			setDeadendMassNH2(deadendMassNH2);
			setDeadendMassTris(deadendMassTris);
		}

		Crosslinker::Crosslinker()
		{
		}

		std::wstring Crosslinker::getCrosslinkerModSites() const
		{
			return privateCrosslinkerModSites;
		}

		void Crosslinker::setCrosslinkerModSites(const std::wstring &value)
		{
			privateCrosslinkerModSites = value;
		}

		std::wstring Crosslinker::getCrosslinkerModSites2() const
		{
			return privateCrosslinkerModSites2;
		}

		void Crosslinker::setCrosslinkerModSites2(const std::wstring &value)
		{
			privateCrosslinkerModSites2 = value;
		}

		std::wstring Crosslinker::getCrosslinkerName() const
		{
			return privateCrosslinkerName;
		}

		void Crosslinker::setCrosslinkerName(const std::wstring &value)
		{
			privateCrosslinkerName = value;
		}

		bool Crosslinker::getCleavable() const
		{
			return privateCleavable;
		}

		void Crosslinker::setCleavable(bool value)
		{
			privateCleavable = value;
		}

		double Crosslinker::getTotalMass() const
		{
			return privateTotalMass;
		}

		void Crosslinker::setTotalMass(double value)
		{
			privateTotalMass = value;
		}

		double Crosslinker::getCleaveMassShort() const
		{
			return privateCleaveMassShort;
		}

		void Crosslinker::setCleaveMassShort(double value)
		{
			privateCleaveMassShort = value;
		}

		double Crosslinker::getCleaveMassLong() const
		{
			return privateCleaveMassLong;
		}

		void Crosslinker::setCleaveMassLong(double value)
		{
			privateCleaveMassLong = value;
		}

		double Crosslinker::getLoopMass() const
		{
			return privateLoopMass;
		}

		void Crosslinker::setLoopMass(double value)
		{
			privateLoopMass = value;
		}

		double Crosslinker::getDeadendMassH2O() const
		{
			return privateDeadendMassH2O;
		}

		void Crosslinker::setDeadendMassH2O(double value)
		{
			privateDeadendMassH2O = value;
		}

		double Crosslinker::getDeadendMassNH2() const
		{
			return privateDeadendMassNH2;
		}

		void Crosslinker::setDeadendMassNH2(double value)
		{
			privateDeadendMassNH2 = value;
		}

		double Crosslinker::getDeadendMassTris() const
		{
			return privateDeadendMassTris;
		}

		void Crosslinker::setDeadendMassTris(double value)
		{
			privateDeadendMassTris = value;
		}

		Crosslinker *Crosslinker::SelectCrosslinker(CrosslinkerType type)
		{
			if (type == CrosslinkerType::DSSO)
			{
				setCrosslinkerName(L"DSSO");
				setCleavable(true);
				setTotalMass(158.0038);
				setCleaveMassShort(54.01056);
				setCleaveMassLong(103.9932);
				setCrosslinkerModSites(L"K");
				setCrosslinkerModSites2(L"K");
				setLoopMass(158.0038);
				setDeadendMassH2O(176.0143);
				setDeadendMassNH2(175.0303);
				setDeadendMassTris(279.0777);
			}
			else if (type == CrosslinkerType::DisulfideBond)
			{
				setCrosslinkerName(L"DisulfideBond");
				setCleavable(true);
				setTotalMass(-2.01565);
				setCleaveMassShort(-33.98772);
				setCleaveMassLong(31.97207);
				setCrosslinkerModSites(L"C");
				setCrosslinkerModSites2(L"C");
			}
			else if (type == CrosslinkerType::DSS)
			{
				setCrosslinkerName(L"DSS");
				setCleavable(false);
				setTotalMass(138.06808);
				setCrosslinkerModSites(L"K");
				setCrosslinkerModSites2(L"K");
				setLoopMass(138.06808);
				setDeadendMassH2O(156.0786);
				setDeadendMassNH2(155.0946);
				setDeadendMassTris(259.142);
			}
			else if (type == CrosslinkerType::DSBU)
			{
				setCrosslinkerName(L"DSBU");
				setCleavable(true);
				setTotalMass(196.0848);
				setCleaveMassShort(85.05276);
				setCleaveMassLong(111.0320);
				setCrosslinkerModSites(L"K");
				setCrosslinkerModSites2(L"K");
				setLoopMass(196.0848);
				setDeadendMassH2O(214.0954);
				setDeadendMassNH2(213.1113);
				setDeadendMassTris(317.1587);
			}

			return this;
		}
	}
}
