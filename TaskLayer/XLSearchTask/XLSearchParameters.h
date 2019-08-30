#pragma once

#include <string>
#include <optional>

using namespace EngineLayer::CrosslinkSearch;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{
	class XlSearchParameters
	{
	private:
		DecoyType *privateDecoyType;
		CrosslinkerType privateCrosslinkerType = static_cast<CrosslinkerType>(0);
		int privateCrosslinkSearchTopNum = 0;
		std::string privateCrosslinkerName;
		std::optional<double> privateCrosslinkerTotalMass;
		std::optional<double> privateCrosslinkerShortMass;
		std::optional<double> privateCrosslinkerLongMass;
		std::optional<double> privateCrosslinkerLoopMass;
		std::string privateCrosslinkerResidues;
		std::string privateCrosslinkerResidues2;
		std::optional<double> privateCrosslinkerDeadEndMassH2O;
		std::optional<double> privateCrosslinkerDeadEndMassNH2;
		std::optional<double> privateCrosslinkerDeadEndMassTris;
		bool privateIsCleavable = false;
		bool privateRestrictToTopNHits = false;
		bool privateWriteOutputForPercolator = false;
		bool privateWritePepXml = false;
		bool privateXlQuench_H2O = false;
		bool privateXlQuench_Tris = false;
		bool privateXlQuench_NH2 = false;

	public:
		XlSearchParameters();

		DecoyType *getDecoyType() const;
		void setDecoyType(DecoyType *value);
		CrosslinkerType getCrosslinkerType() const;
		void setCrosslinkerType(CrosslinkerType value);
		int getCrosslinkSearchTopNum() const;
		void setCrosslinkSearchTopNum(int value);
		std::string getCrosslinkerName() const;
		void setCrosslinkerName(const std::string &value);
		std::optional<double> getCrosslinkerTotalMass() const;
		void setCrosslinkerTotalMass(std::optional<double> value);
		std::optional<double> getCrosslinkerShortMass() const;
		void setCrosslinkerShortMass(std::optional<double> value);
		std::optional<double> getCrosslinkerLongMass() const;
		void setCrosslinkerLongMass(std::optional<double> value);
		std::optional<double> getCrosslinkerLoopMass() const;
		void setCrosslinkerLoopMass(std::optional<double> value);
		std::string getCrosslinkerResidues() const;
		void setCrosslinkerResidues(const std::string &value);
		std::string getCrosslinkerResidues2() const;
		void setCrosslinkerResidues2(const std::string &value);
		std::optional<double> getCrosslinkerDeadEndMassH2O() const;
		void setCrosslinkerDeadEndMassH2O(std::optional<double> value);
		std::optional<double> getCrosslinkerDeadEndMassNH2() const;
		void setCrosslinkerDeadEndMassNH2(std::optional<double> value);
		std::optional<double> getCrosslinkerDeadEndMassTris() const;
		void setCrosslinkerDeadEndMassTris(std::optional<double> value);

		// TODO: 2+ 3+ prime fragments?
		bool getIsCleavable() const;
		void setIsCleavable(bool value);
		bool getRestrictToTopNHits() const;
		void setRestrictToTopNHits(bool value);
		bool getWriteOutputForPercolator() const;
		void setWriteOutputForPercolator(bool value);
		bool getWritePepXml() const;
		void setWritePepXml(bool value);
		bool getXlQuench_H2O() const;
		void setXlQuench_H2O(bool value);
		bool getXlQuench_Tris() const;
		void setXlQuench_Tris(bool value);
		bool getXlQuench_NH2() const;
		void setXlQuench_NH2(bool value);
		//public bool XlCharge_2_3 { get; set; }
	};
}
