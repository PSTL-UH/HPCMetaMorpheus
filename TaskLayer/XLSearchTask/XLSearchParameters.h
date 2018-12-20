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
		std::wstring privateCrosslinkerName;
		Nullable<double> privateCrosslinkerTotalMass;
		Nullable<double> privateCrosslinkerShortMass;
		Nullable<double> privateCrosslinkerLongMass;
		Nullable<double> privateCrosslinkerLoopMass;
		std::wstring privateCrosslinkerResidues;
		std::wstring privateCrosslinkerResidues2;
		Nullable<double> privateCrosslinkerDeadEndMassH2O;
		Nullable<double> privateCrosslinkerDeadEndMassNH2;
		Nullable<double> privateCrosslinkerDeadEndMassTris;
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
		std::wstring getCrosslinkerName() const;
		void setCrosslinkerName(const std::wstring &value);
		Nullable<double> getCrosslinkerTotalMass() const;
		void setCrosslinkerTotalMass(Nullable<double> value);
		Nullable<double> getCrosslinkerShortMass() const;
		void setCrosslinkerShortMass(Nullable<double> value);
		Nullable<double> getCrosslinkerLongMass() const;
		void setCrosslinkerLongMass(Nullable<double> value);
		Nullable<double> getCrosslinkerLoopMass() const;
		void setCrosslinkerLoopMass(Nullable<double> value);
		std::wstring getCrosslinkerResidues() const;
		void setCrosslinkerResidues(const std::wstring &value);
		std::wstring getCrosslinkerResidues2() const;
		void setCrosslinkerResidues2(const std::wstring &value);
		Nullable<double> getCrosslinkerDeadEndMassH2O() const;
		void setCrosslinkerDeadEndMassH2O(Nullable<double> value);
		Nullable<double> getCrosslinkerDeadEndMassNH2() const;
		void setCrosslinkerDeadEndMassNH2(Nullable<double> value);
		Nullable<double> getCrosslinkerDeadEndMassTris() const;
		void setCrosslinkerDeadEndMassTris(Nullable<double> value);

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
