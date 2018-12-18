#include "XLSearchParameters.h"

using namespace EngineLayer::CrosslinkSearch;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{

	XlSearchParameters::XlSearchParameters()
	{
		setDecoyType(getDecoyType()->Reverse);
		setCrosslinkerType(getCrosslinkerType()::DSSO);
		setRestrictToTopNHits(true);
		setCrosslinkSearchTopNum(300);
		setCrosslinkerName(L"");
		setIsCleavable(false);
		setCrosslinkerShortMass(std::nullopt);
		setCrosslinkerLongMass(std::nullopt);
		setCrosslinkerTotalMass(std::nullopt);
		setCrosslinkerDeadEndMassH2O(std::nullopt);
		setCrosslinkerDeadEndMassNH2(std::nullopt);
		setCrosslinkerDeadEndMassTris(std::nullopt);
		setCrosslinkerResidues(L"K");
		setCrosslinkerResidues2(L"K");
		setXlQuench_H2O(true);
		setXlQuench_NH2(false);
		setXlQuench_Tris(true);
		setWriteOutputForPercolator(false);
		setWritePepXml(true);
		//XlCharge_2_3 = true;
	}

	DecoyType *XlSearchParameters::getDecoyType() const
	{
		return privateDecoyType;
	}

	void XlSearchParameters::setDecoyType(DecoyType *value)
	{
		privateDecoyType = value;
	}

	CrosslinkerType XlSearchParameters::getCrosslinkerType() const
	{
		return privateCrosslinkerType;
	}

	void XlSearchParameters::setCrosslinkerType(CrosslinkerType value)
	{
		privateCrosslinkerType = value;
	}

	int XlSearchParameters::getCrosslinkSearchTopNum() const
	{
		return privateCrosslinkSearchTopNum;
	}

	void XlSearchParameters::setCrosslinkSearchTopNum(int value)
	{
		privateCrosslinkSearchTopNum = value;
	}

	std::wstring XlSearchParameters::getCrosslinkerName() const
	{
		return privateCrosslinkerName;
	}

	void XlSearchParameters::setCrosslinkerName(const std::wstring &value)
	{
		privateCrosslinkerName = value;
	}

	Nullable<double> XlSearchParameters::getCrosslinkerTotalMass() const
	{
		return privateCrosslinkerTotalMass;
	}

	void XlSearchParameters::setCrosslinkerTotalMass(Nullable<double> value)
	{
		privateCrosslinkerTotalMass = value;
	}

	Nullable<double> XlSearchParameters::getCrosslinkerShortMass() const
	{
		return privateCrosslinkerShortMass;
	}

	void XlSearchParameters::setCrosslinkerShortMass(Nullable<double> value)
	{
		privateCrosslinkerShortMass = value;
	}

	Nullable<double> XlSearchParameters::getCrosslinkerLongMass() const
	{
		return privateCrosslinkerLongMass;
	}

	void XlSearchParameters::setCrosslinkerLongMass(Nullable<double> value)
	{
		privateCrosslinkerLongMass = value;
	}

	Nullable<double> XlSearchParameters::getCrosslinkerLoopMass() const
	{
		return privateCrosslinkerLoopMass;
	}

	void XlSearchParameters::setCrosslinkerLoopMass(Nullable<double> value)
	{
		privateCrosslinkerLoopMass = value;
	}

	std::wstring XlSearchParameters::getCrosslinkerResidues() const
	{
		return privateCrosslinkerResidues;
	}

	void XlSearchParameters::setCrosslinkerResidues(const std::wstring &value)
	{
		privateCrosslinkerResidues = value;
	}

	std::wstring XlSearchParameters::getCrosslinkerResidues2() const
	{
		return privateCrosslinkerResidues2;
	}

	void XlSearchParameters::setCrosslinkerResidues2(const std::wstring &value)
	{
		privateCrosslinkerResidues2 = value;
	}

	Nullable<double> XlSearchParameters::getCrosslinkerDeadEndMassH2O() const
	{
		return privateCrosslinkerDeadEndMassH2O;
	}

	void XlSearchParameters::setCrosslinkerDeadEndMassH2O(Nullable<double> value)
	{
		privateCrosslinkerDeadEndMassH2O = value;
	}

	Nullable<double> XlSearchParameters::getCrosslinkerDeadEndMassNH2() const
	{
		return privateCrosslinkerDeadEndMassNH2;
	}

	void XlSearchParameters::setCrosslinkerDeadEndMassNH2(Nullable<double> value)
	{
		privateCrosslinkerDeadEndMassNH2 = value;
	}

	Nullable<double> XlSearchParameters::getCrosslinkerDeadEndMassTris() const
	{
		return privateCrosslinkerDeadEndMassTris;
	}

	void XlSearchParameters::setCrosslinkerDeadEndMassTris(Nullable<double> value)
	{
		privateCrosslinkerDeadEndMassTris = value;
	}

	bool XlSearchParameters::getIsCleavable() const
	{
		return privateIsCleavable;
	}

	void XlSearchParameters::setIsCleavable(bool value)
	{
		privateIsCleavable = value;
	}

	bool XlSearchParameters::getRestrictToTopNHits() const
	{
		return privateRestrictToTopNHits;
	}

	void XlSearchParameters::setRestrictToTopNHits(bool value)
	{
		privateRestrictToTopNHits = value;
	}

	bool XlSearchParameters::getWriteOutputForPercolator() const
	{
		return privateWriteOutputForPercolator;
	}

	void XlSearchParameters::setWriteOutputForPercolator(bool value)
	{
		privateWriteOutputForPercolator = value;
	}

	bool XlSearchParameters::getWritePepXml() const
	{
		return privateWritePepXml;
	}

	void XlSearchParameters::setWritePepXml(bool value)
	{
		privateWritePepXml = value;
	}

	bool XlSearchParameters::getXlQuench_H2O() const
	{
		return privateXlQuench_H2O;
	}

	void XlSearchParameters::setXlQuench_H2O(bool value)
	{
		privateXlQuench_H2O = value;
	}

	bool XlSearchParameters::getXlQuench_Tris() const
	{
		return privateXlQuench_Tris;
	}

	void XlSearchParameters::setXlQuench_Tris(bool value)
	{
		privateXlQuench_Tris = value;
	}

	bool XlSearchParameters::getXlQuench_NH2() const
	{
		return privateXlQuench_NH2;
	}

	void XlSearchParameters::setXlQuench_NH2(bool value)
	{
		privateXlQuench_NH2 = value;
	}
}
