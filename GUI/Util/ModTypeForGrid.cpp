#include "ModTypeForGrid.h"

namespace MetaMorpheusGUI
{

	ModTypeForGrid::ModTypeForGrid(const std::wstring &modName)
	{
		setModName(modName);
		setItem2(true);
	}

	std::wstring ModTypeForGrid::getModName() const
	{
		return privateModName;
	}

	void ModTypeForGrid::setModName(const std::wstring &value)
	{
		privateModName = value;
	}

	bool ModTypeForGrid::getItem2() const
	{
		return privateItem2;
	}

	void ModTypeForGrid::setItem2(bool value)
	{
		privateItem2 = value;
	}

	bool ModTypeForGrid::getItem3() const
	{
		return privateItem3;
	}

	void ModTypeForGrid::setItem3(bool value)
	{
		privateItem3 = value;
	}

	bool ModTypeForGrid::getItem4() const
	{
		return privateItem4;
	}

	void ModTypeForGrid::setItem4(bool value)
	{
		privateItem4 = value;
	}

	bool ModTypeForGrid::getItem5() const
	{
		return privateItem5;
	}

	void ModTypeForGrid::setItem5(bool value)
	{
		privateItem5 = value;
	}
}
