#include "ExperimentalDesignForDataGrid.h"


namespace MetaMorpheusGUI
{

	ExperimentalDesignForDataGrid::ExperimentalDesignForDataGrid(const std::wstring &filename)
	{
		setFileName(Path::GetFileNameWithoutExtension(filename));
	}

	std::wstring ExperimentalDesignForDataGrid::getFileName() const
	{
		return privateFileName;
	}

	void ExperimentalDesignForDataGrid::setFileName(const std::wstring &value)
	{
		privateFileName = value;
	}

	std::wstring ExperimentalDesignForDataGrid::getCondition() const
	{
		return privateCondition;
	}

	void ExperimentalDesignForDataGrid::setCondition(const std::wstring &value)
	{
		privateCondition = value;
	}

	std::wstring ExperimentalDesignForDataGrid::getBiorep() const
	{
		return privateBiorep;
	}

	void ExperimentalDesignForDataGrid::setBiorep(const std::wstring &value)
	{
		privateBiorep = value;
	}

	std::wstring ExperimentalDesignForDataGrid::getFraction() const
	{
		return privateFraction;
	}

	void ExperimentalDesignForDataGrid::setFraction(const std::wstring &value)
	{
		privateFraction = value;
	}

	std::wstring ExperimentalDesignForDataGrid::getTechrep() const
	{
		return privateTechrep;
	}

	void ExperimentalDesignForDataGrid::setTechrep(const std::wstring &value)
	{
		privateTechrep = value;
	}

	void ExperimentalDesignForDataGrid::SetQconditionText(const std::wstring &text)
	{
		setCondition(text);
	}

	void ExperimentalDesignForDataGrid::SetQbioRepText(const std::wstring &text)
	{
		setBiorep(text);
	}

	void ExperimentalDesignForDataGrid::SetQfractionText(const std::wstring &text)
	{
		setFraction(text);
	}

	void ExperimentalDesignForDataGrid::SetQtechRepText(const std::wstring &text)
	{
		setTechrep(text);
	}
}
