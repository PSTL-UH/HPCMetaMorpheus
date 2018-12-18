#include "ModListForCalibrationTask.h"

namespace MetaMorpheusGUI
{

	ModListForCalibrationTask::ModListForCalibrationTask(const std::wstring &filePath)
	{
		this->FileName = filePath;
	}

	bool ModListForCalibrationTask::getFixed() const
	{
		return privateFixed;
	}

	void ModListForCalibrationTask::setFixed(bool value)
	{
		privateFixed = value;
	}

	bool ModListForCalibrationTask::getVariable() const
	{
		return privateVariable;
	}

	void ModListForCalibrationTask::setVariable(bool value)
	{
		privateVariable = value;
	}

	bool ModListForCalibrationTask::getLocalize() const
	{
		return privateLocalize;
	}

	void ModListForCalibrationTask::setLocalize(bool value)
	{
		privateLocalize = value;
	}

	std::wstring ModListForCalibrationTask::getFileName() const
	{
		return privateFileName;
	}
}
