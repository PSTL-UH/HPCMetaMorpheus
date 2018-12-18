#include "ModListForSearchTask.h"

namespace MetaMorpheusGUI
{

	ModListForSearchTask::ModListForSearchTask(const std::wstring &filePath)
	{
		this->FileName = filePath;
	}

	bool ModListForSearchTask::getFixed() const
	{
		return privateFixed;
	}

	void ModListForSearchTask::setFixed(bool value)
	{
		privateFixed = value;
	}

	bool ModListForSearchTask::getVariable() const
	{
		return privateVariable;
	}

	void ModListForSearchTask::setVariable(bool value)
	{
		privateVariable = value;
	}

	bool ModListForSearchTask::getLocalize() const
	{
		return privateLocalize;
	}

	void ModListForSearchTask::setLocalize(bool value)
	{
		privateLocalize = value;
	}

	bool ModListForSearchTask::getAlwaysKeep() const
	{
		return privateAlwaysKeep;
	}

	void ModListForSearchTask::setAlwaysKeep(bool value)
	{
		privateAlwaysKeep = value;
	}

	std::wstring ModListForSearchTask::getFileName() const
	{
		return privateFileName;
	}
}
