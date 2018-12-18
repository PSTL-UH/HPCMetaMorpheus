#include "ModListForGPTMDTask.h"

namespace MetaMorpheusGUI
{

	ModListForGptmdTask::ModListForGptmdTask(const std::wstring &filePath)
	{
		this->FileName = filePath;
	}

	bool ModListForGptmdTask::getFixed() const
	{
		return privateFixed;
	}

	void ModListForGptmdTask::setFixed(bool value)
	{
		privateFixed = value;
	}

	bool ModListForGptmdTask::getVariable() const
	{
		return privateVariable;
	}

	void ModListForGptmdTask::setVariable(bool value)
	{
		privateVariable = value;
	}

	bool ModListForGptmdTask::getLocalize() const
	{
		return privateLocalize;
	}

	void ModListForGptmdTask::setLocalize(bool value)
	{
		privateLocalize = value;
	}

	bool ModListForGptmdTask::getGptmd() const
	{
		return privateGptmd;
	}

	void ModListForGptmdTask::setGptmd(bool value)
	{
		privateGptmd = value;
	}

	std::wstring ModListForGptmdTask::getFileName() const
	{
		return privateFileName;
	}
}
