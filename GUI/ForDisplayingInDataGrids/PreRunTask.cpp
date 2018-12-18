#include "PreRunTask.h"
#include "../../TaskLayer/MetaMorpheusTask.h"

using namespace TaskLayer;

namespace MetaMorpheusGUI
{

	PreRunTask::PreRunTask(MetaMorpheusTask *theTask) : metaMorpheusTask(theTask)
	{
	}

	std::wstring PreRunTask::getDisplayName() const
	{
		return privateDisplayName;
	}

	void PreRunTask::setDisplayName(const std::wstring &value)
	{
		privateDisplayName = value;
	}

	PreRunTask *PreRunTask::Clone()
	{
		return static_cast<PreRunTask*>(this->MemberwiseClone());
	}
}
