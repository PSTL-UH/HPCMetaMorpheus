#include "InRunTask.h"
#include "../../TaskLayer/MetaMorpheusTask.h"

using namespace TaskLayer;

namespace MetaMorpheusGUI
{

	InRunTask::InRunTask(const std::wstring &displayName, MetaMorpheusTask *task) : ForTreeView(displayName, displayName), Task(task)
	{
	}
}
