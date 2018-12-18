#pragma once

#include "ForTreeView.h"
#include <string>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class MetaMorpheusTask; }

using namespace TaskLayer;

namespace MetaMorpheusGUI
{
	class InRunTask : public ForTreeView
	{
	public:
		MetaMorpheusTask *const Task;

		virtual ~InRunTask()
		{
			delete Task;
		}

		InRunTask(const std::wstring &displayName, MetaMorpheusTask *task);
	};
}
