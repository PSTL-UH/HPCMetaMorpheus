#pragma once

#include <string>
#include <vector>

#include "EventArgs.h"

namespace EngineLayer
{
	class StringListEventArgs : public EventArgs
	{
	private:
		std::vector<std::string> privateStringList;

	public:
		StringListEventArgs(std::vector<std::string> &stringList);

		std::vector<std::string> getStringList() const;
	};
}
