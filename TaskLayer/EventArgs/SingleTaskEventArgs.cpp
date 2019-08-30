#include "SingleTaskEventArgs.h"


namespace TaskLayer
{

	SingleTaskEventArgs::SingleTaskEventArgs(const std::string &displayName)
	{
		this->setDisplayName(displayName);
	}

	std::string SingleTaskEventArgs::getDisplayName() const
	{
		return privateDisplayName;
	}

	void SingleTaskEventArgs::setDisplayName(const std::string &value)
	{
		privateDisplayName = value;
	}
}
