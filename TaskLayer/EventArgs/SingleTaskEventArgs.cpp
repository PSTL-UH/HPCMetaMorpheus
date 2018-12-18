#include "SingleTaskEventArgs.h"


namespace TaskLayer
{

	SingleTaskEventArgs::SingleTaskEventArgs(const std::wstring &displayName)
	{
		this->setDisplayName(displayName);
	}

	std::wstring SingleTaskEventArgs::getDisplayName() const
	{
		return privateDisplayName;
	}

	void SingleTaskEventArgs::setDisplayName(const std::wstring &value)
	{
		privateDisplayName = value;
	}
}
