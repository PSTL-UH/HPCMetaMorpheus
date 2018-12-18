#include "StringListEventArgs.h"


namespace EngineLayer
{

	StringListEventArgs::StringListEventArgs(std::vector<std::wstring> &stringList)
	{
		StringList = stringList;
	}

	std::vector<std::wstring> StringListEventArgs::getStringList() const
	{
		return privateStringList;
	}
}
