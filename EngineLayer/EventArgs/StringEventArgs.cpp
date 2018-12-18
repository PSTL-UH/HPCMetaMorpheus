#include "StringEventArgs.h"


namespace EngineLayer
{

	StringEventArgs::StringEventArgs(const std::wstring &s, std::vector<std::wstring> &nestedIDs) : MyRecursiveEventArgs(nestedIDs)
	{
		this->S = s;
	}

	std::wstring StringEventArgs::getS() const
	{
		return privateS;
	}
}
