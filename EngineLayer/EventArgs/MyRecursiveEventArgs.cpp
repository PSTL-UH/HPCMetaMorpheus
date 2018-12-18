#include "MyRecursiveEventArgs.h"


namespace EngineLayer
{

	MyRecursiveEventArgs::MyRecursiveEventArgs(std::vector<std::wstring> &nestedIDs) : NestedIDs(nestedIDs)
	{
	}
}
