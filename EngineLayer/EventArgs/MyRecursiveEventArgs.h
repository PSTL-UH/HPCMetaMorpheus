#pragma once

#include <string>
#include <vector>


namespace EngineLayer
{
	class MyRecursiveEventArgs : public EventArgs
	{
	public:
		const std::vector<std::wstring> NestedIDs;

		MyRecursiveEventArgs(std::vector<std::string> &nestedIDs);
	};
}
