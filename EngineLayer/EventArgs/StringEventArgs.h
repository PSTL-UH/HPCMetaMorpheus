#pragma once

#include "MyRecursiveEventArgs.h"
#include <string>
#include <vector>


namespace EngineLayer
{
	class StringEventArgs : public MyRecursiveEventArgs
	{
	private:
		std::wstring privateS;

	public:
		StringEventArgs(const std::wstring &s, std::vector<std::wstring> &nestedIDs);

		std::wstring getS() const;
	};
}
