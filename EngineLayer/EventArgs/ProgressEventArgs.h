#pragma once

#include "MyRecursiveEventArgs.h"
#include <string>
#include <vector>


namespace EngineLayer
{
	class ProgressEventArgs : public MyRecursiveEventArgs

	{
	public:
		int NewProgress = 0;
		std::wstring V;

		ProgressEventArgs(int newProgress, const std::wstring &v, std::vector<std::wstring> &nestedIDs);
	};
}
