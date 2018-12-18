#include "ProgressEventArgs.h"


namespace EngineLayer
{

	ProgressEventArgs::ProgressEventArgs(int newProgress, const std::wstring &v, std::vector<std::wstring> &nestedIDs) : MyRecursiveEventArgs(nestedIDs)
	{
		NewProgress = newProgress;
		V = v;
	}
}
