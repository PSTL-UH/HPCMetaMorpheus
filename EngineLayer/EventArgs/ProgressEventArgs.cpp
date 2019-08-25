#include "ProgressEventArgs.h"


namespace EngineLayer
{

	ProgressEventArgs::ProgressEventArgs(int newProgress, const std::string &v, std::vector<std::string> &nestedIDs) : MyRecursiveEventArgs(nestedIDs)
	{
		NewProgress = newProgress;
		V = v;
	}
}
