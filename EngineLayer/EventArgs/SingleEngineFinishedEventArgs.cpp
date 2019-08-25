#include "SingleEngineFinishedEventArgs.h"
#include "../MetaMorpheusEngineResults.h"


namespace EngineLayer
{

	SingleEngineFinishedEventArgs::SingleEngineFinishedEventArgs(MetaMorpheusEngineResults *myResults) : MyResults(myResults)
	{
	}

	std::string SingleEngineFinishedEventArgs::ToString()
	{
		return MyResults->ToString();
	}
}
