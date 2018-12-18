#include "SingleEngineFinishedEventArgs.h"
#include "../MetaMorpheusEngineResults.h"


namespace EngineLayer
{

	SingleEngineFinishedEventArgs::SingleEngineFinishedEventArgs(MetaMorpheusEngineResults *myResults) : MyResults(myResults)
	{
	}

	std::wstring SingleEngineFinishedEventArgs::ToString()
	{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		return MyResults->ToString();
	}
}
