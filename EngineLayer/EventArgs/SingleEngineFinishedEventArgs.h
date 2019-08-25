#pragma once

#include <string>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MetaMorpheusEngineResults; }


namespace EngineLayer
{
	class SingleEngineFinishedEventArgs : public EventArgs
	{
	public:
		MetaMorpheusEngineResults *const MyResults;

		virtual ~SingleEngineFinishedEventArgs()
		{
			delete MyResults;
		}

		SingleEngineFinishedEventArgs(MetaMorpheusEngineResults *myResults);

		std::string ToString();
	};
}
