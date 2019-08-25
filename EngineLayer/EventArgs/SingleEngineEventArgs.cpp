#include "SingleEngineEventArgs.h"
#include "../MetaMorpheusEngine.h"


namespace EngineLayer
{

	SingleEngineEventArgs::SingleEngineEventArgs(MetaMorpheusEngine *myEngine)
	{
		setMyEngine(myEngine);
	}

	MetaMorpheusEngine *SingleEngineEventArgs::getMyEngine() const
	{
		return privateMyEngine;
	}

	void SingleEngineEventArgs::setMyEngine(MetaMorpheusEngine *value)
	{
		privateMyEngine = value;
	}
}
