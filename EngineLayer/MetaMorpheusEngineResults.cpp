#include "MetaMorpheusEngineResults.h"
#include "MetaMorpheusEngine.h"


namespace EngineLayer
{

	MetaMorpheusEngineResults::MetaMorpheusEngineResults(MetaMorpheusEngine *s)
	{
		MyEngine = s;
	}

	MetaMorpheusEngine *MetaMorpheusEngineResults::getMyEngine() const
	{
		return privateMyEngine;
	}

	std::string MetaMorpheusEngineResults::ToString()
	{
		auto sb = new StringBuilder();
		sb->append("Time to run: " + Time);

		delete sb;
		return sb->toString();
	}
}
