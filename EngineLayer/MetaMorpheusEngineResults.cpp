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

	std::wstring MetaMorpheusEngineResults::ToString()
	{
		auto sb = new StringBuilder();
		sb->append(L"Time to run: " + Time);

		delete sb;
		return sb->toString();
	}
}
