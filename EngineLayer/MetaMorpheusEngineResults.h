#pragma once

#include <string>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MetaMorpheusEngine; }


namespace EngineLayer
{
	class MetaMorpheusEngineResults
	{
	private:
		MetaMorpheusEngine *privateMyEngine;

	public:
		TimeSpan Time;

		MetaMorpheusEngineResults(MetaMorpheusEngine *s);

		MetaMorpheusEngine *getMyEngine() const;

		std::wstring ToString() override;
	};
}
