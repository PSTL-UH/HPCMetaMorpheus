﻿#pragma once

#include "../MetaMorpheusEngineResults.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MetaMorpheusEngine; }

namespace EngineLayer
{
	namespace Localization
	{
		class LocalizationEngineResults : public MetaMorpheusEngineResults
		{
		public:
			LocalizationEngineResults(MetaMorpheusEngine *s);
		};
	}
}