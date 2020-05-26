#pragma once

#include "../MetaMorpheusEngineResults.h"

#include "../MetaMorpheusEngine.h"
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
