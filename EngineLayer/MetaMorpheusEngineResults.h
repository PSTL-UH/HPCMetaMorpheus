#pragma once

#include <string>

#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
//namespace EngineLayer { class MetaMorpheusEngine; }
#include "MetaMorpheusEngine.h"

namespace EngineLayer
{
    class MetaMorpheusEngineResults
    {
    private:
        MetaMorpheusEngine *privateMyEngine;
        
    public:
        //TimeSpan Time;
        double Time;

        MetaMorpheusEngineResults(MetaMorpheusEngine *s);
        
        MetaMorpheusEngine *getMyEngine() const;
        
        std::string ToString();
    };
}
