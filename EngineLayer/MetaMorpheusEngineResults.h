#pragma once

#include <string>
#include <sys/time.h>

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
        //TimeSpan Time;
        struct timeval Time;

        MetaMorpheusEngineResults(MetaMorpheusEngine *s);
        
        MetaMorpheusEngine *getMyEngine() const;
        
        std::string ToString();
    };
}
