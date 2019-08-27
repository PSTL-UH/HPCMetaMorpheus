#pragma once

#include <string>
#include "EventArgs.h"
#include "../MetaMorpheusEngineResults.h"


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
