#pragma once

#include <string>
#include "EventArgs.h"
//#include "../MetaMorpheusEngineResults.h"
namespace EngineLayer { class MetaMorpheusEngineResults; }

namespace EngineLayer
{
    class SingleEngineFinishedEventArgs : public EventArgs
    {
    public:
        MetaMorpheusEngineResults *const MyResults;
        
        virtual ~SingleEngineFinishedEventArgs()
        {
            //delete MyResults;
        }
        
        SingleEngineFinishedEventArgs(MetaMorpheusEngineResults *myResults);
        
        bool Equals( EventArgs *obj) const override;
    
        int GetHashCode() const override;
    
        std::string ToString() const override;
    };
}
