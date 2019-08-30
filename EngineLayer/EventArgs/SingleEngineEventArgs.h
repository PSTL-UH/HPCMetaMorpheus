#pragma once

#include <string>
#include "EventArgs.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MetaMorpheusEngine; }


namespace EngineLayer
{
    class SingleEngineEventArgs : public EventArgs
    {
    private:
        MetaMorpheusEngine *privateMyEngine;
        
    public:
        SingleEngineEventArgs(MetaMorpheusEngine *myEngine);
        
        MetaMorpheusEngine *getMyEngine() const;
        void setMyEngine(MetaMorpheusEngine *value);

        bool Equals( EventArgs *obj) const override;
    
        int GetHashCode() const override;
    
        std::string ToString() const override; 
    };
}
