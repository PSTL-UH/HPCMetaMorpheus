#include "SingleEngineEventArgs.h"
#include "../MetaMorpheusEngine.h"
#include "stringhelper.h"

namespace EngineLayer
{

	SingleEngineEventArgs::SingleEngineEventArgs(MetaMorpheusEngine *myEngine)
	{
		setMyEngine(myEngine);
	}

	MetaMorpheusEngine *SingleEngineEventArgs::getMyEngine() const
	{
		return privateMyEngine;
	}

	void SingleEngineEventArgs::setMyEngine(MetaMorpheusEngine *value)
	{
		privateMyEngine = value;
	}

        bool SingleEngineEventArgs::Equals( EventArgs *obj) const
        {
            //just a temporary implementation to silence the compiler.
            SingleEngineEventArgs *o = dynamic_cast<SingleEngineEventArgs *>(obj);
            return o->getMyEngine()->GetId() == privateMyEngine->GetId();
        }
    
        int SingleEngineEventArgs::GetHashCode() const
        {
            //just a temporary implementation to silence the compiler.
            return StringHelper::GetHashCode(privateMyEngine->GetId());            
        }
    
        std::string SingleEngineEventArgs::ToString() const          
        {
            //just a temporary implementation to silence the compiler.
            return privateMyEngine->GetId();
        }

}
