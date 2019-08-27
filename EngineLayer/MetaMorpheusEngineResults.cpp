#include "MetaMorpheusEngineResults.h"
#include "MetaMorpheusEngine.h"


namespace EngineLayer
{

    MetaMorpheusEngineResults::MetaMorpheusEngineResults(MetaMorpheusEngine *s)
    {
        privateMyEngine = s;
    }
    
    MetaMorpheusEngine *MetaMorpheusEngineResults::getMyEngine() const
    {
        return privateMyEngine;
    }
    
    std::string MetaMorpheusEngineResults::ToString()
	{
            auto sb = new StringBuilder();
            sb->append("Time to run: " + Time);
            
            std::string s =  sb->toString();
            delete sb;
            return s;
	}
}
