#include "MetaMorpheusException.h"


namespace EngineLayer
{
    
    MetaMorpheusException::MetaMorpheusException(const std::string &message) : runtime_error(message)
    {
    }
}
