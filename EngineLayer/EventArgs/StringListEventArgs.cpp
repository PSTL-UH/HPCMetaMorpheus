#include "StringListEventArgs.h"


namespace EngineLayer
{

    StringListEventArgs::StringListEventArgs(std::vector<std::string> &stringList)
    {
        privateStringList = stringList;
    }
    
    std::vector<std::string> StringListEventArgs::getStringList() const
    {
        return privateStringList;
    }
}
