#include "StringEventArgs.h"
#include "stringhelper.h"

namespace EngineLayer
{

    StringEventArgs::StringEventArgs(const std::string &s, std::vector<std::string> nestedIDs) : MyRecursiveEventArgs(nestedIDs)
    {
        this->privateS = s;
    }

    std::string StringEventArgs::getS() const
    {
        return privateS;
    }

    bool StringEventArgs::Equals( EventArgs *obj) const
    {
        StringEventArgs*o = dynamic_cast<StringEventArgs *>(obj);
        if ( o != nullptr ) { 
            return o->getS() == privateS;
        }
        return false;
    }
    
    int StringEventArgs::GetHashCode() const
    {
        return StringHelper::GetHashCode(privateS);            
    }
    
    std::string StringEventArgs::ToString() const          
    {
        return privateS;
    }
    
}
