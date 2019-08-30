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
            //just a temporary implementation to silence the compiler.
            StringEventArgs*o = dynamic_cast<StringEventArgs *>(obj);
            return o->getS() == privateS;
        }
    
        int StringEventArgs::GetHashCode() const
        {
            //just a temporary implementation to silence the compiler.
            return StringHelper::GetHashCode(privateS);            
        }
    
        std::string StringEventArgs::ToString() const          
        {
            //just a temporary implementation to silence the compiler.
            return privateS;
        }

}
