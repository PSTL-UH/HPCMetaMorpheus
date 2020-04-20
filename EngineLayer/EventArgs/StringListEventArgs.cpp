#include "StringListEventArgs.h"
#include "stringhelper.h"

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

    bool StringListEventArgs::Equals( EventArgs *obj) const
    {
        StringListEventArgs*o = dynamic_cast<StringListEventArgs *>(obj);
        if ( o != nullptr ) {
            if (o->getStringList().size() != privateStringList.size() ) {
                return false;
            }
            for ( int i = 0; i < (int)privateStringList.size(); i++ ) {
                if ( o->getStringList()[i] != privateStringList[i] ) {
                    return false;
                }
            }
            return true;
        }
        
        return false;
    }
    
    int StringListEventArgs::GetHashCode() const
    {
        //just a temporary implementation to silence the compiler.
        return StringHelper::GetHashCode(privateStringList);            
    }
    
    std::string StringListEventArgs::ToString() const          
    {
        std::string s;
        for ( auto st : privateStringList ) {
            s += st;
        }

        return s;
    }
    
}
