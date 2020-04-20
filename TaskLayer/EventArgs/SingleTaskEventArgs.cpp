#include "SingleTaskEventArgs.h"
#include "stringhelper.h"

namespace TaskLayer
{

    SingleTaskEventArgs::SingleTaskEventArgs(const std::string &displayName)
    {
        this->setDisplayName(displayName);
    }
    
    std::string SingleTaskEventArgs::getDisplayName() const
    {
        return privateDisplayName;
    }
    
    void SingleTaskEventArgs::setDisplayName(const std::string &value)
    {
        privateDisplayName = value;
    }

    bool SingleTaskEventArgs::Equals( EventArgs *obj) const
    {
        SingleTaskEventArgs*o = dynamic_cast<SingleTaskEventArgs *>(obj);
        if ( o != nullptr ) { 
            return o->getDisplayName() == privateDisplayName;
        }
        return false;
    }
    
    int SingleTaskEventArgs::GetHashCode() const
    {
        return StringHelper::GetHashCode(privateDisplayName);            
    }
    
    std::string SingleTaskEventArgs::ToString() const          
    {
        return privateDisplayName;
    }
}
