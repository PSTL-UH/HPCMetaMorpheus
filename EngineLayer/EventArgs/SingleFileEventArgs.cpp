#include "SingleFileEventArgs.h"
#include "stringhelper.h"


namespace EngineLayer
{

    SingleFileEventArgs::SingleFileEventArgs(const std::string &writtenFile, std::vector<std::string> &nestedIds) : MyRecursiveEventArgs(nestedIds)
    {
        setWrittenFile(writtenFile);
    }
    
    std::string SingleFileEventArgs::getWrittenFile() const
    {
        return privateWrittenFile;
    }

    void SingleFileEventArgs::setWrittenFile(const std::string &value)
    {
        privateWrittenFile = value;
    }

    bool SingleFileEventArgs::Equals( EventArgs *obj) const
    {
        SingleFileEventArgs*o = dynamic_cast<SingleFileEventArgs *>(obj);
        if ( o != nullptr ) { 
            return o->getWrittenFile() == privateWrittenFile;
        }
        return false;
    }
    
    int SingleFileEventArgs::GetHashCode() const
    {
        return StringHelper::GetHashCode(privateWrittenFile);            
    }
    
    std::string SingleFileEventArgs::ToString() const          
    {
        return privateWrittenFile;
    }
}
