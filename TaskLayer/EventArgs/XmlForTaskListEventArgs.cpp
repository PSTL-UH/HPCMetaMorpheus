#include "XmlForTaskListEventArgs.h"
#include "../DbForTask.h"
#include "stringhelper.h"

namespace TaskLayer
{

    XmlForTaskListEventArgs::XmlForTaskListEventArgs(std::vector<DbForTask*> &newDatabases)
    {
        NewDatabases = newDatabases;
    }

    bool XmlForTaskListEventArgs::Equals( EventArgs *obj) const
    {
        XmlForTaskListEventArgs*o = dynamic_cast<XmlForTaskListEventArgs *>(obj);
        if ( o != nullptr ) { 
            if ( o->NewDatabases.size() != NewDatabases.size() ) {
                return false;
            }
            for ( int i=0; i <(int)NewDatabases.size(); i++ ) {
                if (  o->NewDatabases[i]->getFilePath() != NewDatabases[i]->getFilePath()   ||
                      o->NewDatabases[i]->getFileName() != NewDatabases[i]->getFileName()   ||
                      o->NewDatabases[i]->getIsContaminant() == NewDatabases[i]->getIsContaminant() ) {
                    return false;
                }
            }
            return true;
        }
        
        return false;
    }
    
    int XmlForTaskListEventArgs::GetHashCode() const
    {
        std::string s;
        for ( auto db : NewDatabases ) {
            s += db->getFilePath() + db->getFileName();
        }
        return StringHelper::GetHashCode(s);            
    }
    
    std::string XmlForTaskListEventArgs::ToString() const          
    {
        std::string s;
        for ( auto db : NewDatabases ) {
            s += db->getFilePath() + db->getFileName();
        }

        return s;
    }
}
