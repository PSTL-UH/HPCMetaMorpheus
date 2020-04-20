#pragma once

#include "EventArgs.h"
#include <vector>

#include "../DbForTask.h"

namespace TaskLayer
{
    class XmlForTaskListEventArgs : public EventArgs
    {
    public:
        std::vector<DbForTask*> NewDatabases;
        
        XmlForTaskListEventArgs(std::vector<DbForTask*> &newDatabases);
        
        bool Equals( EventArgs *obj) const override;
    
        int GetHashCode() const override;
    
        std::string ToString() const override;

    };
}
