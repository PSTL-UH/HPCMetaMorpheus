#pragma once

#include "EventArgs.h"
#include <string>


namespace TaskLayer
{
    class SingleTaskEventArgs : public EventArgs
    {
    private:
        std::string privateDisplayName;
        
    public:
        SingleTaskEventArgs(const std::string &displayName);
        
        std::string getDisplayName() const;
        void setDisplayName(const std::string &value);
    };
}
