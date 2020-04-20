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

        bool Equals( EventArgs *obj) const override;
    
        int GetHashCode() const override;
    
        std::string ToString() const override;

    };
}
