#pragma once

#include <string>
#include <vector>

#include "EventArgs.h"

namespace EngineLayer
{
    class MyRecursiveEventArgs : public EventArgs
    {
    public:
        const std::vector<std::string> NestedIDs;
        
        MyRecursiveEventArgs(std::vector<std::string> &nestedIDs);
    };
}
