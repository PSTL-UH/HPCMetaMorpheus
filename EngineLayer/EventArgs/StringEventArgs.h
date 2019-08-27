#pragma once

#include "MyRecursiveEventArgs.h"
#include <string>
#include <vector>


namespace EngineLayer
{
    class StringEventArgs : public MyRecursiveEventArgs
    {
    private:
        std::string privateS;
        
	public:
        StringEventArgs(const std::string &s, std::vector<std::string> &nestedIDs);
        
        std::string getS() const;
    };
}
