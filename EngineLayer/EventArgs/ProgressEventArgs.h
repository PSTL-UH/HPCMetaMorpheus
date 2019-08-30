#pragma once

#include "MyRecursiveEventArgs.h"
#include <string>
#include <vector>


namespace EngineLayer
{
    class ProgressEventArgs : public MyRecursiveEventArgs

    {
    public:
        int NewProgress = 0;
        std::string V;
        
        ProgressEventArgs(int newProgress, const std::string &v, std::vector<std::string> &nestedIDs);

        bool Equals( EventArgs *obj) const override;
        
        int GetHashCode() const override;
        
        std::string ToString() const override;
    };
}
