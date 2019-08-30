#include "ProgressEventArgs.h"
#include "stringhelper.h"

namespace EngineLayer
{

    ProgressEventArgs::ProgressEventArgs(int newProgress, const std::string &v, std::vector<std::string> &nestedIDs) : MyRecursiveEventArgs(nestedIDs)
    {
        NewProgress = newProgress;
        V = v;
    }
    bool ProgressEventArgs::Equals( EventArgs *obj) const
    {
        //just a temporary implementation to silence the compiler.
        ProgressEventArgs *o = dynamic_cast<ProgressEventArgs *>(obj);
        return o->V == V;


    }
    int ProgressEventArgs::GetHashCode() const
    {
        //just a temporary implementation to silence the compiler.
        return StringHelper::GetHashCode(V);            

    }
    
    std::string ProgressEventArgs::ToString() const
    {
        return V;
    }
    
}
