#include "MyRecursiveEventArgs.h"


namespace EngineLayer
{
    
    MyRecursiveEventArgs::MyRecursiveEventArgs(std::vector<std::string> &nestedIDs) : NestedIDs(nestedIDs)
    {
    }
}
