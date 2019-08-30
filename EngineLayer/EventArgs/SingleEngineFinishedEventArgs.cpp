#include "SingleEngineFinishedEventArgs.h"
#include "../MetaMorpheusEngineResults.h"
#include "stringhelper.h"

namespace EngineLayer
{

    SingleEngineFinishedEventArgs::SingleEngineFinishedEventArgs(MetaMorpheusEngineResults *myResults) : MyResults(myResults)
    {
    }
    
    bool SingleEngineFinishedEventArgs::Equals( EventArgs *obj) const
    {
        //just a temporary implementation to silence the compiler.
        SingleEngineFinishedEventArgs *o = dynamic_cast<SingleEngineFinishedEventArgs *>(obj);
        return o->MyResults == MyResults;
    }
    
    int SingleEngineFinishedEventArgs::GetHashCode() const
    {
        //just a temporary implementation to silence the compiler.
        return StringHelper::GetHashCode(MyResults->ToString());            
    }
    
    std::string SingleEngineFinishedEventArgs::ToString() const
    {
        return MyResults->ToString();
    }
}
