#include "MassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"


namespace EngineLayer
{

    MassDiffAcceptor::MassDiffAcceptor(const std::string &fileNameAddition)
    {
        privateFileNameAddition = fileNameAddition;
        setNumNotches(1);
    }
    
    int MassDiffAcceptor::getNumNotches() const
    {
        return privateNumNotches;
    }
    
    void MassDiffAcceptor::setNumNotches(int value)
    {
        privateNumNotches = value;
    }
    
    std::string MassDiffAcceptor::getFileNameAddition() const
    {
        return privateFileNameAddition;
    }
}
