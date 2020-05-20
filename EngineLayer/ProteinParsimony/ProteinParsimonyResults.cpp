#include "ProteinParsimonyResults.h"
#include "ProteinGroup.h"
#include "ProteinParsimonyEngine.h"


namespace EngineLayer
{

    ProteinParsimonyResults::ProteinParsimonyResults(ProteinParsimonyEngine *proteinAnalysisEngine) :
        MetaMorpheusEngineResults(proteinAnalysisEngine)
    {
    }
    
    std::vector<ProteinGroup*> ProteinParsimonyResults::getProteinGroups() const
    {
        return privateProteinGroups;
    }
    
    void ProteinParsimonyResults::setProteinGroups(const std::vector<ProteinGroup*> &value)
    {
        privateProteinGroups = value;
    }
    
    std::string ProteinParsimonyResults::ToString()
    {
        auto sb = new StringBuilder();
        sb->appendLine(MetaMorpheusEngineResults::ToString());
        
        std::string s= sb->toString();
        delete sb;
        return s;
    }
}
