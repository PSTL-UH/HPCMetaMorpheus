#pragma once

#include "../MetaMorpheusEngineResults.h"
#include <string>
#include <vector>
#include "stringbuilder.h"

#include "ProteinGroup.h"
namespace EngineLayer { class ProteinParsimonyEngine; }


namespace EngineLayer
{
    class ProteinParsimonyResults : public MetaMorpheusEngineResults
    {
    private:
        std::vector<ProteinGroup*> privateProteinGroups;
        
    public:
        ProteinParsimonyResults(ProteinParsimonyEngine *proteinAnalysisEngine);
        
        std::vector<ProteinGroup*> getProteinGroups() const;
        void setProteinGroups(const std::vector<ProteinGroup*> &value);
        
        std::string ToString();
    };
}
