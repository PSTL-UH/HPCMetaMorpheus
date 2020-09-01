#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "stringhelper.h"


#include "../PeptideSpectralMatch.h"
#include "../ProteinParsimony/ProteinGroup.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"

using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    class ProteinScoringAndFdrEngine : public MetaMorpheusEngine
    {
    private:
        std::vector<PeptideSpectralMatch*> NewPsms;
        const bool NoOneHitWonders;
        const bool TreatModPeptidesAsDifferentPeptides;
        const bool MergeIndistinguishableProteinGroups;
        std::vector<ProteinGroup*> ProteinGroups;
        
    public:
        ProteinScoringAndFdrEngine(std::vector<ProteinGroup*> &proteinGroups, std::vector<PeptideSpectralMatch*> &newPsms,
                                   bool noOneHitWonders, bool treatModPeptidesAsDifferentPeptides,
                                   bool mergeIndistinguishableProteinGroups, CommonParameters *commonParameters,
                                   std::vector<std::string> &nestedIds, int verbosityLevel=0);
        
    protected:
        MetaMorpheusEngineResults *RunSpecific() override;
        
    private:
        static std::string StripDecoyIdentifier(const std::string &proteinGroupName); //we're keeping only the better scoring protein group for each target/decoy pair. to do that we need to strip decoy from the name temporarily. this is the "top-picked" method
        
        void ScoreProteinGroups(std::vector<ProteinGroup*> &proteinGroups, std::vector<PeptideSpectralMatch*> &psmList);
        
        std::vector<ProteinGroup*> DoProteinFdr(std::vector<ProteinGroup*> &proteinGroups);
    };
}
