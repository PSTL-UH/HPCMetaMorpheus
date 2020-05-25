#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include "../MetaMorpheusEngine.h"
#include "../MetaMorpheusEngineResults.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"

#include "ParsimonySequence.h"
#include "ProteinGroup.h"
#include "ProteinParsimonyResults.h"

using namespace EngineLayer;
using namespace EngineLayer::ProteinParsimony;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    class ProteinParsimonyEngine : public MetaMorpheusEngine
    {
        /// <summary>
        /// All peptides meeting the prefiltering criteria for parsimony (e.g., peptides from non-ambiguous
        //  high-confidence PSMs)
        /// </summary>
    private:
        std::unordered_set<PeptideWithSetModifications*> _fdrFilteredPeptides;
        std::vector<PeptideSpectralMatch*> _fdrFilteredPsms;
        std::vector<PeptideSpectralMatch*> _allPsms;
        static constexpr double FdrCutoffForParsimony = 0.01;
        
        /// <summary>
        /// User-selectable option that treats differently-modified forms of a peptide as different peptides
        /// for the purposes of parsimony
        /// </summary>
        const bool _treatModPeptidesAsDifferentPeptides;
        
    public:
        ProteinParsimonyEngine(std::vector<PeptideSpectralMatch*> &allPsms,
                               bool modPeptidesAreDifferent, CommonParameters *commonParameters,
                               std::vector<std::string> &nestedIds);


    protected:
        MetaMorpheusEngineResults* RunSpecific();

    private:
        /// <summary>
        /// TODO: Summarize parsimony;
        /// Parsimony algorithm based on: https://www.ncbi.nlm.nih.gov/pubmed/14632076 Anal Chem. 2003 Sep 1;75(17):4646-58.
        /// TODO: Note describing that peptide objects with the same sequence are associated with different proteins
        /// </summary>
        std::vector<ProteinGroup*> RunProteinParsimonyEngine();
                
        /// <summary>
        /// Builds protein groups after parsimony. Each protein group will only have 1 protein at this point. 
        /// Indistinguishable protein groups will be merged during scoring for computational efficiency reasons 
        /// (it's easier to tell when two protein groups are identical after they're scored)
        /// </summary>
        std::vector<ProteinGroup*> ConstructProteinGroups(std::unordered_set<PeptideWithSetModifications*> &uniquePeptides);
    };
}
