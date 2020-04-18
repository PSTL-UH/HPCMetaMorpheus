#pragma once

#include <string>

#include "Proteomics/ProteolyticDigestion/PeptideWithSetModifications.h"
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    namespace ProteinParsimony
    {
        class ParsimonySequence
        {
        private:
            std::string privateSequence;
            Protease *privateProtease;
            
        public:
            ParsimonySequence(PeptideWithSetModifications *pwsm, bool TreatModPeptidesAsDifferentPeptides);
            
            std::string getSequence() const;
            Protease *getProtease() const;
            
            bool Equals(ParsimonySequence *other);
            
            int GetHashCode();
        };
    }
}
