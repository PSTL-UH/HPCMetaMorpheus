#include "ParsimonySequence.h"
#include "stringhelper.h"

using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
    namespace ProteinParsimony
    {
        
        ParsimonySequence::ParsimonySequence(PeptideWithSetModifications *pwsm, bool TreatModPeptidesAsDifferentPeptides)
        {
            privateSequence = TreatModPeptidesAsDifferentPeptides ? pwsm->getFullSequence() : pwsm->getBaseSequence();
            privateProtease = pwsm->getDigestionParams()->getProtease();
        }
        
        std::string ParsimonySequence::getSequence() const
        {
            return privateSequence;
        }
        
        Protease *ParsimonySequence::getProtease() const
        {
            return privateProtease;
        }
        
        bool ParsimonySequence::Equals(ParsimonySequence *other)
        {
            return other != nullptr &&
                ( (getSequence() == "" && other->getSequence() == "") ||
                 getSequence() == other->getSequence()              ) &&
                ( (getProtease() == nullptr && other->getProtease() == nullptr) ||
                 getProtease()->Equals(other->getProtease()));
        }
        
        int ParsimonySequence::GetHashCode()
        {
            return StringHelper::GetHashCode(getSequence()) ^ getProtease()->GetHashCode();
        }
    }
}
