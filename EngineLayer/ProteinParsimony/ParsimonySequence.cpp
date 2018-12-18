#include "ParsimonySequence.h"

using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
	namespace ProteinParsimony
	{

		ParsimonySequence::ParsimonySequence(PeptideWithSetModifications *pwsm, bool TreatModPeptidesAsDifferentPeptides)
		{
			Sequence = TreatModPeptidesAsDifferentPeptides ? pwsm->FullSequence : pwsm->BaseSequence;
			Protease = pwsm->DigestionParams.Protease;
		}

		std::wstring ParsimonySequence::getSequence() const
		{
			return privateSequence;
		}

		Protease *ParsimonySequence::getProtease() const
		{
			return privateProtease;
		}

		bool ParsimonySequence::Equals(std::any obj)
		{
			ParsimonySequence *other = std::any_cast<ParsimonySequence*>(obj);
			return other != nullptr && (getSequence() == L"" && other->getSequence() == L"" || getSequence() == other->getSequence()) && (getProtease() == nullptr && other->getProtease() == nullptr || getProtease()->Equals(other->getProtease()));
		}

		int ParsimonySequence::GetHashCode()
		{
			return getSequence().GetHashCode() ^ getProtease()->GetHashCode();
		}
	}
}
