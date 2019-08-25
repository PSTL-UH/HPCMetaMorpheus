#include "ProteinParsimonyEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"

using namespace EngineLayer::ProteinParsimony;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

	ProteinParsimonyEngine::ProteinParsimonyEngine(std::vector<PeptideSpectralMatch*> &allPsms, bool modPeptidesAreDifferent, CommonParameters *commonParameters, std::vector<std::string> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds), _fdrFilteredPeptides(std::unordered_set<PeptideWithSetModifications*>()), _allPsms(allPsms), _treatModPeptidesAsDifferentPeptides(modPeptidesAreDifferent)
	{

		if (!allPsms.Any())
		{
			_fdrFilteredPsms = std::vector<PeptideSpectralMatch*>();
		}

		// parsimony will only use non-ambiguous, high-confidence PSMs
		// KEEP decoys and contaminants for use in parsimony!
		if (modPeptidesAreDifferent)
		{
			_fdrFilteredPsms = allPsms.Where([&] (std::any p)
			{
				return p::FullSequence != nullptr && p::FdrInfo::QValue <= FdrCutoffForParsimony && p::FdrInfo::QValueNotch <= FdrCutoffForParsimony;
			}).ToList();
		}
		else
		{
			_fdrFilteredPsms = allPsms.Where([&] (std::any p)
			{
				return p::BaseSequence != nullptr && p::FdrInfo::QValue <= FdrCutoffForParsimony && p::FdrInfo::QValueNotch <= FdrCutoffForParsimony;
			}).ToList();
		}

		// if PSM is a decoy, add only decoy sequences; same for contaminants
		// peptides to use in parsimony = peptides observed in high-confidence PSMs

		for (auto psm : _fdrFilteredPsms)
		{
			if (psm->getIsDecoy())
			{
				for (auto peptide : psm->BestMatchingPeptides->Select([&] (std::any p)
				{
					p::Peptide;
				}).Where([&] (std::any p)
				{
					p::Protein::IsDecoy;
				}))
				{
				}
			}
		}
			else if (psm::IsContaminant)
			{
				for (auto peptide : psm::BestMatchingPeptides->Select([&] (std::any p)
				{
					p::Peptide;
				}).Where([&] (std::any p)
				{
					p::Protein::IsContaminant;
				}))
				{
				}
			}
	}
}
