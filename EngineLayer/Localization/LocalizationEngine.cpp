#include "LocalizationEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "../GlobalVariables.h"
#include "../Ms2ScanWithSpecificMass.h"
#include "LocalizationEngineResults.h"

using namespace MassSpectrometry;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
	namespace Localization
	{

		LocalizationEngine::LocalizationEngine(std::vector<PeptideSpectralMatch*> &allResultingIdentifications, MsDataFile *myMsDataFile, CommonParameters *commonParameters, std::vector<std::string> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds), AllResultingIdentifications(allResultingIdentifications), MyMsDataFile(myMsDataFile)
		{
		}

		MetaMorpheusEngineResults *LocalizationEngine::RunSpecific()
		{
			// don't try to localize mass differences for ambiguous peptides
			for (PeptideSpectralMatch *psm : AllResultingIdentifications.Where([&] (std::any b)
			{
				return b::FullSequence != nullptr;
			}))
			{
				// Stop loop if canceled
				if (GlobalVariables::getStopLoops())
				{
					break;
				}

				MsDataScan *scan = MyMsDataFile->GetOneBasedScan(psm::ScanNumber);
				Ms2ScanWithSpecificMass *scanWithSpecificMass = new Ms2ScanWithSpecificMass(scan, psm::ScanPrecursorMonoisotopicPeakMz, psm::ScanPrecursorCharge, psm::FullFilePath, commonParameters);
				PeptideWithSetModifications *peptide = psm::BestMatchingPeptides::First().Peptide;
				double massDifference = psm::ScanPrecursorMass - peptide->MonoisotopicMass;

				// this section will iterate through all residues of the peptide and try to localize the mass-diff at each residue and report a score for each residue
				auto localizedScores = std::vector<double>();
				for (int r = 0; r < peptide->Length; r++)
				{
					// create new PeptideWithSetMods with unidentified mass difference at the given residue
					PeptideWithSetModifications *peptideWithLocalizedMassDiff = peptide->Localize(r, massDifference);

					// this is the list of theoretical products for this peptide with mass-difference on this residue
					std::vector<Product*> productsWithLocalizedMassDiff = peptideWithLocalizedMassDiff->Fragment(commonParameters->getDissociationType(), commonParameters->getDigestionParams()->FragmentationTerminus).ToList();

					auto matchedIons = MatchFragmentIons(scanWithSpecificMass, productsWithLocalizedMassDiff, commonParameters);

					// score when the mass-diff is on this residue
					double localizedScore = CalculatePeptideScore(scan, matchedIons, 0);

					localizedScores.push_back(localizedScore);
				}

				psm->LocalizedScores = localizedScores;

//C# TO C++ CONVERTER TODO TASK: A 'delete scanWithSpecificMass' statement was not added since scanWithSpecificMass was passed to a method or constructor. Handle memory management manually.
			}

			return new LocalizationEngineResults(this);
		}
	}
}
