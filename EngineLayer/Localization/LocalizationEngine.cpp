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
        
        LocalizationEngine::LocalizationEngine(std::vector<PeptideSpectralMatch*> &allResultingIdentifications,
                                               MsDataFile *myMsDataFile, CommonParameters *commonParameters,
                                               std::vector<std::string> nestedIds, int verbosityLevel) :
            MetaMorpheusEngine(commonParameters, nestedIds, verbosityLevel ),
            AllResultingIdentifications(allResultingIdentifications),
            MyMsDataFile(myMsDataFile)
        {
        }
        
        MetaMorpheusEngineResults *LocalizationEngine::RunSpecific()
        {
            // don't try to localize mass differences for ambiguous peptides
#ifdef ORIG
            //  for (PeptideSpectralMatch *psm : AllResultingIdentifications.Where([&] (std::any b)  {
            //          return b::FullSequence != nullptr;
            //      }))
#endif
            for (PeptideSpectralMatch *psm : AllResultingIdentifications ) 
            {
                if ( psm->getFullSequence().length() == 0 ) {
                    continue;
                }
                // Stop loop if canceled
                if (GlobalVariables::getStopLoops())
                {
                    break;
                }
                
                MsDataScan *scan = MyMsDataFile->GetOneBasedScan(psm->getScanNumber());
                std::vector<MassSpectrometry::IsotopicEnvelope*> v1;
                auto scanWithSpecificMass = new Ms2ScanWithSpecificMass(scan,
                                                                        psm->getScanPrecursorMonoisotopicPeakMz(),
                                                                        psm->getScanPrecursorCharge(),
                                                                        psm->getFullFilePath(),
                                                                        commonParameters,
                                                                        v1);
                PeptideWithSetModifications *peptide = std::get<1>(psm->getBestMatchingPeptides().front());
                double massDifference = psm->getScanPrecursorMass() - peptide->getMonoisotopicMass();

                // this section will iterate through all residues of the peptide and try to localize the mass-diff
                //at each residue and report a score for each residue
                std::vector<double> localizedScores;
                for (int r = 0; r < peptide->getLength(); r++)
                {
                    // create new PeptideWithSetMods with unidentified mass difference at the given residue
                    PeptideWithSetModifications *peptideWithLocalizedMassDiff = peptide->Localize(r, massDifference);

                    // this is the list of theoretical products for this peptide with mass-difference on this residue
                    std::vector<Product*> productsWithLocalizedMassDiff = peptideWithLocalizedMassDiff->Fragment(
                        commonParameters->getDissociationType(),
                        commonParameters->getDigestionParams()->getFragmentationTerminus());
                        
                    
                    auto matchedIons = MatchFragmentIons(scanWithSpecificMass, productsWithLocalizedMassDiff, commonParameters);
                        
                    // score when the mass-diff is on this residue
                    double localizedScore = CalculatePeptideScore(scan, matchedIons, 0);
                    
                    localizedScores.push_back(localizedScore);
                }
                
                psm->setLocalizedScores(localizedScores);
                
                delete scanWithSpecificMass;
            }
            
            return new LocalizationEngineResults(this);
        }
    }
}
