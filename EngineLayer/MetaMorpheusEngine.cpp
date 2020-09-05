#include "MetaMorpheusEngine.h"
#include "CommonParameters.h"
#include "Ms2ScanWithSpecificMass.h"
#include "MetaMorpheusEngineResults.h"
#include "EventArgs/StringEventArgs.h"
#include "EventArgs/ProgressEventArgs.h"
#include "EventArgs/SingleEngineEventArgs.h"
#include "EventArgs/SingleEngineFinishedEventArgs.h"

#include "../TaskLayer/MetaMorpheusTask.h"

#include <ctime>

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics::Fragmentation;

namespace EngineLayer
{

    std::unordered_map<DissociationType, double> MetaMorpheusEngine::complementaryIonConversionDictionary = std::unordered_map<DissociationType, double>
    {
        {DissociationType::HCD, Constants::protonMass},
	{DissociationType::ETD, 2 * Constants::protonMass},
	{DissociationType::CID, Constants::protonMass}
    };

    
    MetaMorpheusEngine::MetaMorpheusEngine(CommonParameters *cmnParameters, std::vector<std::string> nestedIDs,
        int verbosityLevel ) : nestedIds(nestedIDs)
    {
        commonParameters = new CommonParameters(cmnParameters);
        privateVerbosityLevel = verbosityLevel;
    }

    double MetaMorpheusEngine::CalculatePeptideScore(MsDataScan *thisScan, std::vector<MatchedFragmentIon*> &matchedFragmentIons,
                                                     double maximumMassThatFragmentIonScoreIsDoubled)
    {
        double score = 0;
        
        for (auto fragment : matchedFragmentIons)
        {
            double fragmentScore = 1 + (fragment->Intensity / thisScan->getTotalIonCurrent());
            score += fragmentScore;
            
            if (fragment->NeutralTheoreticalProduct->NeutralMass <= maximumMassThatFragmentIonScoreIsDoubled)
            {
                score += fragmentScore;
            }
        }
        
        return score;
    }
    
    std::vector<MatchedFragmentIon*> MetaMorpheusEngine::MatchFragmentIons(Ms2ScanWithSpecificMass *scan,
                                                                           std::vector<Product*> &theoreticalProducts,
                                                                           CommonParameters *commonParameters)
    {
        std::vector<MatchedFragmentIon*> matchedFragmentIons;
        
        // if the spectrum has no peaks
        if (scan->getExperimentalFragments().empty())
        {
            return matchedFragmentIons;
        }
        
        // search for ions in the spectrum
        for (auto product : theoreticalProducts)
        {
            // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
            if (std::isnan(product->NeutralMass))
            {
                continue;
            }
            
            // get the closest peak in the spectrum to the theoretical peak
            auto closestExperimentalMass = scan->GetClosestExperimentalFragmentMass(product->NeutralMass);
            
            // is the mass error acceptable?
            if (commonParameters->getProductMassTolerance()->Within(closestExperimentalMass->monoisotopicMass,
                                                                    product->NeutralMass) &&
                closestExperimentalMass->charge <= scan->getPrecursorCharge())
            {
                auto tempVar = new MatchedFragmentIon(product,
                                                      Chemistry::ClassExtensions::ToMz(closestExperimentalMass->monoisotopicMass,
                                                                                       closestExperimentalMass->charge),
                                                      std::get<1>(closestExperimentalMass->peaks.front()), //intensity
                                                      closestExperimentalMass->charge);
                matchedFragmentIons.push_back(tempVar);
            }
        }
        if (commonParameters->getAddCompIons())
        {
            DissociationType dtype = commonParameters->getDissociationType();
            double d = complementaryIonConversionDictionary[dtype];
            double protonMassShift = Chemistry::ClassExtensions::ToMass( d, 1);
            
            for (auto product : theoreticalProducts)
            {
                // unknown fragment mass or diagnostic ion or precursor; skip those
                if (std::isnan(product->NeutralMass)       ||
                    product->productType == ProductType::D ||
                    product->productType == ProductType::M)
                {
                    continue;
                }
                
                double compIonMass = scan->getPrecursorMass() + protonMassShift - product->NeutralMass;
                
                // get the closest peak in the spectrum to the theoretical peak
                auto closestExperimentalMass = scan->GetClosestExperimentalFragmentMass(compIonMass);
                
                // is the mass error acceptable?
                if (commonParameters->getProductMassTolerance()->Within(closestExperimentalMass->monoisotopicMass,
                                                                        compIonMass)                                &&
                    closestExperimentalMass->charge <= scan->getPrecursorCharge())
                {
                    auto tempVar2 = new MatchedFragmentIon (product,
                                                            Chemistry::ClassExtensions::ToMz(closestExperimentalMass->monoisotopicMass,
                                                                                             closestExperimentalMass->charge),
                                                            closestExperimentalMass->totalIntensity,
                                                            closestExperimentalMass->charge);
                    matchedFragmentIons.push_back(tempVar2);
                }
            }
        }
        
        return matchedFragmentIons;
    }
    
    MetaMorpheusEngineResults *MetaMorpheusEngine::Run()
    {
        StartingSingleEngine();
        time_t start, stop;
        time(&start);
        auto myResults = RunSpecific();
        time(&stop);
        myResults->Time = difftime( stop, start);
        FinishedSingleEngine(myResults);
        
        return myResults;
    }
    
    std::string MetaMorpheusEngine::GetId()
    {
        // return std::string::Join(",", nestedIds);
        std::string s = ",";
        for ( auto n: nestedIds ) {
            s += n;
        }
        return s;
    }
    
    void MetaMorpheusEngine::Warn(const std::string &v)
    {
        TaskLayer::MetaMorpheusTask::Warn(v);
    }

    void MetaMorpheusEngine::Status(const std::string &v)
    {
        TaskLayer::MetaMorpheusTask::Status("", v, privateVerbosityLevel);
    }
    
    void MetaMorpheusEngine::ReportProgress(ProgressEventArgs *v)
    {
        TaskLayer::MetaMorpheusTask::ReportProgress(v, privateVerbosityLevel);
    }

    void MetaMorpheusEngine::ReportEngineProgress(std::string key, int value)
    {
        TaskLayer::MetaMorpheusTask::ReportEngineProgress(key, value, privateVerbosityLevel);
    }
    
    void MetaMorpheusEngine::StartingSingleEngine()
    {
        TaskLayer::MetaMorpheusTask::StartingSingleEngine(nestedIds, privateVerbosityLevel);
    }
    
    void MetaMorpheusEngine::FinishedSingleEngine(MetaMorpheusEngineResults *myResults)
    {
        TaskLayer::MetaMorpheusTask::FinishedSingleEngine(nestedIds, myResults, privateVerbosityLevel);
    }
}
