#include "MetaMorpheusEngine.h"
#include "CommonParameters.h"
#include "Ms2ScanWithSpecificMass.h"
#include "MetaMorpheusEngineResults.h"
#include "EventArgs/StringEventArgs.h"
#include "EventArgs/ProgressEventArgs.h"
#include "EventArgs/SingleEngineEventArgs.h"
#include "EventArgs/SingleEngineFinishedEventArgs.h"

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


    EventHandler<SingleEngineEventArgs> *StartingSingleEngineHandler = new EventHandler<SingleEngineEventArgs>();

    EventHandler<SingleEngineFinishedEventArgs> *FinishedSingleEngineHandler = new EventHandler<SingleEngineFinishedEventArgs>();

    EventHandler<StringEventArgs> *OutLabelStatusHandler = new EventHandler<StringEventArgs>();

    EventHandler<StringEventArgs> *WarnHandler = new EventHandler<StringEventArgs>();

    EventHandler<ProgressEventArgs> *OutProgressHandler = new EventHandler<ProgressEventArgs>();


    
    MetaMorpheusEngine::MetaMorpheusEngine(CommonParameters *commonParameters, std::vector<std::string> &nestedIds) : commonParameters(commonParameters), nestedIds(nestedIds)
    {
    }

    double MetaMorpheusEngine::CalculatePeptideScore(MsDataScan *thisScan, std::vector<MatchedFragmentIon*> &matchedFragmentIons, double maximumMassThatFragmentIonScoreIsDoubled)
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
    
    std::vector<MatchedFragmentIon*> MetaMorpheusEngine::MatchFragmentIons(Ms2ScanWithSpecificMass *scan, std::vector<Product*> &theoreticalProducts, CommonParameters *commonParameters)
    {
        auto matchedFragmentIons = std::vector<MatchedFragmentIon*>();
        
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
            if (commonParameters->getProductMassTolerance()->Within(closestExperimentalMass->monoisotopicMass, product->NeutralMass) && closestExperimentalMass->charge <= scan->getPrecursorCharge())
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
                    MatchedFragmentIon tempVar2(product,
                                                Chemistry::ClassExtensions::ToMz(closestExperimentalMass->monoisotopicMass,
                                                                                 closestExperimentalMass->charge),
                                                closestExperimentalMass->totalIntensity,
                                                closestExperimentalMass->charge);
                    matchedFragmentIons.push_back(&tempVar2);
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
        StringEventArgs tempVar(v, nestedIds);
        //WarnHandler +++ nullptr ? nullptrs : WarnHandler->Invoke(this, &tempVar);
        if ( WarnHandler != nullptr) {
            WarnHandler->Invoke(tempVar);
        }
    }

    void MetaMorpheusEngine::Status(const std::string &v)
    {
        StringEventArgs tempVar(v, nestedIds);
        if ( OutLabelStatusHandler != nullptr ) {
            OutLabelStatusHandler->Invoke(tempVar);
        }
    }
    
    void MetaMorpheusEngine::ReportProgress(ProgressEventArgs *v)
    {
        if ( OutProgressHandler != nullptr ) {
            OutProgressHandler->Invoke(v);
        }
    }
    
    void MetaMorpheusEngine::StartingSingleEngine()
    {
        SingleEngineEventArgs tempVar(this);
        if ( StartingSingleEngineHandler != nullptr ) {
            StartingSingleEngineHandler->Invoke(tempVar);
        }
    }
    
    void MetaMorpheusEngine::FinishedSingleEngine(MetaMorpheusEngineResults *myResults)
    {
        SingleEngineFinishedEventArgs tempVar(myResults);
        if ( FinishedSingleEngineHandler != nullptr ) {
            FinishedSingleEngineHandler->Invoke(tempVar);
        }
    }
}
