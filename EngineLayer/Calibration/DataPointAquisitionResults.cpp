#include "DataPointAquisitionResults.h"
#include "LabeledDataPoint.h"
#include "../MetaMorpheusEngine.h"
#include "../PeptideSpectralMatch.h"
#include "MzLibMath.h"

#include "Chemistry/ClassExtensions.h"
using namespace Chemistry;
//using namespace MathNet::Numerics::Statistics;
namespace EngineLayer
{
	namespace Calibration
	{

		DataPointAquisitionResults::DataPointAquisitionResults(MetaMorpheusEngine *dataPointAcquisitionEngine,
                                                                       std::vector<PeptideSpectralMatch*> &psms,
                                                                       std::vector<LabeledDataPoint*> &ms1List,
                                                                       std::vector<LabeledDataPoint*> &ms2List,
                                                                       int numMs1MassChargeCombinationsConsidered,
                                                                       int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
                                                                       int numMs2MassChargeCombinationsConsidered,
                                                                       int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks) :  MetaMorpheusEngineResults(dataPointAcquisitionEngine)
		{

                    privateMs1List = ms1List;
                    privateMs2List = ms2List;

#ifdef ORIG
                    Ms1InfoTh = getMs1List().Select([&] (std::any b)  {
                            return b::ExperimentalMz - b::TheoreticalMz;
                        }).MeanStandardDeviation();
#endif
                    std::vector<double> dtemp;
                    for ( auto b: privateMs1List ) {
                        dtemp.push_back((b->ExperimentalMz - b->TheoreticalMz));
                    }
                    privateMs1InfoTh = (std::make_tuple(Math::Mean(dtemp), Math::StandardDeviation(dtemp)));

#ifdef ORIG
                    Ms2InfoTh = getMs2List().Select([&] (std::any b)  {
                            return b::ExperimentalMz - b::TheoreticalMz;
			}).MeanStandardDeviation();
#endif
                    dtemp.clear();
                    for ( auto b: privateMs2List ) {
                        dtemp.push_back((b->ExperimentalMz - b->TheoreticalMz));
                    }
                    privateMs2InfoTh = (std::make_tuple(Math::Mean(dtemp), Math::StandardDeviation(dtemp)));
                        
#ifdef ORIG
                    Ms1InfoPpm = getMs1List().Select([&] (std::any b)  	{
                            return (b::ExperimentalMz - b::TheoreticalMz) / b::TheoreticalMz;
			}).MeanStandardDeviation();
#endif
                    dtemp.clear();
                    for ( auto b: privateMs1List ) {
                        dtemp.push_back((b->ExperimentalMz - b->TheoreticalMz)/b->TheoreticalMz);
                    }
                    privateMs1InfoPpm = (std::make_tuple(Math::Mean(dtemp), Math::StandardDeviation(dtemp)));
                    
#ifdef ORIG
                    Ms2InfoPpm = getMs2List().Select([&] (std::any b){
                            return (b::ExperimentalMz - b::TheoreticalMz) / b::TheoreticalMz;
			}).MeanStandardDeviation();
#endif
                    dtemp.clear();
                    for ( auto b: privateMs2List ) {
                        dtemp.push_back((b->ExperimentalMz - b->TheoreticalMz)/b->TheoreticalMz);
                    }
                    privateMs2InfoPpm = (std::make_tuple(Math::Mean(dtemp), Math::StandardDeviation(dtemp)));

                    privateNumMs1MassChargeCombinationsConsidered = numMs1MassChargeCombinationsConsidered;
                    privateNumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
                    privateNumMs2MassChargeCombinationsConsidered = numMs2MassChargeCombinationsConsidered;
                    privateNumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
                    
#ifdef ORIG
                    std::vector<double> precursorErrors = psms.Select([&] (std::any p){
                            return (p::ScanPrecursorMass - p::PeptideMonisotopicMass->Value) / p::PeptideMonisotopicMass->Value * 1e6;
			}).ToList();
#endif
                    std::vector<double> precursorErrors;
                    for ( auto p : psms ) {
                        double d = (p->getScanPrecursorMass() - p->getPeptideMonisotopicMass().value()) / p->getPeptideMonisotopicMass().value() * 1e6;
                        precursorErrors.push_back ( d ) ;
                    }
                    
#ifdef ORIG
                    std::vector<double> productErrors = psms.Where([&] (std::any p){
                            return p::MatchedFragmentIons != nullptr;
			}).SelectMany([&] (std::any p)  {
				p::MatchedFragmentIons;
                            })->Select([&] (std::any p)	{
                                    return (p::Mz::ToMass(p::Charge) - p::NeutralTheoreticalProduct::NeutralMass) / p::NeutralTheoreticalProduct::NeutralMass * 1e6;
                                }).ToList();
#endif
                    std::vector<double> productErrors;
                    for ( auto p : psms ) {
                        if ( !p->getMatchedFragmentIons().empty() ) {
                            for ( auto pp = p->getMatchedFragmentIons().begin(); pp!=p->getMatchedFragmentIons().end(); pp++ ) {
                                productErrors.push_back(((Chemistry::ClassExtensions::ToMass((*pp)->Mz, (*pp)->Charge) - (*pp)->NeutralTheoreticalProduct->NeutralMass) / (*pp)->NeutralTheoreticalProduct->NeutralMass * 1e6));
                            }
                        }
                    }
                    
                    PsmPrecursorMedianPpmError = Math::Median(precursorErrors);
                    PsmProductMedianPpmError   = Math::Median(productErrors);
                    PsmPrecursorIqrPpmError    = Math::InterquartileRange(precursorErrors);
                    PsmProductIqrPpmError      = Math::InterquartileRange(productErrors);
                    Psms = psms;
                    
		}

		std::tuple<double, double> DataPointAquisitionResults::getMs1InfoTh() const
		{
			return privateMs1InfoTh;
		}

		std::tuple<double, double> DataPointAquisitionResults::getMs2InfoTh() const
		{
			return privateMs2InfoTh;
		}

		std::tuple<double, double> DataPointAquisitionResults::getMs1InfoPpm() const
		{
			return privateMs1InfoPpm;
		}

		std::tuple<double, double> DataPointAquisitionResults::getMs2InfoPpm() const
		{
			return privateMs2InfoPpm;
		}

		int DataPointAquisitionResults::getNumMs1MassChargeCombinationsConsidered() const
		{
			return privateNumMs1MassChargeCombinationsConsidered;
		}

		int DataPointAquisitionResults::getNumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks() const
		{
			return privateNumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
		}

		int DataPointAquisitionResults::getNumMs2MassChargeCombinationsConsidered() const
		{
			return privateNumMs2MassChargeCombinationsConsidered;
		}

		int DataPointAquisitionResults::getNumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks() const
		{
			return privateNumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
		}

		std::vector<LabeledDataPoint*> DataPointAquisitionResults::getMs1List() const
		{
			return privateMs1List;
		}

		std::vector<LabeledDataPoint*> DataPointAquisitionResults::getMs2List() const
		{
			return privateMs2List;
		}

		int DataPointAquisitionResults::getCount() const
		{
			return getMs1List().size() + getMs2List().size();
		}

		std::string DataPointAquisitionResults::ToString()
		{
			auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->appendLine(MetaMorpheusEngineResults::ToString());
			sb->appendLine("MS1 calibration datapoint count: " + std::to_string(getMs1List().size()));
			sb->appendLine("MS1 ppm error median: " + std::to_string(std::round(PsmPrecursorMedianPpmError * std::pow(10, 3)) / std::pow(10, 3)));
			sb->appendLine("MS1 ppm error interquartile range: " + std::to_string(std::round(PsmPrecursorIqrPpmError * std::pow(10, 3)) / std::pow(10, 3)));

			sb->appendLine("MS2 calibration datapoint count: " + std::to_string(getMs2List().size()));
			sb->appendLine("MS2 ppm error median: " + std::to_string(std::round(PsmProductMedianPpmError * std::pow(10, 3)) / std::pow(10, 3)));
			sb->appendLine("MS2 ppm error interquartile range: " + std::to_string(std::round(PsmProductIqrPpmError * std::pow(10, 3)) / std::pow(10, 3)));

                        std::string s =  sb->toString();
			delete sb;
                        return s;
		}
	}
}
