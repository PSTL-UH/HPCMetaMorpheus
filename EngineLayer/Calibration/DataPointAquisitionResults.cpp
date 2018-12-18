#include "DataPointAquisitionResults.h"
#include "LabeledDataPoint.h"
#include "../MetaMorpheusEngine.h"
#include "../PeptideSpectralMatch.h"

using namespace Chemistry;
using namespace MathNet::Numerics::Statistics;
namespace EngineLayer
{
	namespace Calibration
	{

		DataPointAquisitionResults::DataPointAquisitionResults(MetaMorpheusEngine *dataPointAcquisitionEngine, std::vector<PeptideSpectralMatch*> &psms, std::vector<LabeledDataPoint*> &ms1List, std::vector<LabeledDataPoint*> &ms2List, int numMs1MassChargeCombinationsConsidered, int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks, int numMs2MassChargeCombinationsConsidered, int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks) : MetaMorpheusEngineResults(dataPointAcquisitionEngine), PsmPrecursorMedianPpmError(Statistics::Median(precursorErrors)), PsmProductMedianPpmError(Statistics::Median(productErrors)), PsmPrecursorIqrPpmError(Statistics::InterquartileRange(precursorErrors)), PsmProductIqrPpmError(Statistics::InterquartileRange(productErrors)), Psms(psms)
		{

			Ms1List = ms1List;
			Ms2List = ms2List;

			Ms1InfoTh = getMs1List().Select([&] (std::any b)
			{
				return b::ExperimentalMz - b::TheoreticalMz;
			}).MeanStandardDeviation();
			Ms2InfoTh = getMs2List().Select([&] (std::any b)
			{
				return b::ExperimentalMz - b::TheoreticalMz;
			}).MeanStandardDeviation();

			Ms1InfoPpm = getMs1List().Select([&] (std::any b)
			{
				return (b::ExperimentalMz - b::TheoreticalMz) / b::TheoreticalMz;
			}).MeanStandardDeviation();
			Ms2InfoPpm = getMs2List().Select([&] (std::any b)
			{
				return (b::ExperimentalMz - b::TheoreticalMz) / b::TheoreticalMz;
			}).MeanStandardDeviation();

			NumMs1MassChargeCombinationsConsidered = numMs1MassChargeCombinationsConsidered;
			NumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
			NumMs2MassChargeCombinationsConsidered = numMs2MassChargeCombinationsConsidered;
			NumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

			auto precursorErrors = psms.Select([&] (std::any p)
			{
				return (p::ScanPrecursorMass - p::PeptideMonisotopicMass->Value) / p::PeptideMonisotopicMass->Value * 1e6;
			}).ToList();

			auto productErrors = psms.Where([&] (std::any p)
			{
				return p::MatchedFragmentIons != nullptr;
			}).SelectMany([&] (std::any p)
			{
				p::MatchedFragmentIons;
			})->Select([&] (std::any p)
			{
				return (p::Mz::ToMass(p::Charge) - p::NeutralTheoreticalProduct::NeutralMass) / p::NeutralTheoreticalProduct::NeutralMass * 1e6;
			}).ToList();
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

		std::wstring DataPointAquisitionResults::ToString()
		{
			auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->appendLine(MetaMorpheusEngineResults::ToString());
			sb->appendLine(L"MS1 calibration datapoint count: " + std::to_wstring(getMs1List().size()));
			sb->appendLine(L"MS1 ppm error median: " + std::round(PsmPrecursorMedianPpmError * std::pow(10, 3)) / std::pow(10, 3));
			sb->appendLine(L"MS1 ppm error interquartile range: " + std::round(PsmPrecursorIqrPpmError * std::pow(10, 3)) / std::pow(10, 3));

			sb->appendLine(L"MS2 calibration datapoint count: " + std::to_wstring(getMs2List().size()));
			sb->appendLine(L"MS2 ppm error median: " + std::round(PsmProductMedianPpmError * std::pow(10, 3)) / std::pow(10, 3));
			sb->appendLine(L"MS2 ppm error interquartile range: " + std::round(PsmProductIqrPpmError * std::pow(10, 3)) / std::pow(10, 3));

			delete sb;
			return sb->toString();
		}
	}
}
