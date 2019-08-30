#include "MetaMorpheusEngine.h"
#include "CommonParameters.h"
#include "Ms2ScanWithSpecificMass.h"
#include "MetaMorpheusEngineResults.h"
#include "EventArgs/StringEventArgs.h"
#include "EventArgs/ProgressEventArgs.h"
#include "EventArgs/SingleEngineEventArgs.h"
#include "EventArgs/SingleEngineFinishedEventArgs.h"

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics::Fragmentation;

namespace EngineLayer
{

const std::unordered_map<DissociationType, double> MetaMorpheusEngine::complementaryIonConversionDictionary = std::unordered_map<DissociationType, double>
{
	{DissociationType::HCD, Constants::ProtonMass},
	{DissociationType::ETD, 2 * Constants::ProtonMass},
	{DissociationType::CID, Constants::ProtonMass}
};

	MetaMorpheusEngine::MetaMorpheusEngine(CommonParameters *commonParameters, std::vector<std::string> &nestedIds) : commonParameters(commonParameters), nestedIds(nestedIds)
	{
	}

	double MetaMorpheusEngine::CalculatePeptideScore(MsDataScan *thisScan, std::vector<MatchedFragmentIon*> &matchedFragmentIons, double maximumMassThatFragmentIonScoreIsDoubled)
	{
		double score = 0;

		for (auto fragment : matchedFragmentIons)
		{
			double fragmentScore = 1 + (fragment->Intensity / thisScan->TotalIonCurrent);
			score += fragmentScore;

			if (fragment->NeutralTheoreticalProduct.NeutralMass <= maximumMassThatFragmentIonScoreIsDoubled)
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
		if (!scan->getExperimentalFragments().Any())
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
				MatchedFragmentIon tempVar(product, closestExperimentalMass->monoisotopicMass.ToMz(closestExperimentalMass->charge), closestExperimentalMass->peaks.First().intensity, closestExperimentalMass->charge);
				matchedFragmentIons.push_back(&tempVar);
			}
		}
		if (commonParameters->getAddCompIons())
		{
			double protonMassShift = complementaryIonConversionDictionary[commonParameters->getDissociationType()].ToMass(1);

			for (auto product : theoreticalProducts)
			{
				// unknown fragment mass or diagnostic ion or precursor; skip those
				if (std::isnan(product->NeutralMass) || product->ProductType == ProductType::D || product->ProductType == ProductType::M)
				{
					continue;
				}

				double compIonMass = scan->getPrecursorMass() + protonMassShift - product->NeutralMass;

				// get the closest peak in the spectrum to the theoretical peak
				auto closestExperimentalMass = scan->GetClosestExperimentalFragmentMass(compIonMass);

				// is the mass error acceptable?
				if (commonParameters->getProductMassTolerance()->Within(closestExperimentalMass->monoisotopicMass, compIonMass) && closestExperimentalMass->charge <= scan->getPrecursorCharge())
				{
					MatchedFragmentIon tempVar2(product, closestExperimentalMass->monoisotopicMass.ToMz(closestExperimentalMass->charge), closestExperimentalMass->totalIntensity, closestExperimentalMass->charge);
					matchedFragmentIons.push_back(&tempVar2);
				}
			}
		}

		return matchedFragmentIons;
	}

	MetaMorpheusEngineResults *MetaMorpheusEngine::Run()
	{
		StartingSingleEngine();
		auto stopWatch = new Stopwatch();
		stopWatch->Start();
		auto myResults = RunSpecific();
		stopWatch->Stop();
		myResults->Time = stopWatch->Elapsed;
		FinishedSingleEngine(myResults);

		delete stopWatch;
		return myResults;
	}

	std::string MetaMorpheusEngine::GetId()
	{
		return std::string::Join(",", nestedIds);
	}

	void MetaMorpheusEngine::Warn(const std::string &v)
	{
		StringEventArgs tempVar(v, nestedIds);
		WarnHandler +== nullptr ? nullptr : WarnHandler::Invoke(this, &tempVar);
	}

	void MetaMorpheusEngine::Status(const std::string &v)
	{
		StringEventArgs tempVar(v, nestedIds);
		OutLabelStatusHandler +== nullptr ? nullptr : OutLabelStatusHandler::Invoke(this, &tempVar);
	}

	void MetaMorpheusEngine::ReportProgress(ProgressEventArgs *v)
	{
		OutProgressHandler +== nullptr ? nullptr : OutProgressHandler::Invoke(this, v);
	}

	void MetaMorpheusEngine::StartingSingleEngine()
	{
		SingleEngineEventArgs tempVar(this);
		StartingSingleEngineHander +== nullptr ? nullptr : StartingSingleEngineHander::Invoke(this, &tempVar);
	}

	void MetaMorpheusEngine::FinishedSingleEngine(MetaMorpheusEngineResults *myResults)
	{
		SingleEngineFinishedEventArgs tempVar(myResults);
		FinishedSingleEngineHandler +== nullptr ? nullptr : FinishedSingleEngineHandler::Invoke(this, &tempVar);
	}
}
