#include "Ms2ScanWithSpecificMass.h"
#include "CommonParameters.h"

using namespace Chemistry;
using namespace MassSpectrometry;

namespace EngineLayer
{

	Ms2ScanWithSpecificMass::Ms2ScanWithSpecificMass(MsDataScan *mzLibScan, double precursorMonoisotopicPeakMz, int precursorCharge, const std::string &fullFilePath, CommonParameters *commonParam, std::vector<IsotopicEnvelope*> &neutralExperimentalFragments)
	{
		PrecursorMonoisotopicPeakMz = precursorMonoisotopicPeakMz;
		PrecursorCharge = precursorCharge;
		PrecursorMass = getPrecursorMonoisotopicPeakMz().ToMass(precursorCharge);
		FullFilePath = fullFilePath;

		TheScan = mzLibScan;

		setExperimentalFragments(neutralExperimentalFragments ? neutralExperimentalFragments : GetNeutralExperimentalFragments(mzLibScan, commonParam));

		if (getExperimentalFragments().Any())
		{
			DeconvolutedMonoisotopicMasses = getExperimentalFragments().Select([&] (std::any p)
			{
				p::monoisotopicMass;
			})->ToArray();
		}
	}

	MsDataScan *Ms2ScanWithSpecificMass::getTheScan() const
	{
		return privateTheScan;
	}

	double Ms2ScanWithSpecificMass::getPrecursorMonoisotopicPeakMz() const
	{
		return privatePrecursorMonoisotopicPeakMz;
	}

	double Ms2ScanWithSpecificMass::getPrecursorMass() const
	{
		return privatePrecursorMass;
	}

	int Ms2ScanWithSpecificMass::getPrecursorCharge() const
	{
		return privatePrecursorCharge;
	}

	std::string Ms2ScanWithSpecificMass::getFullFilePath() const
	{
		return privateFullFilePath;
	}

	std::vector<IsotopicEnvelope*> Ms2ScanWithSpecificMass::getExperimentalFragments() const
	{
		return privateExperimentalFragments;
	}

	void Ms2ScanWithSpecificMass::setExperimentalFragments(const std::vector<IsotopicEnvelope*> &value)
	{
		privateExperimentalFragments = value;
	}

	int Ms2ScanWithSpecificMass::getOneBasedScanNumber() const
	{
		return getTheScan()->OneBasedScanNumber;
	}

	std::optional<int> Ms2ScanWithSpecificMass::getOneBasedPrecursorScanNumber() const
	{
		return getTheScan()->OneBasedPrecursorScanNumber;
	}

	double Ms2ScanWithSpecificMass::getRetentionTime() const
	{
		return getTheScan()->RetentionTime;
	}

	int Ms2ScanWithSpecificMass::getNumPeaks() const
	{
		return getTheScan()->MassSpectrum.Size;
	}

	double Ms2ScanWithSpecificMass::getTotalIonCurrent() const
	{
		return getTheScan()->TotalIonCurrent;
	}

	std::vector<IsotopicEnvelope*> Ms2ScanWithSpecificMass::GetNeutralExperimentalFragments(MsDataScan *scan, CommonParameters *commonParam)
	{
		int minZ = 1;
		int maxZ = 10;

		auto neutralExperimentalFragmentMasses = scan->MassSpectrum.Deconvolute(scan->MassSpectrum.Range, minZ, maxZ, commonParam->getDeconvolutionMassTolerance()->Value, commonParam->getDeconvolutionIntensityRatio()).ToList();

		if (commonParam->getAssumeOrphanPeaksAreZ1Fragments())
		{
			std::unordered_set<double> alreadyClaimedMzs = std::unordered_set<double>(neutralExperimentalFragmentMasses.SelectMany([&] (std::any p)
			{
				p::peaks->Select([&] (std::any v)
				{
					ClassExtensions::RoundedDouble(v::mz)->Value;
				});
			}));

			for (int i = 0; i < scan->MassSpectrum.XArray->Length; i++)
			{
				double mz = scan->MassSpectrum.XArray[i];
				double intensity = scan->MassSpectrum.YArray[i];

				if (!std::find(alreadyClaimedMzs.begin(), alreadyClaimedMzs.end(), ClassExtensions::RoundedDouble(mz)->Value) != alreadyClaimedMzs.end()->Value))
				{
					IsotopicEnvelope tempVar({(mz, intensity)}, mz.ToMass(1), 1, intensity, 0, 0);
					neutralExperimentalFragmentMasses.push_back(&tempVar);
				}
			}
		}

		return neutralExperimentalFragmentMasses.OrderBy([&] (std::any p)
		{
			p::monoisotopicMass;
		})->ToArray();
	}

	IsotopicEnvelope *Ms2ScanWithSpecificMass::GetClosestExperimentalFragmentMass(double theoreticalNeutralMass)
	{
		if (DeconvolutedMonoisotopicMasses.empty())
		{
			return nullptr;
		}
		return getExperimentalFragments()[GetClosestFragmentMass(theoreticalNeutralMass).Value];
	}

	std::optional<int> Ms2ScanWithSpecificMass::GetClosestFragmentMass(double mass)
	{
		if (DeconvolutedMonoisotopicMasses.empty())
		{
			return std::nullopt;
		}
		int index = Array::BinarySearch(DeconvolutedMonoisotopicMasses, mass);
		if (index >= 0)
		{
			return std::make_optional(index);
		}
		index = ~index;

		if (index >= DeconvolutedMonoisotopicMasses.size())
		{
			return std::make_optional(index - 1);
		}
		if (index == 0)
		{
			return std::make_optional(index);
		}

		if (mass - DeconvolutedMonoisotopicMasses[index - 1] > DeconvolutedMonoisotopicMasses[index] - mass)
		{
			return std::make_optional(index);
		}
		return std::make_optional(index - 1);
	}
}
