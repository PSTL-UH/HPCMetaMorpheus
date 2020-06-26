#include <tuple>
#include <algorithm>

#include "Ms2ScanWithSpecificMass.h"
#include "CommonParameters.h"

#include "Search.h"

using namespace Chemistry;
using namespace MassSpectrometry;

namespace EngineLayer
{

	Ms2ScanWithSpecificMass::Ms2ScanWithSpecificMass(MsDataScan *mzLibScan,
                                                         double precursorMonoisotopicPeakMz,
                                                         int precursorCharge,
                                                         const std::string &fullFilePath,
                                                         CommonParameters *commonParam,
                                                         std::vector<IsotopicEnvelope*> &neutralExperimentalFragments)
	{
            privatePrecursorMonoisotopicPeakMz = precursorMonoisotopicPeakMz;
            privatePrecursorCharge = precursorCharge;
            privatePrecursorMass = Chemistry::ClassExtensions::ToMass(getPrecursorMonoisotopicPeakMz(), precursorCharge);
            privateFullFilePath = fullFilePath;

            privateTheScan = mzLibScan;

            setExperimentalFragments(neutralExperimentalFragments.empty() ? neutralExperimentalFragments : GetNeutralExperimentalFragments(mzLibScan, commonParam));
            
#ifdef ORIG
            if (getExperimentalFragments().Any()) {
                DeconvolutedMonoisotopicMasses = getExperimentalFragments().Select([&] (std::any p) {
                        p::monoisotopicMass;
                    })->ToArray();
            }
#endif
            for ( auto p: getExperimentalFragments() ) {
                DeconvolutedMonoisotopicMasses.push_back(p->monoisotopicMass );
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
        return getTheScan()->getOneBasedScanNumber();
    }
    
    std::optional<int> Ms2ScanWithSpecificMass::getOneBasedPrecursorScanNumber() const
    {
        return getTheScan()->getOneBasedPrecursorScanNumber();
    }
    
    double Ms2ScanWithSpecificMass::getRetentionTime() const
    {
        return getTheScan()->getRetentionTime();
    }
    
    int Ms2ScanWithSpecificMass::getNumPeaks() const
    {
        return getTheScan()->getMassSpectrum()->getSize();
    }
    
    double Ms2ScanWithSpecificMass::getTotalIonCurrent() const
    {
        return getTheScan()->getTotalIonCurrent();
    }
    
    std::vector<IsotopicEnvelope*> Ms2ScanWithSpecificMass::GetNeutralExperimentalFragments(MsDataScan *scan, CommonParameters *commonParam)
    {
        int minZ = 1;
        int maxZ = 10;
        
        auto neutralExperimentalFragmentMasses = scan->getMassSpectrum()->Deconvolute(scan->getMassSpectrum()->getRange(),
                                                              minZ, maxZ,
                                                              commonParam->getDeconvolutionMassTolerance()->getValue(),
                                                              commonParam->getDeconvolutionIntensityRatio());
        
        if (commonParam->getAssumeOrphanPeaksAreZ1Fragments())
        {
#ifdef ORIG
            std::unordered_set<double> alreadyClaimedMzs = std::unordered_set<double>(neutralExperimentalFragmentMasses.SelectMany([&] (std::any p)
            {
                p::peaks->Select([&] (std::any v)
            {
                Chemistry::ClassExtensions::RoundedDouble(v::mz)->value();
            });
            }));
#endif
            std::unordered_set<double> alreadyClaimedMzs;
            for ( IsotopicEnvelope* p: neutralExperimentalFragmentMasses ) {
                for ( auto v: p->peaks ) {
                    alreadyClaimedMzs.insert(Chemistry::ClassExtensions::RoundedDouble(std::get<0>(v)) );
                }
            }
            
            for (int i = 0; i < (int) scan->getMassSpectrum()->getXArray().size(); i++)
            {
                double mz = scan->getMassSpectrum()->getXArray()[i];
                double intensity = scan->getMassSpectrum()->getYArray()[i];
                
                if ( std::find(alreadyClaimedMzs.begin(), alreadyClaimedMzs.end(), Chemistry::ClassExtensions::RoundedDouble(mz)) == alreadyClaimedMzs.end() )
                {
                    std::vector<std::tuple<double, double>> dt = {std::make_tuple(mz, intensity)};
                    auto tempVar = new IsotopicEnvelope (dt, Chemistry::ClassExtensions::ToMass(mz, 1), 1, intensity, 0, 0);
                    neutralExperimentalFragmentMasses.push_back(tempVar);
                }
            }
    
        }
#ifdef ORIG
        return neutralExperimentalFragmentMasses.OrderBy([&] (std::any p)  {
                p::monoisotopicMass;
            })->ToArray();
#endif
        std::sort ( neutralExperimentalFragmentMasses.begin(), neutralExperimentalFragmentMasses.end(), [&] (IsotopicEnvelope *e1, IsotopicEnvelope *e2 ) {
                return e1->monoisotopicMass < e2->monoisotopicMass;
            });
        return neutralExperimentalFragmentMasses;
    }

    IsotopicEnvelope *Ms2ScanWithSpecificMass::GetClosestExperimentalFragmentMass(double theoreticalNeutralMass)
    {
        if (DeconvolutedMonoisotopicMasses.empty())
        {
            return nullptr;
        }
        return getExperimentalFragments()[GetClosestFragmentMass(theoreticalNeutralMass).value()];
    }

    std::optional<int> Ms2ScanWithSpecificMass::GetClosestFragmentMass(double mass)
    {
        if (DeconvolutedMonoisotopicMasses.empty())
        {
            return std::nullopt;
        }
        //int index = Array::BinarySearch(DeconvolutedMonoisotopicMasses, mass);
        int index = BinarySearch(DeconvolutedMonoisotopicMasses, mass);
        if (index >= 0)
        {
            return std::make_optional(index);
        }
        index = ~index;
        
        if (index >= (int)DeconvolutedMonoisotopicMasses.size())
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
