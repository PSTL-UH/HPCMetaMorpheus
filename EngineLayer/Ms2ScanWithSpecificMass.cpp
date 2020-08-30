#include <tuple>
#include <algorithm>


#include "Ms2ScanWithSpecificMass.h"
#include "CommonParameters.h"

#include "Search.h"

using namespace Chemistry;
using namespace MassSpectrometry;

#ifdef TIMING_INFO
#include <sys/time.h>

extern double deconvtime;
extern double sorttime;
extern double looptime;

static double timediff (struct timeval t1, struct timeval t2)
{
    double elapsedtime;
    elapsedtime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedtime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

    return elapsedtime/1000;                            //ms to sec
}
#endif

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

            setExperimentalFragments(!neutralExperimentalFragments.empty() ? neutralExperimentalFragments : GetNeutralExperimentalFragments(mzLibScan, commonParam));
            
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
        return privateTheScan->getOneBasedScanNumber();
    }
    
    std::optional<int> Ms2ScanWithSpecificMass::getOneBasedPrecursorScanNumber() const
    {
        return privateTheScan->getOneBasedPrecursorScanNumber();
    }
    
    double Ms2ScanWithSpecificMass::getRetentionTime() const
    {
        return privateTheScan->getRetentionTime();
    }
    
    int Ms2ScanWithSpecificMass::getNumPeaks() const
    {
        return privateTheScan->getMassSpectrum()->getSize();
    }
    
    double Ms2ScanWithSpecificMass::getTotalIonCurrent() const
    {
        return privateTheScan->getTotalIonCurrent();
    }
    
    std::vector<IsotopicEnvelope*> Ms2ScanWithSpecificMass::GetNeutralExperimentalFragments(MsDataScan *scan, CommonParameters *commonParam)
    {
        int minZ = 1;
        int maxZ = 10;
        auto thisScanMassSpec = scan->getMassSpectrum();

#ifdef TIMING_INFO
        struct timeval t1, t1e;
        struct timeval t2, t2e;
        gettimeofday (&t1, NULL);
#endif
        auto neutralExperimentalFragmentMasses = thisScanMassSpec->Deconvolute(thisScanMassSpec->getRange(),
                                                                               minZ, maxZ,
                                                                               commonParam->getDeconvolutionMassTolerance()->getValue(),
                                                                               commonParam->getDeconvolutionIntensityRatio());
#ifdef TIMING_INFO
        gettimeofday (&t1e, NULL);
        deconvtime += timediff(t1, t1e);
#endif
        
        if (commonParam->getAssumeOrphanPeaksAreZ1Fragments())
        {
            std::unordered_set<double> aCMzs;
            for ( IsotopicEnvelope* p: neutralExperimentalFragmentMasses ) {
                for ( auto v: p->peaks ) {
                    aCMzs.insert(Chemistry::ClassExtensions::RoundedDouble(std::get<0>(v)) );
                }
            }
            std::vector<double> alreadyClaimedMzs (aCMzs.begin(), aCMzs.end() );
            std::sort (alreadyClaimedMzs.begin(), alreadyClaimedMzs.end() );

            auto XArray = thisScanMassSpec->getXArray();
            auto YArray = thisScanMassSpec->getYArray();

#ifdef TIMING_INFO
            gettimeofday (&t2, NULL);
#endif
            for (int i = 0; i < (int) XArray.size(); i++)
            {
                double mz = XArray[i];
                double intensity = YArray[i];
                
                if ( !std::binary_search(alreadyClaimedMzs.begin(), alreadyClaimedMzs.end(), Chemistry::ClassExtensions::RoundedDouble(mz)))
                {
                    std::vector<std::tuple<double, double>> dt = {std::make_tuple(mz, intensity)};
                    auto tempVar = new IsotopicEnvelope (dt, Chemistry::ClassExtensions::ToMass(mz, 1), 1, intensity, 0, 0);
                    neutralExperimentalFragmentMasses.push_back(tempVar);
                }
            }
#ifdef TIMING_INFO
            gettimeofday (&t2e, NULL);
            looptime += timediff(t2, t2e);
#endif
        }

#ifdef TIMING_INFO
        gettimeofday (&t2, NULL);
#endif
        std::sort ( neutralExperimentalFragmentMasses.begin(), neutralExperimentalFragmentMasses.end(), [&]
                    (IsotopicEnvelope *e1, IsotopicEnvelope *e2 ) {
                        return e1->monoisotopicMass < e2->monoisotopicMass;
                    });
#ifdef TIMING_INFO
        gettimeofday (&t2e, NULL);
        sorttime += timediff(t2, t2e);
#endif

        return neutralExperimentalFragmentMasses;
    }

    IsotopicEnvelope *Ms2ScanWithSpecificMass::GetClosestExperimentalFragmentMass(double theoreticalNeutralMass)
    {
        if (DeconvolutedMonoisotopicMasses.empty())
        {
            return nullptr;
        }
        return privateExperimentalFragments[GetClosestFragmentMass(theoreticalNeutralMass).value()];
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
