#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <cmath>
#include <any>
#include <mutex>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
namespace EngineLayer { namespace Calibration { class LabeledDataPoint; } }

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics::AminoAcidPolymer;

namespace EngineLayer
{
    namespace Calibration
    {
        class DataPointAcquisitionEngine : public MetaMorpheusEngine
        {
        private:
            static constexpr double FineResolutionForIsotopeDistCalculation = 0.1;
            
            std::vector<PeptideSpectralMatch*> GoodIdentifications;
            MsDataFile *const MyMsDataFile;
            Tolerance *const MzToleranceForMs1Search;
            const int MinMS1isotopicPeaksNeededForConfirmedIdentification;
            
        public:
            virtual ~DataPointAcquisitionEngine()
            {
                //delete MyMsDataFile;
                //delete MzToleranceForMs1Search;
            }
            
            DataPointAcquisitionEngine(std::vector<PeptideSpectralMatch*> &goodIdentifications,
                                       MsDataFile *myMsDataFile,
                                       Tolerance *mzToleranceForMs1Search,
                                       int minMS1isotopicPeaksNeededForConfirmedIdentification,
                                       CommonParameters *commonParameters,
                                       std::vector<std::string> nestedIds, int verbosityLevel=0);
            
        protected:
            MetaMorpheusEngineResults *RunSpecific() override;
            
        private:
            //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
            //private (List<LabeledDataPoint>, int, int) SearchMS1Spectra(double[] theoreticalMasses,
            //                                                            double[] theoreticalIntensities,
            //                                                            int ms2spectrumIndex, int direction,
            //                                                            int peptideCharge,
            //                                                            PeptideSpectralMatch identification)
            
            std::tuple<std::vector<LabeledDataPoint>, int, int> SearchMS1Spectra(std::vector<double> theoreticalMasses,
                                                                                 std::vector<double> theoreticalIntensities,
                                                                                 int ms2spectrumIndex, int direction,
                                                                                 int peptideCharge,
                                                                                 PeptideSpectralMatch identification);

            
            static std::vector<LabeledDataPoint*> SearchMS2Spectrum(MsDataScan *ms2DataScan, PeptideSpectralMatch *identification);
		};
	}
}
