#pragma once

#include "../MetaMorpheusEngineResults.h"
#include <string>
#include <vector>
#include <cmath>
#include <tuple>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { namespace Calibration { class LabeledDataPoint; } }

#include "../MetaMorpheusEngine.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "Chemistry/Chemistry.h"
using namespace Chemistry;
//using namespace MathNet::Numerics::Statistics;

namespace EngineLayer
{
    namespace Calibration
    {
        class DataPointAquisitionResults : public MetaMorpheusEngineResults
        {
        private:
            std::tuple<double, double> privateMs1InfoTh;
            std::tuple<double, double> privateMs2InfoTh;
            std::tuple<double, double> privateMs1InfoPpm;
            std::tuple<double, double> privateMs2InfoPpm;
            int privateNumMs1MassChargeCombinationsConsidered = 0;
            int privateNumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
            int privateNumMs2MassChargeCombinationsConsidered = 0;
            int privateNumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
            std::vector<LabeledDataPoint*> privateMs1List;
            std::vector<LabeledDataPoint*> privateMs2List;
            
        public:
            DataPointAquisitionResults(MetaMorpheusEngine *dataPointAcquisitionEngine, std::vector<PeptideSpectralMatch*> &psms, std::vector<LabeledDataPoint*> &ms1List, std::vector<LabeledDataPoint*> &ms2List, int numMs1MassChargeCombinationsConsidered, int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks, int numMs2MassChargeCombinationsConsidered, int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            
            std::tuple<double, double> getMs1InfoTh() const;
            std::tuple<double, double> getMs2InfoTh() const;
            std::tuple<double, double> getMs1InfoPpm() const;
            std::tuple<double, double> getMs2InfoPpm() const;
                        
            int getNumMs1MassChargeCombinationsConsidered() const;
            int getNumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks() const;
            int getNumMs2MassChargeCombinationsConsidered() const;
            int getNumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks() const;
            
            std::vector<LabeledDataPoint*> getMs1List() const;
            std::vector<LabeledDataPoint*> getMs2List() const;
            
            double PsmPrecursorMedianPpmError;
            double PsmProductMedianPpmError;
            double PsmPrecursorIqrPpmError;
            double PsmProductIqrPpmError;
            std::vector<PeptideSpectralMatch*> Psms;
            
            int getCount() const;
            
            std::string ToString();
        };
    }
}
