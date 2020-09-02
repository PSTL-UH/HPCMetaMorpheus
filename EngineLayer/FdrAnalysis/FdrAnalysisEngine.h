#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <optional>

#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "FdrAnalysisResults.h"


namespace EngineLayer
{
    namespace FdrAnalysis
    {
        class FdrAnalysisEngine : public MetaMorpheusEngine
        {
        private:
            std::vector<PeptideSpectralMatch*>& AllPsms;
            const int MassDiffAcceptorNumNotches;
            const bool UseDeltaScore;
            const bool CalculateEValue;
            const double ScoreCutoff;
            
        public:
            FdrAnalysisEngine(std::vector<PeptideSpectralMatch*> &psms, int massDiffAcceptorNumNotches,
                              CommonParameters *commonParameters, std::vector<std::string> nestedIds,
                              int verbosityLevel=0);
            
        protected:
            MetaMorpheusEngineResults *RunSpecific() override;
            
        private:
            void DoFalseDiscoveryRateAnalysis(FdrAnalysisResults *myAnalysisResults);
            
            static double GetEValue(PeptideSpectralMatch *psm, int globalMeanCount, double globalMeanScore,
                                    double &maximumLikelihood);
            
            static int GetNumPSMsAtqValueCutoff(std::vector<PeptideSpectralMatch*> &psms, double qValueCutoff);
        };
    }
}
