#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <optional>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }
namespace EngineLayer { namespace FdrAnalysis { class FdrAnalysisResults; } }

using namespace MathNet::Numerics;

namespace EngineLayer
{
	namespace FdrAnalysis
	{
		class FdrAnalysisEngine : public MetaMorpheusEngine
		{
		private:
			std::vector<PeptideSpectralMatch*> AllPsms;
			const int MassDiffAcceptorNumNotches;
			const bool UseDeltaScore;
			const bool CalculateEValue;
			const double ScoreCutoff;

		public:
			FdrAnalysisEngine(std::vector<PeptideSpectralMatch*> &psms, int massDiffAcceptorNumNotches, CommonParameters *commonParameters, std::vector<std::wstring> &nestedIds);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;

		private:
			void DoFalseDiscoveryRateAnalysis(FdrAnalysisResults *myAnalysisResults);

			static double GetEValue(PeptideSpectralMatch *psm, int globalMeanCount, double globalMeanScore, double &maximumLikelihood);

			static int GetNumPSMsAtqValueCutoff(std::vector<PeptideSpectralMatch*> &psms, double qValueCutoff);
		};
	}
}
