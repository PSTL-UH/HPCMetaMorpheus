﻿#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>

#include "CommonParameters.h"
#include "Ms2ScanWithSpecificMass.h"

namespace EngineLayer { class MetaMorpheusEngineResults; }

#include "EventArgs/ProgressEventArgs.h"
#include "EventArgs/SingleEngineEventArgs.h"
#include "EventArgs/SingleEngineFinishedEventArgs.h"
#include "EventArgs/StringEventArgs.h"

#include "EventHandler.h"

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include "Proteomics/Fragmentation/Fragmentation.h"
using namespace Proteomics::Fragmentation;

namespace EngineLayer
{
	class MetaMorpheusEngine
	{
	protected:

		CommonParameters * commonParameters;

		std::vector<std::string> nestedIds;

	public:
		virtual ~MetaMorpheusEngine()
		{
			delete commonParameters;
		}

	protected:
		MetaMorpheusEngine(CommonParameters *commonParameters, std::vector<std::string> nestedIds,
                                   int verbosityLevel=0);

	public:
		static std::unordered_map<DissociationType, double> complementaryIonConversionDictionary;

		static double CalculatePeptideScore(MsDataScan *thisScan, std::vector<MatchedFragmentIon*> &matchedFragmentIons,
                                                    double maximumMassThatFragmentIonScoreIsDoubled);

		static std::vector<MatchedFragmentIon*> MatchFragmentIons(Ms2ScanWithSpecificMass *scan,
                                                                          std::vector<Product*> &theoreticalProducts,
                                                                          CommonParameters *commonParameters);

		MetaMorpheusEngineResults *Run();

		std::string GetId();

	protected:
		void Warn(const std::string &v);

		void Status(const std::string &v);

		void ReportProgress(ProgressEventArgs *v);
                
		void ReportEngineProgress(std::string key, int value);

		virtual MetaMorpheusEngineResults *RunSpecific() = 0;

	private:
                int privateVerbosityLevel;
		void StartingSingleEngine();
		void FinishedSingleEngine(MetaMorpheusEngineResults *myResults);
	};
}
