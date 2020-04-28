#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
//namespace EngineLayer { class CommonParameters; }
#include "CommonParameters.h"

//namespace EngineLayer { class Ms2ScanWithSpecificMass; }
#include "Ms2ScanWithSpecificMass.h"

namespace EngineLayer { class MetaMorpheusEngineResults; }

//namespace EngineLayer { class ProgressEventArgs; }
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

		CommonParameters *const commonParameters;

		const std::vector<std::string> nestedIds;

	public:
		virtual ~MetaMorpheusEngine()
		{
			delete commonParameters;
		}

	protected:
		MetaMorpheusEngine(CommonParameters *commonParameters, std::vector<std::string> &nestedIds);

	public:
		static std::unordered_map<DissociationType, double> complementaryIonConversionDictionary;

		static EventHandler<SingleEngineEventArgs> *StartingSingleEngineHandler;

		static EventHandler<SingleEngineFinishedEventArgs> *FinishedSingleEngineHandler;

		static EventHandler<StringEventArgs> *OutLabelStatusHandler;

		static EventHandler<StringEventArgs> *WarnHandler;

		static EventHandler<ProgressEventArgs> *OutProgressHandler;

		static double CalculatePeptideScore(MsDataScan *thisScan, std::vector<MatchedFragmentIon*> &matchedFragmentIons, double maximumMassThatFragmentIonScoreIsDoubled);

		static std::vector<MatchedFragmentIon*> MatchFragmentIons(Ms2ScanWithSpecificMass *scan, std::vector<Product*> &theoreticalProducts, CommonParameters *commonParameters);

		MetaMorpheusEngineResults *Run();

		std::string GetId();

	protected:
		void Warn(const std::string &v);

		void Status(const std::string &v);

		void ReportProgress(ProgressEventArgs *v);

		virtual MetaMorpheusEngineResults *RunSpecific() = 0;

	private:
		void StartingSingleEngine();

		void FinishedSingleEngine(MetaMorpheusEngineResults *myResults);
	};
}
