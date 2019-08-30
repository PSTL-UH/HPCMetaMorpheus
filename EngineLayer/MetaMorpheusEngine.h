#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>
#include "tangible_event.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
//namespace EngineLayer { class CommonParameters; }
#include "CommonParameters.h"

namespace EngineLayer { class Ms2ScanWithSpecificMass; }
namespace EngineLayer { class MetaMorpheusEngineResults; }

//namespace EngineLayer { class ProgressEventArgs; }
#include "EventArgs/ProgressEventArgs.h"
#include "EventArgs/SingleEngineEventArgs.h"
#include "EventArgs/SingleEngineFinishedEventArgs.h"
#include "EventArgs/StringEventArgs.h"


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
		static const std::unordered_map<DissociationType, double> complementaryIonConversionDictionary;

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
		static TangibleEvent<EventHandler<SingleEngineEventArgs>> *StartingSingleEngineHander = new TangibleEvent<EventHandler<SingleEngineEventArgs>>();

		static TangibleEvent<EventHandler<SingleEngineFinishedEventArgs>> *FinishedSingleEngineHandler = new TangibleEvent<EventHandler<SingleEngineFinishedEventArgs>>();

		static TangibleEvent<EventHandler<StringEventArgs>> *OutLabelStatusHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler<StringEventArgs>> *WarnHandler = new TangibleEvent<EventHandler<StringEventArgs>>();

		static TangibleEvent<EventHandler<ProgressEventArgs>> *OutProgressHandler = new TangibleEvent<EventHandler<ProgressEventArgs>>();

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
