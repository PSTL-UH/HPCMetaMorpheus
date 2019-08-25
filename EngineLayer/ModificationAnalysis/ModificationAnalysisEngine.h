#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }

using namespace Chemistry;

namespace EngineLayer
{
	namespace ModificationAnalysis
	{
		class ModificationAnalysisEngine : public MetaMorpheusEngine
		{
		private:
			const std::vector<PeptideSpectralMatch*> NewPsms;

		public:
			ModificationAnalysisEngine(std::vector<PeptideSpectralMatch*> &newPsms, CommonParameters *commonParameters, std::vector<std::string> &nestedIds);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;
		};
	}
}
