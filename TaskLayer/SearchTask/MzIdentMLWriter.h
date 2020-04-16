#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <tuple>

#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../../EngineLayer/ProteinParsimony/ProteinGroup.h"

#include "Chemistry/Chemistry.h"
using namespace Chemistry;
using namespace EngineLayer;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

#include "MzIdentML/mzIdentML110.h"

namespace TaskLayer
{
	class MzIdentMLWriter final
	{
	public:
		static void WriteMzIdentMl(std::vector<PeptideSpectralMatch*> &psms,
                                           std::vector<EngineLayer::ProteinGroup*> &groups,
                                           std::vector<Modification*> &variableMods,
                                           std::vector<Modification*> &fixedMods,
                                           std::vector<Protease*> &proteases,
                                           double qValueFilter,
                                           Tolerance *productTolerance,
                                           Tolerance *parentTolerance, int missedCleavages,
                                           const std::string &outputPath);

	private:
		static mzIdentML110::CVParamType *GetUnimodCvParam(Modification *mod);
	};
}
