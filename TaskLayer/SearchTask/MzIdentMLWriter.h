#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <tuple>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class ProteinGroup; }

using namespace Chemistry;
using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{
	class MzIdentMLWriter final
	{
	public:
		static void WriteMzIdentMl(std::vector<PeptideSpectralMatch*> &psms, std::vector<EngineLayer::ProteinGroup*> &groups, std::vector<Modification*> &variableMods, std::vector<Modification*> &fixedMods, std::vector<Protease*> &proteases, double qValueFilter, Tolerance *productTolerance, Tolerance *parentTolerance, int missedCleavages, const std::string &outputPath);

	private:
		static mzIdentML110::Generated::CVParamType *GetUnimodCvParam(Modification *mod);
	};
}
