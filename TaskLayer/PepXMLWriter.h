#pragma once

#include <string>
#include <vector>
#include <cmath>
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "DbForTask.h"

using namespace EngineLayer;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{
	class PepXMLWriter final
	{
	public:
		static void WritePepXml(std::vector<PeptideSpectralMatch*> &psms,
                                        std::vector<DbForTask*> &database,
                                        std::vector<Modification*> &variableModifications,
                                        std::vector<Modification*> &fixedModifications,
                                        CommonParameters *CommonParameters,
                                        const std::string &outputPath,
                                        double qValueFilter);
	};
}
