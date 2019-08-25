#pragma once

#include "../MetaMorpheusEngineResults.h"
#include <string>
#include <vector>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class ProteinGroup; }
namespace EngineLayer { class ProteinScoringAndFdrEngine; }


namespace EngineLayer
{
	class ProteinScoringAndFdrResults : public MetaMorpheusEngineResults
	{
	public:
		std::vector<ProteinGroup*> SortedAndScoredProteinGroups;

		ProteinScoringAndFdrResults(ProteinScoringAndFdrEngine *proteinAnalysisEngine);

		std::string ToString();
	};
}
