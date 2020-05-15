#pragma once

#include "../MetaMorpheusEngineResults.h"
#include <string>
#include <vector>
#include "stringbuilder.h"

#include "../ProteinParsimony/ProteinGroup.h"
#include "ProteinScoringAndFdrEngine.h"


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
