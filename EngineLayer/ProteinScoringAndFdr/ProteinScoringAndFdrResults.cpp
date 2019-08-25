#include "ProteinScoringAndFdrResults.h"
#include "../ProteinParsimony/ProteinGroup.h"
#include "ProteinScoringAndFdrEngine.h"


namespace EngineLayer
{

	ProteinScoringAndFdrResults::ProteinScoringAndFdrResults(ProteinScoringAndFdrEngine *proteinAnalysisEngine) : MetaMorpheusEngineResults(proteinAnalysisEngine)
	{
	}

        std::string ProteinScoringAndFdrResults::ToString()
	{
		auto sb = new StringBuilder();

		sb->appendLine(MetaMorpheusEngineResults::ToString());
		sb->append("Number of proteins within 1% FDR: " + SortedAndScoredProteinGroups.size()([&] (std::any b)
		{
		delete sb;
			return b::QValue < 0.01;
		}));

		delete sb;
		return sb->toString();
	}
}
