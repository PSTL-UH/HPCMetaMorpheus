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

#ifdef ORIG
        sb->append("Number of proteins within 1% FDR: " + SortedAndScoredProteinGroups.size()([&] (std::any b) {
                    delete sb;
                    return b::QValue < 0.01;
                }));
#endif
        int count=0;
        for ( auto b: SortedAndScoredProteinGroups ) {
            if ( b->getQValue() < 0.01 ) {
                count++;
            }
        }
        sb->append("Number of proteins within 1% FDR: " + std::to_string(count));
                   
        std::string s =  sb->toString();
        delete sb;
        return s;
	}
}
