#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Test
{
    class SearchEngineTests final
    {
	public:
		static void TestClassicSearchEngine();

		static void TestClassicSearchEngineWithWeirdPeptide();

		static void TestModernSearchEngine();

		static void TestModernSearchFragmentEdges();

		static void TestNonViablePSM();

		static void TestModernSearchEngineWithWeirdPeptide();

		static void TestNonSpecificEnzymeSearchEngineSingleN();

		static void TestNonSpecificEnzymeSearchEngineSingleCModifications();

		static void TestNonSpecificEnzymeSearchEngineSingleNModifications();

		static void TestNonSpecificEnzymeSearchEngineSingleC();

		static void TestNonSpecificEnzymeVariableModificationHandlingNTerm();

		static void TestNonSpecificEnzymeVariableModificationHandlingCTerm();

		static void TestSemiSpecificEnzymeEngineSingleN();

		static void TestSemiSpecificEnzymeEngineSingleC();

		static void TestClassicSemiProtease();

		static void TestClassicSemiProteolysis();

		static void TestClassicSearchOneNterminalModifiedPeptideOneScan();
	};
}
