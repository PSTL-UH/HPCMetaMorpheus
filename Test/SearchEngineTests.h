#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::Indexing;
using namespace EngineLayer::ModernSearch;
using namespace EngineLayer::NonSpecificEnzymeSearch;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class SearchEngineTests
	class SearchEngineTests final
	{
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestClassicSearchEngine()
		static void TestClassicSearchEngine();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestClassicSearchEngineWithWeirdPeptide()
		static void TestClassicSearchEngineWithWeirdPeptide();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestModernSearchEngine()
		static void TestModernSearchEngine();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestModernSearchFragmentEdges()
		static void TestModernSearchFragmentEdges();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestNonViablePSM()
		static void TestNonViablePSM();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestModernSearchEngineWithWeirdPeptide()
		static void TestModernSearchEngineWithWeirdPeptide();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestNonSpecificEnzymeSearchEngineSingleN()
		static void TestNonSpecificEnzymeSearchEngineSingleN();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestNonSpecificEnzymeSearchEngineSingleCModifications()
		static void TestNonSpecificEnzymeSearchEngineSingleCModifications();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestNonSpecificEnzymeSearchEngineSingleNModifications()
		static void TestNonSpecificEnzymeSearchEngineSingleNModifications();


//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestNonSpecificEnzymeSearchEngineSingleC()
		static void TestNonSpecificEnzymeSearchEngineSingleC();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestNonSpecificEnzymeVariableModificationHandlingNTerm()
		static void TestNonSpecificEnzymeVariableModificationHandlingNTerm();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestNonSpecificEnzymeVariableModificationHandlingCTerm()
		static void TestNonSpecificEnzymeVariableModificationHandlingCTerm();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestSemiSpecificEnzymeEngineSingleN()
		static void TestSemiSpecificEnzymeEngineSingleN();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestSemiSpecificEnzymeEngineSingleC()
		static void TestSemiSpecificEnzymeEngineSingleC();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestClassicSemiProtease()
		static void TestClassicSemiProtease();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestClassicSemiProteolysis()
		static void TestClassicSemiProteolysis();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestClassicSearchOneNterminalModifiedPeptideOneScan()
		static void TestClassicSearchOneNterminalModifiedPeptideOneScan();
	};
}
