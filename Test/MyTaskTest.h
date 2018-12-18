#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <tuple>
#include "tangible_filesystem.h"

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Nett;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class MyTaskTest
	class MyTaskTest final
	{
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestEverythingRunner()
		static void TestEverythingRunner();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestMultipleFilesRunner()
		static void TestMultipleFilesRunner();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MakeSureFdrDoesntSkip()
		static void MakeSureFdrDoesntSkip();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MakeSureGptmdTaskMatchesExactMatches()
		static void MakeSureGptmdTaskMatchesExactMatches();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestPeptideCount()
		static void TestPeptideCount();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestFileOutput()
		static void TestFileOutput();

		/// <summary>
		/// This tests for a bug in annotating mods in the search task. The situation is that if you search with a fasta database (no mods annotated),
		/// and then do GPTMD, then search with the GPTMD database, the resulting PSM will have a UniProt mod annotated on it.
		/// Also, if GPTMD has a mod with the same name as a UniProt mod, the annotated PSM will be ambiguous between
		/// the UniProt and the MetaMorpheus modification.
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestUniprotNamingConflicts()
		static void TestUniprotNamingConflicts();

		/// <summary>
		/// Tests that pepXML is written
		/// 
		/// TODO: Assert pepXML properties
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestPepXmlOutput()
		static void TestPepXmlOutput();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestModernAndClassicSearch()
		static void TestModernAndClassicSearch();
	};
}
