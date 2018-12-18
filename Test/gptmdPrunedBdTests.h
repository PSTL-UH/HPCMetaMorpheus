#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <tuple>
#include "tangible_filesystem.h"

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class GptmdPrunedDbTests
	class GptmdPrunedDbTests final
	{
		// want a psm whose base sequence is not ambigous but full sequence is (ptm is not localized): make sure this does not make it in DB

	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestPrunedGeneration()
		static void TestPrunedGeneration();

		//test if prunedDatabase matches expected output
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestPrunedDatabase()
		static void TestPrunedDatabase();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestUserModSelectionInPrunedDB()
		static void TestUserModSelectionInPrunedDB();
	};
}
