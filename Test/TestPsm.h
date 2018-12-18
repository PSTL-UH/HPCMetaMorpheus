#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <tuple>
#include "tangible_filesystem.h"

using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::Localization;
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
//ORIGINAL LINE: [TestFixture] public static class TestPsm
	class TestPsm final
	{
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestPsmHeader()
		static void TestPsmHeader();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestQValueFilter()
		static void TestQValueFilter();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestDecoyContaminantsFilter()
		static void TestDecoyContaminantsFilter();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestPsmMatchingToTargetAndDecoyWithSameSequence()
		static void TestPsmMatchingToTargetAndDecoyWithSameSequence();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestPsmMatchingToTargetAndDecoyWithDifferentSequences()
		static void TestPsmMatchingToTargetAndDecoyWithDifferentSequences();
	};
}
