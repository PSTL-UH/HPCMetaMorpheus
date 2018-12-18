#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <limits>

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::Indexing;
using namespace EngineLayer::ModernSearch;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class MyPeptideTest
	class MyPeptideTest final
	{
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestIdenticalPeaks()
		static void TestIdenticalPeaks();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestLastPeaks()
		static void TestLastPeaks();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestVeryCloseExperimentalsClassic()
		static void TestVeryCloseExperimentalsClassic();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestVeryCloseExperimentalsModern()
		static void TestVeryCloseExperimentalsModern();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestAllNaN()
		static void TestAllNaN();
	};
}
