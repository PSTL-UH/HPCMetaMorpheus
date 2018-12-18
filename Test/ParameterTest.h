#pragma once

#include <string>
#include <vector>

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class ParameterTest
	class ParameterTest final
	{
		//This test exists because file specific parameters code has a tendency to overwrite common parameters and make life horrible 
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestFileSpecifcParameterOverwrite()
		static void TestFileSpecifcParameterOverwrite();
	};
}
