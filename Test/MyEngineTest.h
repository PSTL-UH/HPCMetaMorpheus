#pragma once

#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/MetaMorpheusEngineResults.h"
#include <string>
#include <vector>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MetaMorpheusEngine; }

using namespace EngineLayer;
using namespace NUnit::Framework;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class MyEngineTest
	class MyEngineTest final
	{

	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestMyEngine()
		static void TestMyEngine();

	private:
		class TestEngine : public MetaMorpheusEngine
		{

		public:
			TestEngine(int level);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;

		private:
			class TestResults : public MetaMorpheusEngineResults
			{

			public:
				TestResults(MetaMorpheusEngine *e);

				std::wstring ToString() override;

			};

		};

	};
}
