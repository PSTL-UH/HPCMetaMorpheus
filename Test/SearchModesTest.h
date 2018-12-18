#pragma once

#include "../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include <string>
#include <vector>
#include <limits>
#include "exceptionhelper.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class AllowedIntervalWithNotch; }

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace NUnit::Framework;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class SearchModesTest
	class SearchModesTest final
	{
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestSearchModeTest()
		static void TestSearchModeTest();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestDotSearchMode()
		static void TestDotSearchMode();

		// Accept if scanPrecursorMass*peptideMass>=1.
	private:
		class TestSearchMode : public MassDiffAcceptor
		{
		public:
			TestSearchMode(const std::wstring &fileNameAddition);

			int Accepts(double scanPrecursorMass, double peptideMass) override;

			std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) override;

			std::vector<AllowedIntervalWithNotch*> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) override;

			std::wstring ToProseString() override;
		};
	};
}
