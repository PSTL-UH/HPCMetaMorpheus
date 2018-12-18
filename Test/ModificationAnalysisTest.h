#pragma once

#include "../EngineLayer/IScan.h"
#include <string>
#include <unordered_map>
#include <vector>
#include <limits>
#include <optional>

using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::ModificationAnalysis;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class ModificationAnalysisTest
	class ModificationAnalysisTest final
	{
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestModificationAnalysis()
		static void TestModificationAnalysis();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestModificationAnalysisWithNonLocalizedPtms()
		static void TestModificationAnalysisWithNonLocalizedPtms();
	};

	class ThisTestScan : public IScan
	{
	public:
		std::wstring getFullFilePath() const override;

		int getOneBasedScanNumber() const override;

		std::optional<int> getOneBasedPrecursorScanNumber() const override;

		double getRetentionTime() const override;

		int getNumPeaks() const override;

		double getTotalIonCurrent() const override;

		int getPrecursorCharge() const override;

		double getPrecursorMonoisotopicPeakMz() const override;

		double getPrecursorMass() const override;
	};
}
