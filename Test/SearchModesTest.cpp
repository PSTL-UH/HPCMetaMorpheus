#include "SearchModesTest.h"
#include "../EngineLayer/PrecursorSearchModes/AllowedIntervalWithNotch.h"
#include "../EngineLayer/PrecursorSearchModes/DotMassDiffAcceptor.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace NUnit::Framework;

namespace Test
{

	void SearchModesTest::TestSearchModeTest()
	{
		MassDiffAcceptor *sm = new TestSearchMode(L"My custom");
		Assert::IsTrue(sm->Accepts(2, 2) >= 0);
		Assert::IsTrue(sm->Accepts(0.5, 4) >= 0);
		Assert::IsFalse(sm->Accepts(0.5, 0.5) >= 0);
		Assert::IsTrue(sm->Accepts(1, 1) >= 0);
		Assert::AreEqual(2, sm->GetAllowedPrecursorMassIntervalsFromTheoreticalMass(0.5).front().AllowedInterval.Minimum);
		Assert::AreEqual(0.5, sm->GetAllowedPrecursorMassIntervalsFromTheoreticalMass(2).front().AllowedInterval.Minimum);

		delete sm;
	}

	void SearchModesTest::TestDotSearchMode()
	{
		AbsoluteTolerance tempVar(0.1);
		auto dsm1 = new DotMassDiffAcceptor(L"test1", std::vector<double> {0, 1}, &tempVar);

		Assert::IsTrue(dsm1->Accepts(1000, 1000) >= 0);
		Assert::IsTrue(dsm1->Accepts(1000, 1000 + 0.1 / 2) >= 0);
		Assert::IsFalse(dsm1->Accepts(1000, 1000 + 0.1 * 2) >= 0);
		Assert::IsTrue(dsm1->Accepts(1000 + 0.1 / 2, 1000) >= 0);
		Assert::IsFalse(dsm1->Accepts(1000 + 0.1 * 2, 1000) >= 0);

		Assert::IsTrue(dsm1->Accepts(1000 + 1, 1000 + 0.1 / 2) >= 0);
		Assert::IsFalse(dsm1->Accepts(1000 + 1, 1000 + 0.1 * 2) >= 0);
		Assert::IsTrue(dsm1->Accepts(1000 + 1 + 0.1 / 2, 1000) >= 0);
		Assert::IsFalse(dsm1->Accepts(1000 + 1 + 0.1 * 2, 1000) >= 0);

		auto theList = dsm1->GetAllowedPrecursorMassIntervalsFromTheoreticalMass(100).ToList();

		Assert::AreEqual(99.9, theList[0].AllowedInterval::Minimum);
		Assert::AreEqual(100.1, theList[0].AllowedInterval::Maximum);
		Assert::AreEqual(100.9, theList[1].AllowedInterval::Minimum);
		Assert::AreEqual(101.1, theList[1].AllowedInterval::Maximum);

		PpmTolerance tempVar2(5);
		auto dsm2 = new DotMassDiffAcceptor(L"test2", std::vector<double> {0, 1}, &tempVar2);

		Assert::IsTrue(dsm2->Accepts(1000, 1000) >= 0);

		Assert::IsTrue(dsm2->Accepts(1000 * (1 + 5.0 / 1e6 / 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND
		Assert::IsTrue(dsm2->Accepts(1000 * (1 - 5.0 / 1e6 / 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND

		Assert::IsFalse(dsm2->Accepts(1000, 1000 * (1 - 5.0 / 1e6 / 1.0000001)) >= 0); // VERY CAREFUL

		Assert::IsFalse(dsm2->Accepts(1000 * (1 + 5.0 / 1e6 * 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND
		Assert::IsFalse(dsm2->Accepts(1000 * (1 - 5.0 / 1e6 * 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND

		Assert::IsTrue(dsm2->Accepts(1000, 1000 * (1 + 5.0 / 1e6 * 1.0000001)) >= 0); // VERY CAREFUL

		auto theList2 = dsm2->GetAllowedPrecursorMassIntervalsFromTheoreticalMass(1000).ToList();

		Assert::IsTrue(theList2[0].AllowedInterval->Contains(1000));

		Assert::IsTrue(1000 * (1 + 5.0 / 1e6 / 1.0000001) < theList2[0].AllowedInterval::Maximum);
		Assert::IsTrue(1000 * (1 - 5.0 / 1e6 / 1.0000001) > theList2[0].AllowedInterval::Minimum);
		Assert::IsTrue(1000 * (1 + 5.0 / 1e6 * 1.0000001) > theList2[0].AllowedInterval::Maximum);
		Assert::IsTrue(1000 * (1 - 5.0 / 1e6 * 1.0000001) < theList2[0].AllowedInterval::Minimum);

		Assert::IsTrue(theList2[1].AllowedInterval->Contains(1001));

		delete dsm2;
		delete dsm1;
	}

	SearchModesTest::TestSearchMode::TestSearchMode(const std::wstring &fileNameAddition) : MassDiffAcceptor(fileNameAddition)
	{
	}

	int SearchModesTest::TestSearchMode::Accepts(double scanPrecursorMass, double peptideMass)
	{
		return scanPrecursorMass * peptideMass >= 1 ? 1 : -1;
	}

	std::vector<AllowedIntervalWithNotch*> SearchModesTest::TestSearchMode::GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
	{
		DoubleRange tempVar(1 / peptideMonoisotopicMass, std::numeric_limits<double>::max());
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 1);
	}

	std::vector<AllowedIntervalWithNotch*> SearchModesTest::TestSearchMode::GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
	{
		DoubleRange tempVar(std::numeric_limits<double>::lowest(), 1 / peptideMonoisotopicMass);
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return new AllowedIntervalWithNotch(&tempVar, 1);
	}

	std::wstring SearchModesTest::TestSearchMode::ToProseString()
	{
		throw NotImplementedException();
	}
}
