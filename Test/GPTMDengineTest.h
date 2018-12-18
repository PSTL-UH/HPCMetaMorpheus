#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <tuple>

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::Gptmd;
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
//ORIGINAL LINE: [TestFixture] public static class GptmdEngineTest
	class GptmdEngineTest final
	{
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test][TestCase("NNNNN", "accession", "TangibleTempVerbatimOpenTagnot appliedTangibleTempVerbatimCloseTag", 5)][TestCase("NNNNN", "accession", "TangibleTempVerbatimOpenTag1TangibleTempVerbatimBackslasht50000000TangibleTempVerbatimBackslasht.TangibleTempVerbatimBackslashtATangibleTempVerbatimBackslashtGTangibleTempVerbatimBackslasht.TangibleTempVerbatimBackslashtPASSTangibleTempVerbatimBackslashtANN=G||||||||||||||||TangibleTempVerbatimBackslashtGT:AD:DPTangibleTempVerbatimBackslasht1/1:30,30:30TangibleTempVerbatimCloseTag", 4)] public static void TestGptmdEngine(string proteinSequence, string accession, string sequenceVariantDescription, int numModifiedResidues)
		static void TestGptmdEngine(const std::wstring &proteinSequence, const std::wstring &accession, const std::wstring &sequenceVariantDescription, int numModifiedResidues);

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test][TestCase("NNNPPP", "accession", "A", "TangibleTempVerbatimOpenTagnot appliedTangibleTempVerbatimCloseTag", 1, 6, 3, 3, 0)][TestCase("NNNPPP", "accession", "A", "TangibleTempVerbatimOpenTag1TangibleTempVerbatimBackslasht50000000TangibleTempVerbatimBackslasht.TangibleTempVerbatimBackslashtATangibleTempVerbatimBackslashtGTangibleTempVerbatimBackslasht.TangibleTempVerbatimBackslashtPASSTangibleTempVerbatimBackslashtANN=G||||||||||||||||TangibleTempVerbatimBackslashtGT:AD:DPTangibleTempVerbatimBackslasht1/1:30,30:30TangibleTempVerbatimCloseTag", 1, 5, 2, 3, 0)][TestCase("NNNPPP", "accession", "P", "TangibleTempVerbatimOpenTag1TangibleTempVerbatimBackslasht50000000TangibleTempVerbatimBackslasht.TangibleTempVerbatimBackslashtATangibleTempVerbatimBackslashtGTangibleTempVerbatimBackslasht.TangibleTempVerbatimBackslashtPASSTangibleTempVerbatimBackslashtANN=G||||||||||||||||TangibleTempVerbatimBackslashtGT:AD:DPTangibleTempVerbatimBackslasht1/1:30,30:30TangibleTempVerbatimCloseTag", 2, 5, 2, 3, 1)] public static void TestCombos(string proteinSequence, string accession, string variantAA, string sequenceVariantDescription, int numModHashes, int numModifiedResidues, int numModifiedResiduesN, int numModifiedResiduesP, int numModifiedResiduesNP)
		static void TestCombos(const std::wstring &proteinSequence, const std::wstring &accession, const std::wstring &variantAA, const std::wstring &sequenceVariantDescription, int numModHashes, int numModifiedResidues, int numModifiedResiduesN, int numModifiedResiduesP, int numModifiedResiduesNP);

		//[Test]
		//public static void LoadOriginalMismatchedModifications()
		//{
		//    var protein = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "oblm.xml"), true, DecoyType.Reverse, null, false, null, out var unknownModifications);
		//    Assert.AreEqual(0, protein[0].OneBasedPossibleLocalizedModifications.Count);
		//    protein[0].RestoreUnfilteredModifications();
		//    Assert.AreEqual(1, protein[0].OneBasedPossibleLocalizedModifications.Count);
		//}

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestSearchPtmVariantDatabase()
		static void TestSearchPtmVariantDatabase();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test][TestCase("P", "PETID", "junk", 1, 5, 1, false)][TestCase("P", "PETID", "Unassigned.", 1, 5, 1, false)][TestCase("P", "PETID", "Anywhere.", 1, 5, 1, true)][TestCase("P", "PETID", "N-terminal.", 1, 5, 1, true)][TestCase("P", "PETID", "Peptide N-terminal.", 1, 5, 1, true)][TestCase("P", "PETID", "C-terminal.", 1, 5, 1, false)][TestCase("P", "PETID", "Peptide C-terminal.", 1, 5, 1, false)][TestCase("E", "PETID", "Anywhere.", 2, 5, 2, true)][TestCase("E", "PETID", "N-terminal.", 2, 5, 2, true)][TestCase("E", "PETID", "Peptide N-terminal.", 2, 5, 2, false)][TestCase("E", "PETID", "C-terminal.", 2, 5, 2, false)][TestCase("E", "PETID", "Peptide C-terminal.", 2, 5, 2, false)][TestCase("D", "PETID", "Anywhere.", 5, 5, 5, true)][TestCase("D", "PETID", "N-terminal.", 5, 5, 5, false)][TestCase("D", "PETID", "Peptide N-terminal.", 5, 5, 5, false)][TestCase("D", "PETID", "C-terminal.", 5, 5, 5, true)][TestCase("D", "PETID", "Peptide C-terminal.", 5, 5, 5, true)] public static void Test_GptmdEngineModFits(string targetAminoAcid, string proteinSequence, string locationRestriction, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex, bool result)
		static void Test_GptmdEngineModFits(const std::wstring &targetAminoAcid, const std::wstring &proteinSequence, const std::wstring &locationRestriction, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex, bool result);
	};
}
