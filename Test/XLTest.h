#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class CommonParameters; }
namespace TaskLayer { class XlSearchParameters; }

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::CrosslinkSearch;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Nett;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::AminoAcidPolymer;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class XLTest
	class XLTest final
	{
	private:
		static IndexingResults *privateindexResults;
		static CommonParameters *privatecommonParameters;
		static XlSearchParameters *privatexlSearchParameters;
		static std::vector<Protein*> privateproteinList;
		static std::vector<Modification*> privatevariableModifications;
		static std::vector<Modification*> privatefixedModifications;
		static Crosslinker *privatecrosslinker;
		static std::vector<PeptideWithSetModifications*> privatedigestedList;

		static IndexingResults *getindexResults();
		static void setindexResults(IndexingResults *value);
		static CommonParameters *getcommonParameters();
		static void setcommonParameters(CommonParameters *value);
		static XlSearchParameters *getxlSearchParameters();
		static void setxlSearchParameters(XlSearchParameters *value);
		static std::vector<Protein*> getproteinList();
		static void setproteinList(const std::vector<Protein*> &value);
		static std::vector<Modification*> getvariableModifications();
		static void setvariableModifications(const std::vector<Modification*> &value);
		static std::vector<Modification*> getfixedModifications();
		static void setfixedModifications(const std::vector<Modification*> &value);
		static Crosslinker *getcrosslinker();
		static void setcrosslinker(Crosslinker *value);
		static std::vector<PeptideWithSetModifications*> getdigestedList();
		static void setdigestedList(const std::vector<PeptideWithSetModifications*> &value);

	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void XlTestXlPosCal()
		static void XlTestXlPosCal();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void XlTestGenerateIntensityRanks()
		static void XlTestGenerateIntensityRanks();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void XlTest_BSA_DSSO()
		static void XlTest_BSA_DSSO();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void XlTest_BSA_DSS_file()
		static void XlTest_BSA_DSS_file();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void XlTest_GenerateUserDefinedCrosslinker()
		static void XlTest_GenerateUserDefinedCrosslinker();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void XlTest_DiffCrosslinkSites()
		static void XlTest_DiffCrosslinkSites();

		/// <summary>
		/// Verifies that crosslinker is generated properly
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void CrosslinkCreateTest()
		static void CrosslinkCreateTest();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void DeadendPeptideTest()
		static void DeadendPeptideTest();

		/// <summary>
		/// Makes sure helper methods that generate indices function properly
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void XLSearchWithGeneratedIndices()
		static void XLSearchWithGeneratedIndices();

		/// <summary>
		/// Tests that a crosslinker that links at proline ("P") generates the correct indices
		/// of potential crosslink sites in the sequence PEPTIDE. The indices should be positions 1 and 3
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestGetPossibleCrosslinkerSites()
		static void TestGetPossibleCrosslinkerSites();

		/// <summary>
		/// Generate and test fragments for this loop-linked peptide:
		///    _
		///   | |
		/// PEPTIDE
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestTheoreticalFragmentsLoop()
		static void TestTheoreticalFragmentsLoop();

		/// <summary>
		/// Generate and test fragments for loop-linked peptides with a mod placed at several different locations
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestTheoreticalLoopFragmentsWithMod()
		static void TestTheoreticalLoopFragmentsWithMod();

		/// <summary>
		/// Generate and test fragments for this dead-end peptide (quenched with tris):
		///   
		///   |
		/// PEPTIDE[Methyl]
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestDeadendTris()
		static void TestDeadendTris();

		/// <summary>
		/// Generate and test fragments for this crosslinked peptide pair that was crosslinked
		/// with a non-cleavable crosslinker:
		///  PRLTEIN
		///   |
		/// PEPTIDE
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestTheoreticalFragmentsNonCleavableCrosslink()
		static void TestTheoreticalFragmentsNonCleavableCrosslink();

		/// <summary>
		/// Generate and test fragments for this crosslinked peptide pair that was crosslinked
		/// with a cleavable crosslinker:
		///   PROTEIN
		///   |
		/// PEPTIDE
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestTheoreticalFragmentsCleavableCrosslink()
		static void TestTheoreticalFragmentsCleavableCrosslink();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestWriteToPercolator()
		static void TestWriteToPercolator();

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestWriteNonSingleCross()
		static void TestWriteNonSingleCross();
	};

	class XLTestDataFile : public MsDataFile
	{
		//Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
	public:
		XLTestDataFile();

		std::wstring getFilePath() const;

		std::wstring getName() const;

		void ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities);
	};

	class XLTestDataFileDiffSite : public MsDataFile
	{
		//Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
	public:
		XLTestDataFileDiffSite();

		std::wstring getFilePath() const;

		std::wstring getName() const;

		void ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities);
	};
}
