#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <tuple>
#include "tangible_filesystem.h"

using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class MultiProteaseParsimonyTest
	class MultiProteaseParsimonyTest final
	{
		/// <summary>
		/// three proteins proteases will be trypsina and argC two proteins will produce the peptide by both arge C and trypsin (cleave at R)
		/// the  last one can only make it by trypsin. There is unique peptide also for the trypsin only protein.
		/// </summary>
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseTest()
		static void MultiProteaseTest();

		/// <summary>
		/// In this test, we are simulating a single peptide (pepA) with the same base sequence coming from two different protease digestions.
		/// The proteases cleave at the same spots (mimicking the overlap of peptides between Trypsin and either ArgC or LysC).
		/// So pepA is a shared peptide between protein 1 and protein 2 for both proteases.
		/// With only pepA being observed a protein group containign protien1|protein2 would exist.
		/// If for only one protease we observe a unique peptide (pepB) for protein 2 we would then know for certain that protein2 exists in our sample.
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseSamePeptideSameProteinsDifferentProteases()
		static void MultiProteaseSamePeptideSameProteinsDifferentProteases();

		/// <summary>
		/// In this test, we are ensuring that although two peptides may have the same base sequence (ABC) if they result from only a single protein in the "proteome"
		/// when digested with a  protease they should be considered unique.
		/// Therefore, ABC should be a unique peptide for protein 1 with protease A and for protein 2 with protease B.
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseParsimony_SharedSequenceCanBeUniquePeptide()
		static void MultiProteaseParsimony_SharedSequenceCanBeUniquePeptide();

		/// <summary>
		/// In this test, we want to ensure that protein groups that are actually distinguishable becasue of multiprotease data are not being merged.
		/// Without taking into account the protease peptides would result from, these two proteins (1 and2) would have the same peptide base sequences supporting them.
		/// If that was all that was considered for merging indistinguishable protein groups then we would end up with "1|2", but this would be incorrect.
		/// Becasue of multiple proteases, these two proteins are actually distinguishable ABC can only come from protein 1 with protease A and only from protein2 wiht proteaseB.
		/// The same can be said for EFG. Therefore, we should end up with a protein list contianing both protein 1 and protein 2 supported by 2 unique peptides for each!
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseParsimony_IndistringuishableProteinsNowDistinguishable()
		static void MultiProteaseParsimony_IndistringuishableProteinsNowDistinguishable();

		/// <summary>
		/// In this test, we are showing that peptide ABC although coming from the same protein, same location can be 2 separate unique peptides
		/// because of the digestion of the two proteases. The resultant Protein 1 protien group should contian 2 unique peptides both of which have sequence ABC.
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseParsimony_SameAminoAcidsResultInTwoUniquePeptidesForOneProtein()
		static void MultiProteaseParsimony_SameAminoAcidsResultInTwoUniquePeptidesForOneProtein();

		/// <summary>
		/// In this test, the peptide sequence ABC  results in a unique peptide for protein 1 when the sample is digested with protease test5.
		/// But when the sample is digested with protease test6 the base sequece ABC is a shared peptide between protein 2 and 3.
		/// The protein list should contain protein 1 and the protein group protein2|protein3
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseParsimony_TestingPeptideBaseSequenceCanBeBothSharedAndUnique()
		static void MultiProteaseParsimony_TestingPeptideBaseSequenceCanBeBothSharedAndUnique();

		/// <summary>
		/// In this test, like the previous test, ABC base sewuence can be either unique or shared depending on the protease.
		/// Unlike the previous test, only a psm corresponding to the test3 protease digesting, producing the unique peptide occurs.
		/// Therefore, only protein 1 is listed in the protein list. This test ensures in another manner that unique peptide assignment
		/// is operating correctly.
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseParsimony_BaseSequenceCanBeSharedOrUniqueButOnlyUnqiuePSMSeen()
		static void MultiProteaseParsimony_BaseSequenceCanBeSharedOrUniqueButOnlyUnqiuePSMSeen();

		/// <summary>
		/// In this test, generated PSMs are set to various good or bad FDR levels (above or below 1% FDR) and the filterng
		/// of PSMs is tested prior to parsimony. Only PSMs with FDR at or below 1% should make it through to parsimony.
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestPSMFdrFiltering_Simulated()
		static void TestPSMFdrFiltering_Simulated();

		/// <summary>
		/// This test again check the PSM FDR filtering capabilities prior to parsimony, unlike the previous test,
		/// this test utilizes a sample file with real speectra generating FDR from target decoy analysis instead of simulated values.
		/// This test actually checks the code inside the post search analysis test.
		/// </summary>

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void TestPSMFdrFiltering_RealFile()
		static void TestPSMFdrFiltering_RealFile();

		/// <summary>
		/// In this test, the peptide sequence ABC  results in a unique peptide for protein 1 when the sample is digested with protease alpha.
		/// But when the sample is digested with protease beta the base sequece ABC is a shared peptide between protein 2 and 4.
		/// Peptide EFG is shared between protein 3 and 4. This is a more complex testing set to ensure that Parsing of shared peptides when unique proteins
		/// are present is being handled correctly.
		/// The protein list should contain protein 1 and the protein 4.
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseParsimony_TestingSameBaseSequenceSharedandUniqueMoreComplexSample()
		static void MultiProteaseParsimony_TestingSameBaseSequenceSharedandUniqueMoreComplexSample();

		/// <summary>
		/// This test is to ensure that proteins, even with multiprotease are truly indistinguishable.
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseParsimony_TestingActuallyIndistinguisableProteins()
		static void MultiProteaseParsimony_TestingActuallyIndistinguisableProteins();

		/// <summary>
		/// This test ensures that having the main portion of the greedy algorithm be agnostic of protease is valid, and ensures that when needed the algorithm uses protease to break ties
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseParsimony_TestingGreedyAlgorithm()
		static void MultiProteaseParsimony_TestingGreedyAlgorithm();

		/// <summary>
		/// This test ensures that FDR for each psm is calculated accoriding to its protease
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MultiProteaseParsimony_TestingProteaseSpecificFDRCalculations()
		static void MultiProteaseParsimony_TestingProteaseSpecificFDRCalculations();
	};
}
