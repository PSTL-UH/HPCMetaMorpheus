#pragma once

#include <string>
#include <vector>
#include "tangible_filesystem.h"

using namespace EngineLayer;
using namespace NUnit::Framework;
using namespace TaskLayer;
using namespace Proteomics::ProteolyticDigestion;
using namespace Proteomics::Fragmentation;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] public static class SearchTaskTest
	class SearchTaskTest final
	{
		/// <summary>
		/// Tests each type of mass difference acceptor type to make sure values are assigned properly
		/// </summary>
	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void MassDiffAceptorTest()
		static void MassDiffAceptorTest();

		/// <summary>
		/// Tests to make sure custom mass difference acceptor inputs are parsed properly
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void ParseSearchModeTest()
		static void ParseSearchModeTest();

		/// <summary>
		/// Ensures that the minimum peptide length is observed (KLEDHPK)
		/// Ensures semispecific search finds peptides that were cleaved correctly during the first digestion (precursor index is made and searched correctly) (KEDEEDKFDAMGNK)
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void SemiSpecificFullAndSmallMatches()
		static void SemiSpecificFullAndSmallMatches();

		/// <summary>
		/// Ensures semispecific search runs and outputs properly
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void SemiSpecificTest()
		static void SemiSpecificTest();

		/// <summary>
		/// Tests that normalization in a search task works properly with an Experimental Design file read in,
		/// and crashes when that file is absent
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void PostSearchNormalizeTest()
		static void PostSearchNormalizeTest();

		/// <summary>
		/// Test that we don't get a crash if protein groups are not constructed
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void ProteinGroupsNoParsimonyTest()
		static void ProteinGroupsNoParsimonyTest();

		/// <summary>
		/// Test ensures pruned databases are written when contaminant DB is searched
		/// </summary>
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test] public static void PrunedDbWithContaminantsTest()
		static void PrunedDbWithContaminantsTest();
	};
}
