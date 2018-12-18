#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include "stringhelper.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class CommonParameters; }

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace MzLibUtil;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [TestFixture] internal class VariantSearchTests
	class VariantSearchTests
	{
	private:
		static CommonParameters *CommonParameters;

	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test][TestCase(0, 0, true, "P4V")][TestCase(1, 1, false, "PT4KT")][TestCase(2, 0, true, "P4PPP")][TestCase(3, 0, true, "PPP4P")][TestCase(4, 0, false, "PKPK4PK")][TestCase(5, 1, true, "PTA4KT")][TestCase(6, 0, false, "KKA4K")][TestCase(7, 1, true, "P4V[type:mod on V]")][TestCase(8, 1, true, "P4PP[type:mod on P]P")][TestCase(0, 0, true, "P6V", DecoyType.Reverse)][TestCase(2, 0, true, "P6PPP", DecoyType.Reverse)][TestCase(3, 0, true, "PPP6P", DecoyType.Reverse)][TestCase(7, 1, true, "P6V[type:mod on V]", DecoyType.Reverse)][TestCase(8, 1, true, "P6PP[type:mod on P]P", DecoyType.Reverse)][TestCase(9, 0, true, "PTIDEPEPTIDE4PPP")][TestCase(9, 0, true, "EDITPEPEDITP2PPP", DecoyType.Reverse)] public static void SearchTests(int proteinIdx, int peptideIdx, bool containsVariant, string variantPsmShort, DecoyType decoyType = DecoyType.None)
		static void SearchTests(int proteinIdx, int peptideIdx, bool containsVariant, const std::wstring &variantPsmShort, DecoyType *decoyType = DecoyType::None);

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Test][TestCase("frameshift.xml")][TestCase("frameshift.xml", DecoyType.Reverse)] public void MoreTests(string filename, DecoyType decoyType = DecoyType.None)
		void MoreTests(const std::wstring &filename, DecoyType *decoyType = DecoyType::None);
	};
}
