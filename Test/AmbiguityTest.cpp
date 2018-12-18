#include "AmbiguityTest.h"
#include "../EngineLayer/CommonParameters.h"
#include "TestDataFile.h"
#include "../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../EngineLayer/PeptideSpectralMatch.h"

using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;

namespace Test
{

	void AmbiguityTest::TestResolveAmbiguities()
	{
		Protease *protease = new Protease(L"Custom Protease4", CleavageSpecificity::Full, nullptr, nullptr, std::vector<DigestionMotif*> {new DigestionMotif(L"K", nullptr, 1, L"")});
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		DigestionParams tempVar(protease: protease->Name, minPeptideLength: 1);
		CommonParameters *CommonParameters_t = new CommonParameters(, MassSpectrometry::DissociationType::HCD, , , , , true, , , 1, , , , , , , , , , , &tempVar);

		DigestionParams tempVar2(protease: protease->Name, minPeptideLength: 1);
		CommonParameters *CommonParameters_f = new CommonParameters(, MassSpectrometry::DissociationType::HCD, , , , , false, , , 1, , , , , , , , , , , &tempVar2);

		auto myMsDataFile = new TestDataFile();
		auto variableModifications = std::vector<Modification*>();
		auto fixedModifications = std::vector<Modification*>();
		auto proteinList = std::vector<Protein*>
		{
			new Protein(L"MNNKNKNKQQQ", L"Prot1"),
			new Protein(L"MNNNKQQQ", L"Prot2")
		};
		auto searchModes = new SinglePpmAroundZeroSearchMode(5);

		bool DoPrecursorDeconvolution = true;
		bool UseProvidedPrecursorInfo = true;
		double DeconvolutionIntensityRatio = 4;
		int DeconvolutionMaxAssumedChargeState = 10;
		Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);

		CommonParameters tempVar3();
		auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, L"", &tempVar3).OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		std::vector<PeptideSpectralMatch*> allPsmsArray_withAmbiguity(listOfSortedms2Scans.size());

		std::vector<PeptideSpectralMatch*> allPsmsArray_withOutAmbiguity(listOfSortedms2Scans.size());

		ClassicSearchEngine tempVar4(allPsmsArray_withAmbiguity, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters_t, new std::vector<std::wstring>());
		(&tempVar4)->Run(); //report all ambiguity TRUE
		ClassicSearchEngine tempVar5(allPsmsArray_withOutAmbiguity, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters_f, new std::vector<std::wstring>());
		(&tempVar5)->Run(); //report all ambiguity FALSE

		Assert::AreEqual(L"QQQ", allPsmsArray_withAmbiguity[0]->getBaseSequence());
		Assert::AreEqual(L"QQQ", allPsmsArray_withOutAmbiguity[0]->getBaseSequence());
		Assert::IsTrue(allPsmsArray_withAmbiguity[0]->getProteinLength() == nullptr);
		Assert::IsTrue(allPsmsArray_withOutAmbiguity[0]->getProteinLength() != nullptr);
		Assert::IsTrue(allPsmsArray_withAmbiguity[0]->getOneBasedStartResidueInProtein() == nullptr);
		Assert::IsTrue(allPsmsArray_withOutAmbiguity[0]->getOneBasedStartResidueInProtein() != nullptr);

		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters_f' statement was not added since CommonParameters_f was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters_t' statement was not added since CommonParameters_t was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}
}
