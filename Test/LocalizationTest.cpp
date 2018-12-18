#include "LocalizationTest.h"
#include "TestDataFile.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/PeptideSpectralMatch.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::Localization;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace Test
{

	void LocalizationTest::TestNonSpecific()
	{
		Protease *p = ProteaseDictionary::Dictionary[L"non-specific"];
		Protein *prot = new Protein(L"MABCDEFGH", nullptr);

		DigestionParams *digestionParams = new DigestionParams(protease: p->Name, maxMissedCleavages: 8, minPeptideLength: 1, maxPeptideLength: 9, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain);
		auto peps = prot->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()).ToList();

		 Assert::AreEqual(1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9, peps.size());

//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
		delete prot;
	}

	void LocalizationTest::TestLocalization()
	{
		Protein *parentProteinForMatch = new Protein(L"MEK", nullptr);
		DigestionParams *digestionParams = new DigestionParams(minPeptideLength: 1);
		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"E", motif);
		std::vector<Modification*> variableModifications = {new Modification(_originalId: L"21", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 21.981943)};

		std::vector<PeptideWithSetModifications*> allPeptidesWithSetModifications = parentProteinForMatch->Digest(digestionParams, std::vector<Modification*>(), variableModifications).ToList();
		Assert::AreEqual(4, allPeptidesWithSetModifications.size()());
		PeptideWithSetModifications *ps = allPeptidesWithSetModifications.front();

		PeptideWithSetModifications *pepWithSetModsForSpectrum = allPeptidesWithSetModifications[1];
		MsDataFile *myMsDataFile = new TestDataFile(std::vector<PeptideWithSetModifications*> {pepWithSetModsForSpectrum});
		Tolerance *fragmentTolerance = new AbsoluteTolerance(0.01);

		CommonParameters tempVar();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(myMsDataFile->GetAllScansList().Last(), pepWithSetModsForSpectrum->MonoisotopicMass.ToMz(1), 1, L"", &tempVar);

		auto theoreticalProducts = ps->Fragment(DissociationType::HCD, FragmentationTerminus::Both).ToList();
		CommonParameters tempVar2();
		auto matchedIons = MetaMorpheusEngine::MatchFragmentIons(scan, theoreticalProducts, &tempVar2);
		PeptideSpectralMatch *newPsm = new PeptideSpectralMatch(ps, 0, 0, 2, scan, digestionParams, matchedIons);
		newPsm->ResolveAllAmbiguities();

		CommonParameters *commonParameters = new CommonParameters(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, fragmentTolerance);

		LocalizationEngine *f = new LocalizationEngine({newPsm}, myMsDataFile, commonParameters, std::vector<std::wstring>());
		f->Run();

		// single peak matches
		Assert::AreEqual(1, newPsm->MatchedFragmentIons.Where([&] (std::any p)
		{
		delete f;
//C# TO C++ CONVERTER TODO TASK: A 'delete commonParameters' statement was not added since commonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete newPsm' statement was not added since newPsm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
		delete parentProteinForMatch;
			return p::NeutralTheoreticalProduct->ProductType == ProductType::b;
		})->Count()); //including b1 now
		Assert::AreEqual(1, newPsm->MatchedFragmentIons.Where([&] (std::any p)
		{
		delete f;
//C# TO C++ CONVERTER TODO TASK: A 'delete commonParameters' statement was not added since commonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete newPsm' statement was not added since newPsm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
		delete parentProteinForMatch;
			return p::NeutralTheoreticalProduct->ProductType == ProductType::y;
		})->Count());

		// when localizing, three peaks match
		Assert::IsTrue(newPsm->getLocalizedScores()[1] > 4 && newPsm->getLocalizedScores()[1] < 5); //we have another matched ion

		delete f;
//C# TO C++ CONVERTER TODO TASK: A 'delete commonParameters' statement was not added since commonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete newPsm' statement was not added since newPsm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
		delete parentProteinForMatch;
	}

	void LocalizationTest::TestSeparateIonsByTerminus()
	{
		std::vector<ProductType*> allIonTypes = {ProductType::b, ProductType::c, ProductType::zPlusOne, ProductType::y};
		std::vector<std::vector<ProductType*>> separated = ProductTypeMethods::SeparateIonsByTerminus(allIonTypes);
		Assert::IsTrue(separated.size() == 2);
		Assert::IsTrue(separated[0].size() == 2);
		Assert::IsTrue(separated[1].size() == 2);
		Assert::IsTrue(std::find(separated[0].begin(), separated[0].end(), ProductType::b) != separated[0].end());
		Assert::IsTrue(std::find(separated[0].begin(), separated[0].end(), ProductType::c) != separated[0].end());
		Assert::IsTrue(std::find(separated[1].begin(), separated[1].end(), ProductType::y) != separated[1].end());
		Assert::IsTrue(std::find(separated[1].begin(), separated[1].end(), ProductType::zPlusOne) != separated[1].end());
		std::vector<std::vector<ProductType*>> nOnly = ProductTypeMethods::SeparateIonsByTerminus(separated[0]);
		Assert::IsTrue(nOnly.size() == 1);
		Assert::IsTrue(nOnly[0].size() == 2);
		Assert::IsTrue(std::find(nOnly[0].begin(), nOnly[0].end(), ProductType::b) != nOnly[0].end());
		Assert::IsTrue(std::find(nOnly[0].begin(), nOnly[0].end(), ProductType::c) != nOnly[0].end());
		std::vector<std::vector<ProductType*>> cOnly = ProductTypeMethods::SeparateIonsByTerminus(separated[1]);
		Assert::IsTrue(cOnly.size() == 1);
		Assert::IsTrue(cOnly[0].size() == 2);
		Assert::IsTrue(std::find(cOnly[0].begin(), cOnly[0].end(), ProductType::y) != cOnly[0].end());
		Assert::IsTrue(std::find(cOnly[0].begin(), cOnly[0].end(), ProductType::zPlusOne) != cOnly[0].end());
	}
}
