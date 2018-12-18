#include "ModificationAnalysisTest.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"

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

	void ModificationAnalysisTest::TestModificationAnalysis()
	{
		IScan *scan = new ThisTestScan();

		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"N", motif1);
		Modification *mod1 = new Modification(_originalId: L"mod1", _modificationType: L"myModType", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);

		ModificationMotif motif2;
		ModificationMotif::TryGetMotif(L"L", motif2);
		Modification *mod2 = new Modification(_originalId: L"mod2", _modificationType: L"myModType", _target: motif2, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);

		std::unordered_map<int, std::vector<Modification*>> oneBasedModifications =
		{
			{
				2, {mod1}
			},
			{
				5, {mod2}
			},
			{
				7, {mod1}
			}
		};
		Protein *protein1 = new Protein(L"MNLDLDNDL", L"prot1", oneBasedModifications: oneBasedModifications);

		std::unordered_map<int, Modification*> allModsOneIsNterminus1 =
		{
			{2, mod1}
		};
		DigestionParams tempVar();
		PeptideWithSetModifications *pwsm1 = new PeptideWithSetModifications(protein1, &tempVar, 2, 9, CleavageSpecificity::Unknown, nullptr, 0, allModsOneIsNterminus1, 0);

		std::unordered_map<int, Modification*> allModsOneIsNterminus2 =
		{
			{2, mod1},
			{7, mod1}
		};
		DigestionParams tempVar2();
		PeptideWithSetModifications *pwsm2 = new PeptideWithSetModifications(protein1, &tempVar2, 2, 9, CleavageSpecificity::Unknown, nullptr, 0, allModsOneIsNterminus2, 0);

		std::unordered_map<int, Modification*> allModsOneIsNterminus3 =
		{
			{7, mod1}
		};
		DigestionParams tempVar3();
		PeptideWithSetModifications *pwsm3 = new PeptideWithSetModifications(protein1, &tempVar3, 2, 9, CleavageSpecificity::Unknown, nullptr, 0, allModsOneIsNterminus3, 0);

		std::unordered_map<int, Modification*> allModsOneIsNterminus4 =
		{
			{8, mod1}
		};
		DigestionParams tempVar4();
		PeptideWithSetModifications *pwsm4 = new PeptideWithSetModifications(protein1, &tempVar4, 1, 9, CleavageSpecificity::Unknown, nullptr, 0, allModsOneIsNterminus4, 0);

		DigestionParams tempVar5(maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: std::numeric_limits<int>::max());
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar5, scoreCutoff: 1);

		auto newPsms = std::vector<PeptideSpectralMatch*>
		{
			new PeptideSpectralMatch(pwsm1, 0, 10, 0, scan, CommonParameters->getDigestionParams(), std::vector<MatchedFragmentIon*>()),
			new PeptideSpectralMatch(pwsm1, 0, 10, 0, scan, CommonParameters->getDigestionParams(), std::vector<MatchedFragmentIon*>()),
			new PeptideSpectralMatch(pwsm2, 0, 10, 0, scan, CommonParameters->getDigestionParams(), std::vector<MatchedFragmentIon*>()),
			new PeptideSpectralMatch(pwsm3, 0, 10, 0, scan, CommonParameters->getDigestionParams(), std::vector<MatchedFragmentIon*>()),
			new PeptideSpectralMatch(pwsm4, 0, 10, 0, scan, CommonParameters->getDigestionParams(), std::vector<MatchedFragmentIon*>())
		};

		for (auto psm : newPsms)
		{
			psm->ResolveAllAmbiguities();
		}

		MassDiffAcceptor *searchMode = new SinglePpmAroundZeroSearchMode(5);
		std::vector<Protein*> proteinList = {protein1};

		FdrAnalysisEngine *fdrAnalysisEngine = new FdrAnalysisEngine(newPsms, searchMode->NumNotches, CommonParameters, std::vector<std::wstring>());
		fdrAnalysisEngine->Run();
		CommonParameters tempVar6();
		ModificationAnalysisEngine *modificationAnalysisEngine = new ModificationAnalysisEngine(newPsms, &tempVar6, std::vector<std::wstring>());
		auto res = static_cast<ModificationAnalysisResults*>(modificationAnalysisEngine->Run());

		Assert::AreEqual(2, res->getCountOfEachModSeenOnProteins().size()());
		Assert::AreEqual(2, res->getCountOfEachModSeenOnProteins()[mod1->IdWithMotif]);
		Assert::AreEqual(1, res->getCountOfEachModSeenOnProteins()[mod2->IdWithMotif]);

		Assert::AreEqual(1, res->getCountOfModsSeenAndLocalized().size()());
		Assert::AreEqual(2, res->getCountOfModsSeenAndLocalized()[mod1->IdWithMotif]);

		Assert::AreEqual(0, res->getCountOfAmbiguousButLocalizedModsSeen().size()());

		Assert::AreEqual(0, res->getCountOfUnlocalizedMods().size()());

		Assert::AreEqual(0, res->getCountOfUnlocalizedFormulas().size()());

		delete modificationAnalysisEngine;
		delete fdrAnalysisEngine;
		delete searchMode;
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pwsm4' statement was not added since pwsm4 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pwsm3' statement was not added since pwsm3 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pwsm2' statement was not added since pwsm2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pwsm1' statement was not added since pwsm1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protein1' statement was not added since protein1 was passed to a method or constructor. Handle memory management manually.
		delete mod2;
		delete mod1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
	}

	void ModificationAnalysisTest::TestModificationAnalysisWithNonLocalizedPtms()
	{
		IScan *scan = new ThisTestScan();

		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"N", motif1);
		Modification *mod1 = new Modification(_originalId: L"mod1", _modificationType: L"mt", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10, _neutralLosses: std::unordered_map<DissociationType*, std::vector<double>>
		{
			{
				MassSpectrometry::DissociationType::AnyActivationType, {10}
			}
		});

		std::unordered_map<int, std::vector<Modification*>> oneBasedModifications =
		{
			{
				2, {mod1}
			},
			{
				7, {mod1}
			}
		};
		Protein *protein1 = new Protein(L"MNLDLDNDL", L"prot1", oneBasedModifications: oneBasedModifications);

		std::unordered_map<int, Modification*> allModsOneIsNterminus1 =
		{
			{2, mod1}
		};
		DigestionParams tempVar();
		PeptideWithSetModifications *pwsm1 = new PeptideWithSetModifications(protein1, &tempVar, 2, 9, CleavageSpecificity::Unknown, nullptr, 0, allModsOneIsNterminus1, 0);

		std::unordered_map<int, Modification*> allModsOneIsNterminus3 =
		{
			{7, mod1}
		};
		DigestionParams tempVar2();
		PeptideWithSetModifications *pwsm2 = new PeptideWithSetModifications(protein1, &tempVar2, 2, 9, CleavageSpecificity::Unknown, nullptr, 0, allModsOneIsNterminus3, 0);

		DigestionParams tempVar3(maxMissedCleavages: 0, minPeptideLength: 1);
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar3, scoreCutoff: 1);

		DigestionParams tempVar4();
		PeptideSpectralMatch *myPsm = new PeptideSpectralMatch(pwsm1, 0, 10, 0, scan, &tempVar4, std::vector<MatchedFragmentIon*>());
		myPsm->AddOrReplace(pwsm2, 10, 0, true, std::vector<MatchedFragmentIon*>());

		myPsm->ResolveAllAmbiguities();

		MassDiffAcceptor *searchMode = new SinglePpmAroundZeroSearchMode(5);
		std::vector<Protein*> proteinList = {protein1};

		FdrAnalysisEngine *fdrAnalysisEngine = new FdrAnalysisEngine({myPsm}, searchMode->NumNotches, CommonParameters, std::vector<std::wstring>());
		fdrAnalysisEngine->Run();
		CommonParameters tempVar5();
		ModificationAnalysisEngine *modificationAnalysisEngine = new ModificationAnalysisEngine({myPsm}, &tempVar5, std::vector<std::wstring>());
		auto res = static_cast<ModificationAnalysisResults*>(modificationAnalysisEngine->Run());

		Assert::AreEqual(1, res->getCountOfEachModSeenOnProteins().size()());
		Assert::AreEqual(2, res->getCountOfEachModSeenOnProteins()[mod1->IdWithMotif]);
		Assert::AreEqual(0, res->getCountOfModsSeenAndLocalized().size()());
		Assert::AreEqual(0, res->getCountOfAmbiguousButLocalizedModsSeen().size());
		Assert::AreEqual(1, res->getCountOfUnlocalizedMods()[mod1->IdWithMotif]); // Saw it, but not sure where!
		Assert::AreEqual(0, res->getCountOfUnlocalizedFormulas().size()());

		delete modificationAnalysisEngine;
		delete fdrAnalysisEngine;
		delete searchMode;
//C# TO C++ CONVERTER TODO TASK: A 'delete myPsm' statement was not added since myPsm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pwsm2' statement was not added since pwsm2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pwsm1' statement was not added since pwsm1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protein1' statement was not added since protein1 was passed to a method or constructor. Handle memory management manually.
		delete mod1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
	}

	std::wstring ThisTestScan::getFullFilePath() const
	{
		return L"";
	}

	int ThisTestScan::getOneBasedScanNumber() const
	{
		return 0;
	}

	std::optional<int> ThisTestScan::getOneBasedPrecursorScanNumber() const
	{
		return std::make_optional(0);
	}

	double ThisTestScan::getRetentionTime() const
	{
		return 0;
	}

	int ThisTestScan::getNumPeaks() const
	{
		return 0;
	}

	double ThisTestScan::getTotalIonCurrent() const
	{
		return 0;
	}

	int ThisTestScan::getPrecursorCharge() const
	{
		return 0;
	}

	double ThisTestScan::getPrecursorMonoisotopicPeakMz() const
	{
		return 0;
	}

	double ThisTestScan::getPrecursorMass() const
	{
		return 0;
	}
}
