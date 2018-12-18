#include "GPTMDengineTest.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../EngineLayer/GlobalVariables.h"
#include "TestDataFile.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"

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

	void GptmdEngineTest::TestGptmdEngine(const std::wstring &proteinSequence, const std::wstring &accession, const std::wstring &sequenceVariantDescription, int numModifiedResidues)
	{
		std::vector<PeptideSpectralMatch*> allResultingIdentifications;
		ModificationMotif motifN;
		ModificationMotif::TryGetMotif(L"N", motifN);
		auto gptmdModifications = std::vector<Modification*> {new Modification(_originalId: L"21", _modificationType: L"mt", _target: motifN, _locationRestriction: L"Anywhere.", _monoisotopicMass: 21.981943)};
		std::vector<std::tuple<double, double>> combos;
		Tolerance *precursorMassTolerance = new PpmTolerance(10);

		allResultingIdentifications = std::vector<PeptideSpectralMatch*>();
		CommonParameters tempVar();
		auto engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, std::unordered_map<std::wstring, Tolerance*>
		{
			{L"filepath", precursorMassTolerance}
		},
		&tempVar, std::vector<std::wstring>());
		auto res = static_cast<GptmdResults*>(engine->Run());
		Assert::AreEqual(0, res->getMods().size());

		auto parentProtein = new Protein(proteinSequence, accession, sequenceVariations: std::vector<SequenceVariation*> {new SequenceVariation(1, L"N", L"A", sequenceVariantDescription)});
		auto variantProteins = parentProtein->GetVariantProteins();

		DigestionParams *digestionParams = new DigestionParams(minPeptideLength: 5);
		std::vector<Modification*> variableModifications;
		auto modPep = variantProteins->SelectMany([&] (std::any p)
		{
			p::Digest(digestionParams, std::vector<Modification*>(), variableModifications);
		}).First();

		//PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
		MsDataScan tempVar2(new MzSpectrum(new double[] {1}, new double[] {1}, false), 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar4();
		Proteomics::AminoAcidPolymer::Peptide tempVar3(modPep->BaseSequence);
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(&tempVar2, ((&tempVar3)->MonoisotopicMass + 21.981943).ToMz(1), 1, L"filepath", &tempVar4);

		auto peptidesWithSetModifications = std::vector<PeptideWithSetModifications*> {modPep};
		PeptideSpectralMatch *newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.front(), 0, 0, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());

		Tolerance *fragmentTolerance = new AbsoluteTolerance(0.01);

		newPsm->SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0, 0, false);
		allResultingIdentifications.push_back(newPsm);

		CommonParameters tempVar5();
		engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, std::unordered_map<std::wstring, Tolerance*>
		{
			{L"filepath", precursorMassTolerance}
		},
		&tempVar5, std::vector<std::wstring>());
		res = static_cast<GptmdResults*>(engine->Run());
		Assert::AreEqual(1, res->getMods().size());
		Assert::AreEqual(numModifiedResidues, res->getMods()[L"accession"].size());

		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete newPsm' statement was not added since newPsm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
		delete parentProtein;
		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete precursorMassTolerance' statement was not added since precursorMassTolerance was passed to a method or constructor. Handle memory management manually.
	}

	void GptmdEngineTest::TestCombos(const std::wstring &proteinSequence, const std::wstring &accession, const std::wstring &variantAA, const std::wstring &sequenceVariantDescription, int numModHashes, int numModifiedResidues, int numModifiedResiduesN, int numModifiedResiduesP, int numModifiedResiduesNP)
	{
		std::vector<PeptideSpectralMatch*> allIdentifications;
		ModificationMotif motifN;
		ModificationMotif::TryGetMotif(L"N", motifN);
		ModificationMotif motifP;
		ModificationMotif::TryGetMotif(L"P", motifP);
		auto gptmdModifications = std::vector<Modification*>
		{
			new Modification(_originalId: L"21", _modificationType: L"mt", _target: motifN, _locationRestriction: L"Anywhere.", _monoisotopicMass: 21.981943),
			new Modification(_originalId: L"16", _modificationType: L"mt", _target: motifP, _locationRestriction: L"Anywhere.", _monoisotopicMass: 15.994915)
		};
		std::vector<std::tuple<double, double>> combos = {std::tuple<double, double>(21.981943, 15.994915)};
		Tolerance *precursorMassTolerance = new PpmTolerance(10);

		auto parentProtein = new Protein(proteinSequence, accession, sequenceVariations: std::vector<SequenceVariation*> {new SequenceVariation(1, L"N", variantAA, sequenceVariantDescription)});
		auto variantProteins = parentProtein->GetVariantProteins();

		DigestionParams *digestionParams = new DigestionParams(minPeptideLength: 5);
		std::vector<Modification*> variableModifications;
		auto modPep = variantProteins->SelectMany([&] (std::any p)
		{
			p::Digest(digestionParams, std::vector<Modification*>(), variableModifications);
		}).First();

		MzSpectrum tempVar(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfd = new MsDataScan(&tempVar, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		Proteomics::AminoAcidPolymer::Peptide tempVar2(modPep->BaseSequence);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfd, ((&tempVar2)->MonoisotopicMass + 21.981943 + 15.994915).ToMz(1), 1, L"filepath", &tempVar3);

		auto peptidesWithSetModifications = std::vector<PeptideWithSetModifications*> {modPep};
		PeptideSpectralMatch *match = new PeptideSpectralMatch(peptidesWithSetModifications.front(), 0, 0, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.front(), 0, 0, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());

		Tolerance *fragmentTolerance = new AbsoluteTolerance(0.01);

		match->SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0, 0, false);
		allIdentifications = {match};

		CommonParameters tempVar4();
		auto engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, std::unordered_map<std::wstring, Tolerance*>
		{
			{L"filepath", precursorMassTolerance}
		},
		&tempVar4, std::vector<std::wstring>());
		auto res = static_cast<GptmdResults*>(engine->Run());
		Assert::AreEqual(numModHashes, res->getMods().size());
		Assert::AreEqual(numModifiedResidues, res->getMods()[L"accession"].size());
		Assert::AreEqual(numModifiedResiduesN, res->getMods()[L"accession"].Where([&] (std::any b)
		{
			b::Item2->OriginalId->Equals(L"21");
		})->Count());
		Assert::AreEqual(numModifiedResiduesP, res->getMods()[L"accession"].Where([&] (std::any b)
		{
			b::Item2->OriginalId->Equals(L"16");
		})->Count());
		TValue hash;
		std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>::const_iterator res->Mods_iterator = res->getMods().find(L"accession_N1P");
		hash = res->Mods_iterator->second;
		Assert::AreEqual(numModifiedResiduesNP, ((hash != nullptr) ? hash : std::unordered_set<std::tuple<int, Modification*>>())->Count);

		delete engine;
		delete fragmentTolerance;
		delete newPsm;
		delete match;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfd' statement was not added since dfd was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
		delete parentProtein;
//C# TO C++ CONVERTER TODO TASK: A 'delete precursorMassTolerance' statement was not added since precursorMassTolerance was passed to a method or constructor. Handle memory management manually.
	}

	void GptmdEngineTest::TestSearchPtmVariantDatabase()
	{
		//Create Search Task
		SearchTask *task1 = new SearchTask();
		SearchParameters tempVar();
		task1->setSearchParameters(&tempVar);
		task1->getSearchParameters()->setSearchTarget(true);
		task1->getSearchParameters()->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
		CommonParameters tempVar2(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, new DigestionParams(minPeptideLength: 5));
		task1->setCommonParameters(&tempVar2);

		//add task to task list
		auto taskList = std::vector<(std::wstring, MetaMorpheusTask)*> {(L"task1", task1)};

		//create modification lists
		std::vector<Modification*> variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			task1->getCommonParameters()->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();

		//protein Creation (One with mod and one without)
		ModificationMotif motifP;
		ModificationMotif::TryGetMotif(L"P", motifP);
		ModificationMotif motifK;
		ModificationMotif::TryGetMotif(L"K", motifK);
		auto variant = new SequenceVariation(3, L"P", L"K", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G|||||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)");
		Protein *testProteinWithMod = new Protein(L"PEPTID", L"accession1", sequenceVariations: std::vector<SequenceVariation*> {variant});
		std::wstring variantAcc = VariantApplication::GetAccession(testProteinWithMod, std::vector<SequenceVariation*> {variant});
		//First Write XML Database
		std::wstring xmlName = L"oblm.xml";

		//Add Mod to list and write XML input database
		auto modList = std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>();
		Modification tempVar3(_originalId: L"acetyl on P", _modificationType: L"type", _target: motifP, _monoisotopicMass: 42, _locationRestriction: L"Anywhere.");
		auto hash = std::unordered_set<std::tuple<int, Modification*>> {std::tuple<int, Modification*>(1, &tempVar3)};
		Modification tempVar4(_originalId: L"acetyl on K", _modificationType: L"type", _target: motifK, _monoisotopicMass: 42, _locationRestriction: L"Anywhere.");
		auto hashVar = std::unordered_set<std::tuple<int, Modification*>> {std::tuple<int, Modification*>(3, &tempVar4)};
		modList.emplace(testProteinWithMod->Accession, hash);
		modList.emplace(variantAcc, hashVar);
		ProteinDbWriter::WriteXmlDatabase(modList, {testProteinWithMod}, xmlName);

		//now write MZML file
		std::any unknownModifications;
		auto variantProteins = ProteinDbLoader::LoadProteinXML(xmlName, true, DecoyType::Reverse, nullptr, false, nullptr, unknownModifications);
		auto variantProtein = variantProteins[0];
		auto variantDecoy = variantProteins[1];
		Assert::AreEqual(0, unknownModifications->Count);

		Assert::AreEqual(2, variantProteins->Count); // target & decoy
		Assert::AreEqual(2, variantProteins[0].OneBasedPossibleLocalizedModifications->Count);
		std::vector<int> foundResidueIndicies = variantProtein->OneBasedPossibleLocalizedModifications->Select([&] (std::any k)
		{
			k::Key;
		}).ToList();
		std::vector<int> expectedResidueIndices = {1, 3};
		Assert::That(foundResidueIndicies, Is::EquivalentTo(expectedResidueIndices));
		Assert::AreEqual(2, variantDecoy->OneBasedPossibleLocalizedModifications->Count);
		foundResidueIndicies = variantDecoy->OneBasedPossibleLocalizedModifications->Select([&] (std::any k)
		{
			k::Key;
		}).ToList();
		expectedResidueIndices = {4, 6}; //originally modified residues are now at the end in the decoy
		Assert::That(foundResidueIndicies, Is::EquivalentTo(expectedResidueIndices));

		auto thisOk = unknownModifications; //for debugging
		auto commonParamsAtThisPoint = task1->getCommonParameters()->getDigestionParams(); //for debugging

		auto digestedList = variantProteins[0].GetVariantProteins()[0].Digest(task1->getCommonParameters()->getDigestionParams(), std::vector<Modification*>(), variableModifications).ToList();
		Assert::AreEqual(4, digestedList.size());

		//Set Peptide with 1 mod at position 3
		PeptideWithSetModifications *pepWithSetMods1 = digestedList[1];

		//Finally Write MZML file
		Assert::AreEqual(L"PEK[type:acetyl on K]TID", pepWithSetMods1->FullSequence); //this might be base sequence
		MsDataFile *myMsDataFile = new TestDataFile(std::vector<PeptideWithSetModifications*> {pepWithSetMods1});
		std::wstring mzmlName = LR"(hello.mzML)";
		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

		//run!
		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mzmlName}, std::vector<DbForTask*> {new DbForTask(xmlName, false)}, Environment::CurrentDirectory);
		engine->Run();

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete testProteinWithMod' statement was not added since testProteinWithMod was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete variant' statement was not added since variant was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task1' statement was not added since task1 was passed to a method or constructor. Handle memory management manually.
	}

	void GptmdEngineTest::Test_GptmdEngineModFits(const std::wstring &targetAminoAcid, const std::wstring &proteinSequence, const std::wstring &locationRestriction, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex, bool result)
	{
		ModificationMotif motif;
		ModificationMotif::TryGetMotif(targetAminoAcid, motif);
		Modification *attemptToLocalize = new Modification(nullptr, nullptr, nullptr, nullptr, _target: motif, _locationRestriction: locationRestriction, _chemicalFormula: nullptr, _monoisotopicMass: 1, _databaseReference: nullptr, _taxonomicRange: nullptr, _keywords: nullptr, _neutralLosses: nullptr, _diagnosticIons: nullptr, _fileOrigin: nullptr);
		std::unordered_map<int, std::vector<Modification*>> oneBasedModifications;
		oneBasedModifications.emplace(proteinOneBasedIndex, std::vector<Modification*> {attemptToLocalize});
		Protein *protein = new Protein(proteinSequence, nullptr, nullptr, nullptr, oneBasedModifications, nullptr, nullptr, nullptr, false, false, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, L"");

		Assert::AreEqual(result, GptmdEngine::ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

//C# TO C++ CONVERTER TODO TASK: A 'delete protein' statement was not added since protein was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete attemptToLocalize' statement was not added since attemptToLocalize was passed to a method or constructor. Handle memory management manually.
	}
}
