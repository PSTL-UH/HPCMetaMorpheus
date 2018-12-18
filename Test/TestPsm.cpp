#include "TestPsm.h"
#include "TestDataFile.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"
#include "../TaskLayer/MetaMorpheusTask.h"

using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::Localization;
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

	void TestPsm::TestPsmHeader()
	{
		DigestionParams *digestionParams = new DigestionParams();
		Protein tempVar(L"MQQQQQQQ", L"accession1", L"org", {new std::tuple<std::wstring, std::wstring>(L"geneNameType", L"geneName")},
		new std::unordered_map<int, std::vector<Modification*>>
		{
			{
				2, {new Modification(L"mod", L"mod")}
			}
		},
		name: L"name", fullName: L"fullName", sequenceVariations: new std::vector<SequenceVariation*> {new SequenceVariation(2, L"P", L"Q", L"changed this sequence")});
		PeptideWithSetModifications *pepWithSetMods = (&tempVar)->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()).First();
		MsDataFile *myMsDataFile = new TestDataFile(pepWithSetMods, L"quadratic");
		MsDataScan *scann = myMsDataFile->GetOneBasedScan(2);
		CommonParameters tempVar2();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(scann, 4, 1, L"", &tempVar2);

		auto theoreticalIons = pepWithSetMods->Fragment(DissociationType::HCD, FragmentationTerminus::Both).ToList();
		CommonParameters tempVar3();
		auto matchedIons = MetaMorpheusEngine::MatchFragmentIons(scan, theoreticalIons, &tempVar3);
		PeptideSpectralMatch *psm = new PeptideSpectralMatch(pepWithSetMods, 1, 2, 3, scan, digestionParams, matchedIons);
		psm->ResolveAllAmbiguities();

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		auto t = psm->ToString();
		auto tabsepheader = PeptideSpectralMatch::GetTabSeparatedHeader();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		Assert::AreEqual(psm->ToString()->Count([&] (std::any f)
		{
		delete psm;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
		delete myMsDataFile;
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
			return f == L'\t';
		}), PeptideSpectralMatch::GetTabSeparatedHeader().Count([&] (std::any f)
		{
		delete psm;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
		delete myMsDataFile;
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
			return f == L'\t';
		}));

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		Assert::AreEqual(psm->ToString()->Count([&] (std::any f)
		{
		delete psm;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
		delete myMsDataFile;
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
			return f == L'\t';
		}), PeptideSpectralMatch::GetTabSeparatedHeader().Count([&] (std::any f)
		{
		delete psm;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
		delete myMsDataFile;
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
			return f == L'\t';
		}));

		Tolerance *fragmentTolerance = new PpmTolerance(10);
		LocalizationEngine tempVar4({psm}, myMsDataFile, new CommonParameters(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, fragmentTolerance), new std::vector<std::wstring>());
		(&tempVar4)->Run();

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		Assert::AreEqual(psm->ToString()->Count([&] (std::any f)
		{
		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete psm' statement was not added since psm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
			return f == L'\t';
		}), PeptideSpectralMatch::GetTabSeparatedHeader().Count([&] (std::any f)
		{
		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete psm' statement was not added since psm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
			return f == L'\t';
		}));

		psm->SetFdrValues(6, 6, 6, 6, 6, 6, 0, 0, 0, true);

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		Assert::AreEqual(psm->ToString()->Count([&] (std::any f)
		{
		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete psm' statement was not added since psm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
			return f == L'\t';
		}), PeptideSpectralMatch::GetTabSeparatedHeader().Count([&] (std::any f)
		{
		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete psm' statement was not added since psm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
			return f == L'\t';
		}));

		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete psm' statement was not added since psm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
	}

	void TestPsm::TestQValueFilter()
	{
		SearchTask *searchTask = new SearchTask();
		CommonParameters tempVar(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, = nullptr, = nullptr, = nullptr, 1);
		searchTask->setCommonParameters(&tempVar);


		SearchTask *searchTask2 = new SearchTask();
		CommonParameters tempVar2(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, = nullptr, = nullptr, = nullptr, 0);
		searchTask2->setCommonParameters(&tempVar2);

		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestPSMOutput)");
		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\PrunedDbSpectra.mzml)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\DbForPrunedDb.fasta)");

		auto engine = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"QValueTest", searchTask)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		engine->Run();

		std::wstring psmFile = FileSystem::combine(outputFolder, LR"(QValueTest\AllPSMs.psmtsv)");
		auto lines = File::ReadAllLines(psmFile);
		Assert::That(lines.size() == 12);

		auto engine2 = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"QValueTest", searchTask2)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		engine2->Run();

		auto lines2 = File::ReadAllLines(psmFile);
		Assert::That(lines2.size() == 7);
		Directory::Delete(outputFolder, true);

		delete engine2;
		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTask2' statement was not added since searchTask2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTask' statement was not added since searchTask was passed to a method or constructor. Handle memory management manually.
	}

	void TestPsm::TestDecoyContaminantsFilter()
	{
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestPSMOutput)");
		std::wstring myFile = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\PrunedDbSpectra.mzml)");
		std::wstring myDatabase = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestData\DbForPrunedDb.fasta)");

		//Filter decoys
		SearchTask *searchTaskDecoy = new SearchTask();

		searchTaskDecoy->getSearchParameters()->setWriteDecoys(false);

		auto engineDecoy = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"DecoyTest", searchTaskDecoy)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		engineDecoy->Run();

		std::wstring psmFileDecoy = FileSystem::combine(outputFolder, LR"(DecoyTest\AllPSMs.psmtsv)");
		auto linesDecoy = File::ReadAllLines(psmFileDecoy);
		Assert::That(linesDecoy.size() == 9);

		//Filter contaminants
		SearchTask *searchTaskContaminant = new SearchTask();

		searchTaskContaminant->getSearchParameters()->setWriteContaminants(false);

		auto engineContaminant = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"ContaminantTest", searchTaskContaminant)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, false)}, outputFolder);
		engineContaminant->Run();

		std::wstring psmFileContaminant = FileSystem::combine(outputFolder, LR"(ContaminantTest\AllPSMs.psmtsv)");
		auto linesContaminant = File::ReadAllLines(psmFileContaminant);
		Assert::That(linesContaminant.size() == 12);

		std::wstring proteinFileContaminant = FileSystem::combine(outputFolder, LR"(ContaminantTest\AllProteinGroups.tsv)");
		auto linesContaminantProtein = File::ReadAllLines(proteinFileContaminant);
		Assert::That(linesContaminantProtein.size() == 7);

		auto engineContaminant2 = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"ContaminantTest", searchTaskContaminant)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, true)}, outputFolder);
		engineContaminant2->Run();

		auto linesContaminant2 = File::ReadAllLines(psmFileContaminant);
		Assert::That(linesContaminant2.size() == 1);

		auto linesContaminantProtein2 = File::ReadAllLines(proteinFileContaminant);
		Assert::That(linesContaminantProtein2.size() == 1);

		//Filter contaminants and decoys
		SearchTask *searchTaskDecoyContaminant = new SearchTask();
		SearchTask *searchTaskDecoyContaminant2 = new SearchTask();

		searchTaskDecoyContaminant2->getSearchParameters()->setWriteContaminants(false);
		searchTaskDecoyContaminant2->getSearchParameters()->setWriteDecoys(false);

		auto engineDecoyContaminant = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"DecoyContaminantTest", searchTaskDecoyContaminant)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, true)}, outputFolder);
		engineContaminant->Run();

		std::wstring psmFileDecoyContaminant = FileSystem::combine(outputFolder, LR"(DecoyContaminantTest\AllPSMs.psmtsv)");
		auto linesDecoyContaminant = File::ReadAllLines(psmFileContaminant);
		Assert::That(linesContaminant.size() == 12);

		auto engineDecoyContaminant2 = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"DecoyContaminantTest", searchTaskDecoyContaminant2)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, true)}, outputFolder);
		engineContaminant2->Run();

		std::wstring psmFileDecoyContaminant2 = FileSystem::combine(outputFolder, LR"(DecoyContaminantTest\AllPSMs.psmtsv)");
		auto linesDecoyContaminant2 = File::ReadAllLines(psmFileContaminant);
		Assert::That(linesContaminant2.size() == 1);

		//No filter
		SearchTask *searchTask = new SearchTask();

		auto engine = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"NoFilterTest", searchTask)}, std::vector<std::wstring> {myFile}, std::vector<DbForTask*> {new DbForTask(myDatabase, true)}, outputFolder);
		engine->Run();

		std::wstring psmFile = FileSystem::combine(outputFolder, LR"(NoFilterTest\AllPSMs.psmtsv)");
		auto lines = File::ReadAllLines(psmFile);
		Assert::That(lines.size() == 12);
		Directory::Delete(outputFolder, true);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTask' statement was not added since searchTask was passed to a method or constructor. Handle memory management manually.
		delete engineDecoyContaminant2;
		delete engineDecoyContaminant;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTaskDecoyContaminant2' statement was not added since searchTaskDecoyContaminant2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTaskDecoyContaminant' statement was not added since searchTaskDecoyContaminant was passed to a method or constructor. Handle memory management manually.
		delete engineContaminant2;
		delete engineContaminant;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTaskContaminant' statement was not added since searchTaskContaminant was passed to a method or constructor. Handle memory management manually.
		delete engineDecoy;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchTaskDecoy' statement was not added since searchTaskDecoy was passed to a method or constructor. Handle memory management manually.
	}

	void TestPsm::TestPsmMatchingToTargetAndDecoyWithSameSequence()
	{
		DigestionParams *digest = new DigestionParams();
		std::vector<Modification*> mods;

		Protein tempVar(L"PEPTIDE", L"");
		PeptideWithSetModifications *target = (&tempVar)->Digest(digest, mods, mods).First();
		Protein tempVar2(L"PEPTIDE", L"", isDecoy: true);
		PeptideWithSetModifications *decoy = (&tempVar2)->Digest(digest, mods, mods).First();

		MsDataFile *msDataFile = new TestDataFile(target);
		MsDataScan *msDataScan = msDataFile->GetOneBasedScan(2);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scanWithMass = new Ms2ScanWithSpecificMass(msDataScan, 4, 1, L"", &tempVar3);

		PeptideSpectralMatch *psm = new PeptideSpectralMatch(target, 0, 1, 1, scanWithMass, digest, nullptr);
		psm->AddOrReplace(decoy, 1, 0, true, nullptr);

		Assert::AreEqual(2, psm->BestMatchingPeptides->Count());
		Assert::That(psm->BestMatchingPeptides.Any([&] (std::any p)
		{
			p::Peptide::Protein::IsDecoy;
		}));

		psm->ResolveAllAmbiguities();

		Assert::AreEqual(1, psm->BestMatchingPeptides->Count());
		Assert::That(psm->BestMatchingPeptides.All([&] (std::any p)
		{
			!p::Peptide::Protein::IsDecoy;
		}));
		Assert::That(!psm->getIsDecoy());

		delete psm;
//C# TO C++ CONVERTER TODO TASK: A 'delete scanWithMass' statement was not added since scanWithMass was passed to a method or constructor. Handle memory management manually.
		delete msDataFile;
//C# TO C++ CONVERTER TODO TASK: A 'delete digest' statement was not added since digest was passed to a method or constructor. Handle memory management manually.
	}

	void TestPsm::TestPsmMatchingToTargetAndDecoyWithDifferentSequences()
	{
		DigestionParams *digest = new DigestionParams();
		std::vector<Modification*> mods;

		Protein tempVar(L"PEPTIDE", L"");
		PeptideWithSetModifications *target = (&tempVar)->Digest(digest, mods, mods).First();
		Protein tempVar2(L"PEPTIDEL", L"", isDecoy: true);
		PeptideWithSetModifications *decoy = (&tempVar2)->Digest(digest, mods, mods).First();

		MsDataFile *msDataFile = new TestDataFile(target);
		MsDataScan *msDataScan = msDataFile->GetOneBasedScan(2);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scanWithMass = new Ms2ScanWithSpecificMass(msDataScan, 4, 1, L"", &tempVar3);

		PeptideSpectralMatch *psm = new PeptideSpectralMatch(target, 0, 1, 1, scanWithMass, digest, nullptr);
		psm->AddOrReplace(decoy, 1, 0, true, nullptr);

		Assert::AreEqual(2, psm->BestMatchingPeptides->Count());
		Assert::That(psm->BestMatchingPeptides.Any([&] (std::any p)
		{
			p::Peptide::Protein::IsDecoy;
		}));

		psm->ResolveAllAmbiguities();

		Assert::AreEqual(2, psm->BestMatchingPeptides->Count());
		Assert::That(psm->getIsDecoy());

		FdrAnalysisEngine tempVar4({psm}, 1, new CommonParameters(), new std::vector<std::wstring>());
		(&tempVar4)->Run();
		Assert::AreEqual(0.5, psm->getFdrInfo()->getCumulativeDecoy());

//C# TO C++ CONVERTER TODO TASK: A 'delete psm' statement was not added since psm was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scanWithMass' statement was not added since scanWithMass was passed to a method or constructor. Handle memory management manually.
		delete msDataFile;
//C# TO C++ CONVERTER TODO TASK: A 'delete digest' statement was not added since digest was passed to a method or constructor. Handle memory management manually.
	}
}
