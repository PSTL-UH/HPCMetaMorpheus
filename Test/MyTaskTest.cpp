#include "MyTaskTest.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "TestDataFile.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Nett;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{

	void MyTaskTest::TestEverythingRunner()
	{
		for (auto modFile : Directory::GetFiles(LR"(Mods)"))
		{
			std::any fmww;
			GlobalVariables::AddMods(PtmListLoader::ReadModsFromFile(modFile, fmww), false);
		}

		CalibrationTask *task1 = new CalibrationTask();
		CommonParameters tempVar(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain));
		task1->setCommonParameters(&tempVar);
		CalibrationParameters tempVar2();
		task1->setCalibrationParameters(&tempVar2);
		task1->getCalibrationParameters()->setWriteIntermediateFiles(true);
		task1->getCalibrationParameters()->setNumFragmentsNeededForEveryIdentification(6);

		GptmdTask *task2 = new GptmdTask();
		CommonParameters tempVar3();
		task2->setCommonParameters(&tempVar3);

		SearchTask *task3 = new SearchTask();
		CommonParameters tempVar4();
		task3->setCommonParameters(&tempVar4);
		SearchParameters tempVar5();
		task3->setSearchParameters(&tempVar5);
		task3->getSearchParameters()->setDoParsimony(true);
		task3->getSearchParameters()->setSearchTarget(true);
		task3->getSearchParameters()->setSearchType(SearchType::Modern);

		SearchTask *task4 = new SearchTask();
		CommonParameters tempVar6();
		task4->setCommonParameters(&tempVar6);
		SearchParameters tempVar7();
		task4->setSearchParameters(&tempVar7);
		task4->getSearchParameters()->setSearchType(SearchType::Modern);

		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"task1", task1), (L"task2", task2), (L"task3", task3), (L"task4", task4)};

		std::vector<Modification*> variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			task1->getCommonParameters()->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();
		std::vector<Modification*> fixedModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			task1->getCommonParameters()->ListOfModsFixed->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();

		// Generate data for files
		Protein *ParentProtein = new Protein(L"MPEPTIDEKANTHE", L"accession1");

		auto digestedList = ParentProtein->Digest(task1->getCommonParameters()->getDigestionParams(), fixedModifications, variableModifications).ToList();

		Assert::AreEqual(3, digestedList.size());

		PeptideWithSetModifications *pepWithSetMods1 = digestedList[0];

		PeptideWithSetModifications *pepWithSetMods2 = digestedList[2];

		auto dictHere = std::unordered_map<int, std::vector<Modification*>>();
		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"E", motif);
		dictHere.emplace(3, std::vector<Modification*> {new Modification(_originalId: L"21", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 21.981943)});
		Protein *ParentProteinToNotInclude = new Protein(L"MPEPTIDEK", L"accession2", L"organism", std::vector<std::tuple<std::wstring, std::wstring>>(), dictHere);
		digestedList = ParentProteinToNotInclude->Digest(task1->getCommonParameters()->getDigestionParams(), fixedModifications, variableModifications).ToList();

		MsDataFile *myMsDataFile = new TestDataFile(std::vector<std::vector<PeptideWithSetModifications*>>(1) });

		Protein *proteinWithChain = new Protein(L"MAACNNNCAA", L"accession3", L"organism", std::vector<std::tuple<std::wstring, std::wstring>>(), std::unordered_map<int, std::vector<Modification*>>(), {new ProteolysisProduct(4, 8, L"chain")}, L"name2", L"fullname2");

		std::wstring mzmlName = LR"(ok.mzML)";
		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
		std::wstring xmlName = L"okk.xml";
		ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>(), {ParentProtein, proteinWithChain}, xmlName);

		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestEverythingRunner)");
		// RUN!
		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mzmlName}, std::vector<DbForTask*> {new DbForTask(xmlName, false)}, outputFolder);
		engine->Run();
		File::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, mzmlName));
		File::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, xmlName));
		Directory::Delete(outputFolder, true);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete proteinWithChain' statement was not added since proteinWithChain was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
		delete ParentProteinToNotInclude;
//C# TO C++ CONVERTER TODO TASK: A 'delete ParentProtein' statement was not added since ParentProtein was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task4' statement was not added since task4 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task3' statement was not added since task3 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task2' statement was not added since task2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task1' statement was not added since task1 was passed to a method or constructor. Handle memory management manually.
	}

	void MyTaskTest::TestMultipleFilesRunner()
	{
		for (auto modFile : Directory::GetFiles(LR"(Mods)"))
		{
			std::any fmww;
			GlobalVariables::AddMods(PtmListLoader::ReadModsFromFile(modFile, fmww), false);
		}

		CalibrationTask *task1 = new CalibrationTask();
		CommonParameters tempVar(, , , , = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, new AbsoluteTolerance(0.01), , , , new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain), {(L"Common Variable", L"Oxidation on M")}, new std::vector<(std::wstring, std::wstring)*> {(L"Common Fixed", L"Carbamidomethyl on C")});
		task1->setCommonParameters(&tempVar);
		CalibrationParameters tempVar2();
		task1->setCalibrationParameters(&tempVar2);
		task1->getCalibrationParameters()->setNumFragmentsNeededForEveryIdentification(6);
		GptmdTask *task2 = new GptmdTask();
		CommonParameters tempVar3(, , = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, new AbsoluteTolerance(0.01), , , , new DigestionParams());
		task2->setCommonParameters(&tempVar3);

		SearchTask *task3 = new SearchTask();
		CommonParameters tempVar4();
		task3->setCommonParameters(&tempVar4);
		SearchParameters tempVar5();
		task3->setSearchParameters(&tempVar5);
		task3->getSearchParameters()->setDoParsimony(true);
		task3->getSearchParameters()->setSearchTarget(true);
		task3->getSearchParameters()->setSearchType(SearchType::Modern);
		SearchTask *task4 = new SearchTask();
		CommonParameters tempVar6();
		task4->setCommonParameters(&tempVar6);
		SearchParameters tempVar7();
		task4->setSearchParameters(&tempVar7);
		task4->getSearchParameters()->setSearchType(SearchType::Modern);
		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"task1", task1), (L"task2", task2), (L"task3", task3), (L"task4", task4)};

		std::vector<Modification*> variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			task1->getCommonParameters()->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();
		std::vector<Modification*> fixedModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			task1->getCommonParameters()->ListOfModsFixed->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();

		// Generate data for files
		Protein *ParentProtein = new Protein(L"MPEPTIDEKANTHE", L"accession1");

		auto digestedList = ParentProtein->Digest(task1->getCommonParameters()->getDigestionParams(), fixedModifications, variableModifications).ToList();

		Assert::AreEqual(3, digestedList.size());

		PeptideWithSetModifications *pepWithSetMods1 = digestedList[0];

		PeptideWithSetModifications *pepWithSetMods2 = digestedList[2];

		auto dictHere = std::unordered_map<int, std::vector<Modification*>>();
		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"E", motif);
		dictHere.emplace(3, std::vector<Modification*> {new Modification(_originalId: L"21", _modificationType: L"myModType", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 21.981943)});
		Protein *ParentProteinToNotInclude = new Protein(L"MPEPTIDEK", L"accession2", L"organism", std::vector<std::tuple<std::wstring, std::wstring>>(), dictHere);
		digestedList = ParentProteinToNotInclude->Digest(task1->getCommonParameters()->getDigestionParams(), fixedModifications, variableModifications).ToList();
		Assert::AreEqual(4, digestedList.size());

		MsDataFile *myMsDataFile1 = new TestDataFile(std::vector<std::vector<PeptideWithSetModifications*>>(1) });

		std::wstring mzmlName1 = LR"(ok1.mzML)";
		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName1, false);

		MsDataFile *myMsDataFile2 = new TestDataFile(std::vector<std::vector<PeptideWithSetModifications*>>(1) });

		std::wstring mzmlName2 = LR"(ok2.mzML)";
		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile2, mzmlName2, false);

		Protein *proteinWithChain1 = new Protein(L"MAACNNNCAA", L"accession3", L"organism", std::vector<std::tuple<std::wstring, std::wstring>>(), std::unordered_map<int, std::vector<Modification*>>(), {new ProteolysisProduct(4, 8, L"chain")}, L"name2", L"fullname2", false, false, std::vector<DatabaseReference*>(), std::vector<SequenceVariation*>(), nullptr);
		Protein *proteinWithChain2 = new Protein(L"MAACNNNCAA", L"accession3", L"organism", std::vector<std::tuple<std::wstring, std::wstring>>(), std::unordered_map<int, std::vector<Modification*>>(), {new ProteolysisProduct(4, 8, L"chain")}, L"name2", L"fullname2", false, false, std::vector<DatabaseReference*>(), std::vector<SequenceVariation*>(), nullptr);

		std::wstring xmlName = L"okk.xml";
		ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>(), {ParentProtein, proteinWithChain1, proteinWithChain2}, xmlName);

		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestMultipleFilesRunner)");
		// RUN!
		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mzmlName1, mzmlName2}, std::vector<DbForTask*> {new DbForTask(xmlName, false)}, outputFolder);
		engine->Run();
		Directory::Delete(outputFolder, true);
		File::Delete(xmlName);
		File::Delete(mzmlName1);
		File::Delete(mzmlName2);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete proteinWithChain2' statement was not added since proteinWithChain2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete proteinWithChain1' statement was not added since proteinWithChain1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile2' statement was not added since myMsDataFile2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile1' statement was not added since myMsDataFile1 was passed to a method or constructor. Handle memory management manually.
		delete ParentProteinToNotInclude;
//C# TO C++ CONVERTER TODO TASK: A 'delete ParentProtein' statement was not added since ParentProtein was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task4' statement was not added since task4 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task3' statement was not added since task3 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task2' statement was not added since task2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task1' statement was not added since task1 was passed to a method or constructor. Handle memory management manually.
	}

	void MyTaskTest::MakeSureFdrDoesntSkip()
	{
		MetaMorpheusTask *task = new SearchTask();
		CommonParameters tempVar(, , , , 999, , , , , 1, , , , , , , , , new PpmTolerance(50), , new DigestionParams(minPeptideLength: 2));
		task->setCommonParameters(&tempVar);
		task->SearchParameters = new SearchParameters();
		task->SearchParameters->DecoyType = DecoyType::None;
		task->SearchParameters->MassDiffAcceptorType = MassDiffAcceptorType::Open;

		std::wstring xmlName = L"MakeSureFdrDoesntSkip.xml";

		{
			Protein *theProtein = new Protein(L"MG", L"accession1");
			ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>(), {theProtein}, xmlName);

//C# TO C++ CONVERTER TODO TASK: A 'delete theProtein' statement was not added since theProtein was passed to a method or constructor. Handle memory management manually.
		}

		std::wstring mzmlName = LR"(MakeSureFdrDoesntSkip.mzML)";


		std::unordered_map<std::wstring, Modification*> ok;
		auto theProteins = ProteinDbLoader::LoadProteinXML(xmlName, true, DecoyType::Reverse, std::vector<Modification*>(), false, std::vector<std::wstring>(), ok);

		std::vector<Modification*> fixedModifications;

		auto targetDigested = theProteins[0].Digest(task->getCommonParameters()->getDigestionParams(), fixedModifications, GlobalVariables::getAllModsKnown().OfType<Modification*>().ToList()).ToList();

		PeptideWithSetModifications *targetGood = targetDigested.front();

		TestDataFile *myMsDataFile = new TestDataFile({targetGood}, true);

		auto ii = myMsDataFile->GetOneBasedScan(1)->MassSpectrum.YArray.ToList();

		ii.push_back(1);
		ii.push_back(1);
		ii.push_back(1);
		ii.push_back(1);

		auto intensities = ii.ToArray();

		auto mm = myMsDataFile->GetOneBasedScan(1)->MassSpectrum.XArray.ToList();

		auto hah = 104.35352;
		mm.push_back(hah);
		mm.push_back(hah + 1);
		mm.push_back(hah + 2);

		auto mz = mm.ToArray();

		Array::Sort(mz, intensities);

		myMsDataFile->ReplaceFirstScanArrays(mz, intensities);

		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestMakeSureFdrDoesntSkip)");
		FileSystem::createDirectory(outputFolder);

		// RUN!
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		auto theStringResult = task->RunTask(outputFolder, std::vector<DbForTask*> {new DbForTask(xmlName, false)}, std::vector<std::wstring> {mzmlName}, L"taskId1")->ToString();
		Assert::IsTrue(theStringResult.find(L"All target PSMS within 1% FDR: 1") != std::wstring::npos);
		Directory::Delete(outputFolder, true);
		File::Delete(xmlName);
		File::Delete(mzmlName);

//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
		delete task;
	}

	void MyTaskTest::MakeSureGptmdTaskMatchesExactMatches()
	{
		MetaMorpheusTask *task1;

		{
			ModificationMotif motif;
			ModificationMotif::TryGetMotif(L"T", motif);
			Modification *myNewMod = new Modification(_originalId: L"ok", _modificationType: L"okType", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 229);

			GlobalVariables::AddMods({myNewMod}, false);
			task1 = new GptmdTask();
			CommonParameters tempVar(, , , , , = 12, = true, = false, = 1, 1, , , , , , , , new AbsoluteTolerance(1), , , new DigestionParams(initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain), new std::vector<(std::wstring, std::wstring)*>(), new std::vector<(std::wstring, std::wstring)*>());
			task1->setCommonParameters(&tempVar);
			task1->GptmdParameters = new GptmdParameters();
			task1->GptmdParameters->ListOfModsGptmd = {(L"okType", L"ok on T")};

//C# TO C++ CONVERTER TODO TASK: A 'delete myNewMod' statement was not added since myNewMod was passed to a method or constructor. Handle memory management manually.
		}

		std::wstring xmlName = L"sweetness.xml";

		{
			Protein *theProtein = new Protein(L"MPEPTIDEKANTHE", L"accession1");
			ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>(), {theProtein}, xmlName);

//C# TO C++ CONVERTER TODO TASK: A 'delete theProtein' statement was not added since theProtein was passed to a method or constructor. Handle memory management manually.
		}

		std::wstring mzmlName = LR"(ok.mzML)";

		{
			std::unordered_map<std::wstring, Modification*> ok;
			auto theProteins = ProteinDbLoader::LoadProteinXML(xmlName, true, DecoyType::Reverse, std::vector<Modification*>(), false, std::vector<std::wstring>(), ok);

			std::vector<Modification*> fixedModifications;

			auto targetDigested = theProteins[0].Digest(task1->getCommonParameters()->getDigestionParams(), fixedModifications, GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
			{
				b::OriginalId->Equals(L"ok");
			}).ToList()).ToList();

			ModificationMotif motif;
			ModificationMotif::TryGetMotif(L"T", motif);
			PeptideWithSetModifications *targetGood = targetDigested[0];

			PeptideWithSetModifications *targetWithUnknownMod = targetDigested[1];
			MsDataFile *myMsDataFile = new TestDataFile({targetGood, targetWithUnknownMod}, true);

			IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
		}
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestMakeSureGptmdTaskMatchesExactMatchesTest)");
		FileSystem::createDirectory(outputFolder);

		// RUN!
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		auto theStringResult = task1->RunTask(outputFolder, std::vector<DbForTask*> {new DbForTask(xmlName, false)}, std::vector<std::wstring> {mzmlName}, L"taskId1")->ToString();
		Assert::IsTrue(theStringResult.find(L"Modifications added: 1") != std::wstring::npos);
		Directory::Delete(outputFolder, true);
		File::Delete(xmlName);
		File::Delete(mzmlName);
		Directory::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(Task Settings)"), true);

		delete task1;
	}

	void MyTaskTest::TestPeptideCount()
	{
		SearchTask *testPeptides = new SearchTask();
		CommonParameters tempVar(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, new DigestionParams(minPeptideLength: 5));
		testPeptides->setCommonParameters(&tempVar);
		SearchParameters tempVar2();
		testPeptides->setSearchParameters(&tempVar2);
		testPeptides->getSearchParameters()->setWritePrunedDatabase(true);
		testPeptides->getSearchParameters()->setSearchTarget(true);
		testPeptides->getSearchParameters()->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);

		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"TestPeptides", testPeptides)};

		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"P", motif);

		auto testUniqeMod = new Modification(_originalId: L"testPeptideMod", _modificationType: L"mt", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);
		GlobalVariables::AddMods({testUniqeMod}, false);

		//create modification lists

		std::vector<Modification*> variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			testPeptides->getCommonParameters()->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();

		//add modification to Protein object
		auto modDictionary = std::unordered_map<int, std::vector<Modification*>>();
		Modification *modToAdd = testUniqeMod;
		modDictionary.emplace(1, std::vector<Modification*> {modToAdd});
		modDictionary.emplace(3, std::vector<Modification*> {modToAdd});

		//protein Creation (One with mod and one without)
		Protein *TestProtein = new Protein(L"PEPTID", L"accession1", L"organism", std::vector<std::tuple<std::wstring, std::wstring>>(), modDictionary);

		//First Write XML Database

		std::wstring xmlName = L"singleProteinWithTwoMods.xml";

		//Add Mod to list and write XML input database
		std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> modList;
		auto Hash = std::unordered_set<std::tuple<int, Modification*>> {std::tuple<int, Modification*>(3, modToAdd)};
		modList.emplace(L"test", Hash);
		ProteinDbWriter::WriteXmlDatabase(modList, {TestProtein}, xmlName);

		//now write MZML file
		std::unordered_map<std::wstring, Modification*> ok;
		auto protein = ProteinDbLoader::LoadProteinXML(xmlName, true, DecoyType::Reverse, std::vector<Modification*>(), false, std::vector<std::wstring>(), ok);
		auto setList1 = protein[0].Digest(testPeptides->getCommonParameters()->getDigestionParams(), std::vector<Modification*> { }, variableModifications).ToList();
		Assert::AreEqual(4, setList1.size());

		//Finally Write MZML file
		MsDataFile *myMsDataFile = new TestDataFile(std::vector<std::vector<PeptideWithSetModifications*>>(0), setList1[1], setList1[2], setList1[3], setList1[0], setList1[1] });
		std::wstring mzmlName = LR"(singleProteinWithRepeatedMods.mzML)";
		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestMultipleFilesRunner)");
		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mzmlName}, std::vector<DbForTask*> {new DbForTask(xmlName, false)}, outputFolder);
		engine->Run();

		std::wstring line;

		bool foundD = false;
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamReader file = new StreamReader(Path.Combine(MySetUpClass.outputFolder, "TestPeptides", "results.txt")))
		{
			StreamReader file = StreamReader(FileSystem::combine(MySetUpClass::outputFolder, L"TestPeptides", L"results.txt"));
			while ((line = file.ReadLine()) != L"")
			{
				if (line.find(L"All target peptides within 1% FDR: 4") != std::wstring::npos)
				{
					foundD = true;
				}
			}
		}
		Assert::IsTrue(foundD);
		Directory::Delete(outputFolder, true);
		File::Delete(mzmlName);
		File::Delete(xmlName);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete TestProtein' statement was not added since TestProtein was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete testUniqeMod' statement was not added since testUniqeMod was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete testPeptides' statement was not added since testPeptides was passed to a method or constructor. Handle memory management manually.
	}

	void MyTaskTest::TestFileOutput()
	{
		std::wstring thisTaskOutputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestFileOutput)");

		SearchTask *task = Toml::ReadFile<SearchTask*>(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(SlicedSearchTaskConfig.toml)"), MetaMorpheusTask::tomlConfig);
		task->getSearchParameters()->setDecoyType(DecoyType::None);

		DbForTask *db = new DbForTask(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(sliced-db.fasta)"), false);
		DbForTask *db2 = new DbForTask(FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestData", LR"(DbForPrunedDb.fasta)"), false);
		std::wstring raw = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(sliced-raw.mzML)");
		std::wstring raw2 = FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestData", LR"(PrunedDbSpectra.mzml)");
		EverythingRunnerEngine *singleMassSpectraFile = new EverythingRunnerEngine({(L"SingleMassSpectraFileOutput", task)}, {raw}, {db}, thisTaskOutputFolder);
		EverythingRunnerEngine *multipleMassSpectraFiles = new EverythingRunnerEngine({(L"MultipleMassSpectraFileOutput", task)}, {raw, raw2}, {db, db2}, thisTaskOutputFolder);

		singleMassSpectraFile->Run();
		multipleMassSpectraFiles->Run();

		// test single file output
		std::unordered_set<std::wstring> expectedFiles = {L"AllPeptides.psmtsv", L"AllProteinGroups.tsv", L"AllPSMs.psmtsv", L"AllPSMs_FormattedForPercolator.tsv", L"AllQuantifiedPeaks.tsv", L"AllQuantifiedPeptides.tsv", L"prose.txt", L"results.txt"};

		std::unordered_set<std::wstring> files = std::unordered_set<std::wstring>(Directory::GetFiles(FileSystem::combine(thisTaskOutputFolder, L"SingleMassSpectraFileOutput")).Select([&] (std::any v)
		{
			FileSystem::getFileName(v);
		}));

		// these 2 lines are for debug purposes, so you can see which files you're missing (if any)
		auto missingFiles = expectedFiles.Except(files);
		auto extraFiles = files.Except(expectedFiles);

		// test that output is what's expected
		Assert::That(files.SetEquals(expectedFiles));

		// test multi file output
		files = std::unordered_set<std::wstring>(Directory::GetFiles(FileSystem::combine(thisTaskOutputFolder, L"MultipleMassSpectraFileOutput")).Select([&] (std::any v)
		{
			FileSystem::getFileName(v);
		}));
		missingFiles = expectedFiles.Except(files);
		extraFiles = files.Except(expectedFiles);

		Assert::That(files.SetEquals(expectedFiles));

		expectedFiles = {L"PrunedDbSpectra.mzID", L"PrunedDbSpectra_PSMs.psmtsv", L"PrunedDbSpectra_PSMsFormattedForPercolator.tsv", L"PrunedDbSpectra_Peptides.psmtsv", L"PrunedDbSpectra_ProteinGroups.tsv", L"PrunedDbSpectra_QuantifiedPeaks.tsv", L"sliced-raw.mzID", L"sliced-raw_PSMs.psmtsv", L"sliced-raw_PSMsFormattedForPercolator.tsv", L"sliced-raw_Peptides.psmtsv", L"sliced-raw_ProteinGroups.tsv", L"sliced-raw_QuantifiedPeaks.tsv"};

		std::wstring individualFilePath = FileSystem::combine(thisTaskOutputFolder, L"MultipleMassSpectraFileOutput", L"Individual File Results");
		Assert::That(FileSystem::directoryExists(individualFilePath));

		files = std::unordered_set<std::wstring>(Directory::GetFiles(individualFilePath).Select([&] (std::any v)
		{
			FileSystem::getFileName(v);
		}));
		missingFiles = expectedFiles.Except(files);
		extraFiles = files.Except(expectedFiles);

		Assert::That(files.SetEquals(expectedFiles));

		files = std::unordered_set<std::wstring>(Directory::GetFiles(FileSystem::combine(thisTaskOutputFolder, L"Task Settings")).Select([&] (std::any v)
		{
			FileSystem::getFileName(v);
		}));
		expectedFiles = {L"MultipleMassSpectraFileOutputconfig.toml", L"SingleMassSpectraFileOutputconfig.toml"};
		Assert::That(files.SetEquals(expectedFiles));
		Directory::Delete(thisTaskOutputFolder, true);

		delete multipleMassSpectraFiles;
		delete singleMassSpectraFile;
//C# TO C++ CONVERTER TODO TASK: A 'delete db2' statement was not added since db2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete db' statement was not added since db was passed to a method or constructor. Handle memory management manually.
	}

	void MyTaskTest::TestUniprotNamingConflicts()
	{
		// write the mod
		auto outputDir = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestUniprotNamingConflicts)");
		FileSystem::createDirectory(outputDir);
		std::wstring modToWrite = std::wstring(L"Custom List\nID   Hydroxyproline\nTG   P\nPP   Anywhere.\nMT   Biological\nCF   O1\n") + LR"(//)";
		auto filePath = FileSystem::combine(GlobalVariables::getDataDir(), LR"(Mods)", LR"(hydroxyproline.txt)");
		File::WriteAllLines(filePath, std::vector<std::wstring> {modToWrite});

		// read the mod
		std::any fmww;
		GlobalVariables::AddMods(PtmListLoader::ReadModsFromFile(filePath, fmww), false);
		Assert::That(GlobalVariables::getAllModsKnown().Where([&] (std::any v)
		{
			return v->IdWithMotif == L"Hydroxyproline on P";
		})->Count() == 1);

		// should have an error message...
		Assert::That(GlobalVariables::ErrorsReadingMods.Where([&] (std::any v)
		{
			v->Contains(L"Hydroxyproline");
		})->Count() > 0);
		Directory::Delete(outputDir, true);
	}

	void MyTaskTest::TestPepXmlOutput()
	{
		SearchTask *search = new SearchTask();
		SearchParameters tempVar();
		search->setSearchParameters(&tempVar);
		search->getSearchParameters()->setWritePepXml(true);

		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"TestPepXmlOutput", search)};

		std::wstring mzmlName = LR"(TestData\PrunedDbSpectra.mzml)";
		std::wstring fastaName = LR"(TestData\DbForPrunedDb.fasta)";
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestPepXmlOutput)");

		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mzmlName}, std::vector<DbForTask*> {new DbForTask(fastaName, false)}, outputFolder);
		engine->Run();

		std::wstring outputPepXmlPath = FileSystem::combine(outputFolder, LR"(TestPepXmlOutput\Individual File Results\PrunedDbSpectra.pep.XML)");
		Assert::That(FileSystem::fileExists(outputPepXmlPath));
		Directory::Delete(outputFolder, true);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete search' statement was not added since search was passed to a method or constructor. Handle memory management manually.
	}

	void MyTaskTest::TestModernAndClassicSearch()
	{
		SearchTask *classicSearch = new SearchTask();

		SearchTask *modernSearch = new SearchTask();
		SearchParameters tempVar();
		modernSearch->setSearchParameters(&tempVar);
		modernSearch->getSearchParameters()->setSearchType(SearchType::Modern);
		std::vector<int> counts;

		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"ClassicSearch", classicSearch), (L"ModernSearch", modernSearch)};

		std::wstring mzmlName = LR"(TestData\PrunedDbSpectra.mzml)";
		std::wstring fastaName = LR"(TestData\DbForPrunedDb.fasta)";
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestPepXmlOutput)");

		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mzmlName}, std::vector<DbForTask*> {new DbForTask(fastaName, false)}, outputFolder);
		engine->Run();

		std::wstring classicPath = FileSystem::combine(outputFolder, LR"(ClassicSearch\AllPSMs.psmtsv)");
		auto classicPsms = File::ReadAllLines(classicPath).ToList();

		std::wstring modernPath = FileSystem::combine(outputFolder, LR"(ModernSearch\AllPSMs.psmtsv)");
		auto modernPsms = File::ReadAllLines(modernPath).ToList();
		counts.push_back(modernPsms.size());

		Assert::That(modernPsms.SequenceEqual(classicPsms));
		Directory::Delete(outputFolder, true);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete modernSearch' statement was not added since modernSearch was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete classicSearch' statement was not added since classicSearch was passed to a method or constructor. Handle memory management manually.
	}
}
