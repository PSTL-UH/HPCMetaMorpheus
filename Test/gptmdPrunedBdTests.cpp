#include "gptmdPrunedBdTests.h"
#include "../TaskLayer/GPTMDTask/GPTMDTask.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"
#include "SetUpTests.h"
#include "../EngineLayer/GlobalVariables.h"
#include "TestDataFile.h"

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{

	void GptmdPrunedDbTests::TestPrunedGeneration()
	{
		//Create GPTMD Task
		//Create Search Task
		GptmdTask *task1 = new GptmdTask();
		CommonParameters tempVar();
		task1->setCommonParameters(&tempVar);

		SearchTask *task2 = new SearchTask();
		CommonParameters tempVar2();
		task2->setCommonParameters(&tempVar2);
		SearchParameters tempVar3();
		task2->setSearchParameters(&tempVar3);
		task2->getSearchParameters()->setDoParsimony(true);
		task2->getSearchParameters()->setSearchTarget(true);
		task2->getSearchParameters()->setWritePrunedDatabase(true);
		task2->getSearchParameters()->setSearchType(SearchType::Classic);
		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"task1", task1), (L"task2", task2)};
		std::wstring mzmlName = LR"(TestData\PrunedDbSpectra.mzml)";
		std::wstring fastaName = LR"(TestData\DbForPrunedDb.fasta)";
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestPrunedGeneration)");
		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mzmlName}, std::vector<DbForTask*> {new DbForTask(fastaName, false)}, outputFolder);
		engine->Run();
		std::wstring final_Renamed = FileSystem::combine(MySetUpClass::outputFolder, L"task2", L"DbForPrunedDbGPTMDproteinPruned.xml");
		std::any ok;
		std::vector<Protein*> proteins = ProteinDbLoader::LoadProteinXML(final_Renamed, true, DecoyType::Reverse, std::vector<Modification*>(), false, std::vector<std::wstring>(), ok);
		//ensures that protein out put contins the correct number of proteins to match the folowing conditions. 
			// all proteins in DB have baseSequence!=null (not ambiguous)
			// all proteins that belong to a protein group are written to DB
		Assert::AreEqual(20, proteins.size());
		int totalNumberOfMods = proteins.Sum([&] (std::any p)
		{
		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete task2' statement was not added since task2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task1' statement was not added since task1 was passed to a method or constructor. Handle memory management manually.
			return p::OneBasedPossibleLocalizedModifications->Count + p::SequenceVariations::Sum([&] (std::any sv)
			{
				sv::OneBasedModifications->Count;
			});
		});

		//tests that modifications are being done correctly
		Assert::AreEqual(0, totalNumberOfMods);
		Directory::Delete(outputFolder, true);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete task2' statement was not added since task2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task1' statement was not added since task1 was passed to a method or constructor. Handle memory management manually.
	}

	void GptmdPrunedDbTests::TestPrunedDatabase()
	{
		//Create Search Task
		SearchTask *task1 = new SearchTask();
		SearchParameters tempVar();
		task1->setSearchParameters(&tempVar);
		task1->getSearchParameters()->setWritePrunedDatabase(true);
		task1->getSearchParameters()->setSearchTarget(true);
		task1->getSearchParameters()->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
		task1->getSearchParameters()->setModsToWriteSelection(std::unordered_map<std::wstring, int>
		{
			{L"ConnorModType", 1}
		});
		CommonParameters tempVar2(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, new DigestionParams(minPeptideLength: 5));
		task1->setCommonParameters(&tempVar2);

		//add task to task list
		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"task1", task1)};

		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"P", motif);

		auto connorMod = new Modification(_originalId: L"ConnorMod on P", _modificationType: L"ConnorModType", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);

		GlobalVariables::AddMods({connorMod}, false);

		//create modification lists
		std::vector<Modification*> variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			task1->getCommonParameters()->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();

		//add modification to Protein object
		auto dictHere = std::unordered_map<int, std::vector<Modification*>>();
		Modification *modToAdd = connorMod;
		Modification *modToAdd2 = connorMod;
		dictHere.emplace(1, std::vector<Modification*> {modToAdd});
		dictHere.emplace(3, std::vector<Modification*> {modToAdd2});

		//protein Creation (One with mod and one without)
		Protein *TestProteinWithMod = new Protein(L"PEPTID", L"accession1", L"organism", std::vector<std::tuple<std::wstring, std::wstring>>(), dictHere);

		//First Write XML Database
		std::wstring xmlName = L"okkk.xml";

		//Add Mod to list and write XML input database
		std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> modList;
		auto Hash = std::unordered_set<std::tuple<int, Modification*>> {std::tuple<int, Modification*>(3, modToAdd)};
		modList.emplace(L"test", Hash);
		ProteinDbWriter::WriteXmlDatabase(modList, {TestProteinWithMod}, xmlName);

		//now write MZML file
		std::unordered_map<std::wstring, Modification*> ok;
		auto protein = ProteinDbLoader::LoadProteinXML(xmlName, true, DecoyType::Reverse, std::vector<Modification*>(), false, std::vector<std::wstring>(), ok);

		//Dictionary 'ok' contains unknown modifications. There are no unknown modifications in this test.
		Assert::AreEqual(0, ok->Count);
		//One protein is read from the .xml database and one decoy is created. Therefore, the list of proteins contains 2 entries.
		Assert::AreEqual(2, protein->Count);
		//The original database had two localized mods on the protein. Therefore. both protein and decoy should have two mods.
		Assert::AreEqual(2, protein[0].OneBasedPossibleLocalizedModifications->Count);
		std::vector<int> foundResidueIndicies = protein[0].OneBasedPossibleLocalizedModifications->Select([&] (std::any k)
		{
			k::Key;
		}).ToList();
		std::vector<int> expectedResidueIndices = {1, 3};
		Assert::That(foundResidueIndicies, Is::EquivalentTo(expectedResidueIndices));
		Assert::AreEqual(2, protein[1].OneBasedPossibleLocalizedModifications->Count);
		foundResidueIndicies = protein[1].OneBasedPossibleLocalizedModifications->Select([&] (std::any k)
		{
			k::Key;
		}).ToList();
		expectedResidueIndices = {4, 6}; //originally modified residues are now at the end in the decoy
		Assert::That(foundResidueIndicies, Is::EquivalentTo(expectedResidueIndices));


		auto thisOk = ok; //for debugging
		auto commonParamsAtThisPoint = task1->getCommonParameters()->getDigestionParams(); //for debugging

		auto digestedList = protein[0].Digest(task1->getCommonParameters()->getDigestionParams(), std::vector<Modification*> { }, variableModifications).ToList();
		Assert::AreEqual(4, digestedList.size());

		//Set Peptide with 1 mod at position 3
		PeptideWithSetModifications *pepWithSetMods1 = digestedList[1];

		//Finally Write MZML file
		Assert::AreEqual(L"PEP[ConnorModType:ConnorMod on P]TID", pepWithSetMods1->FullSequence); //this might be base sequence
		MsDataFile *myMsDataFile = new TestDataFile(std::vector<PeptideWithSetModifications*> {pepWithSetMods1});
		std::wstring mzmlName = LR"(hello.mzML)";
		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

		//run!
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestPrunedDatabase)");
		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mzmlName}, std::vector<DbForTask*> {new DbForTask(xmlName, false)}, outputFolder);
		engine->Run();

		std::wstring final_Renamed = FileSystem::combine(MySetUpClass::outputFolder, L"task1", L"okkkpruned.xml");

		auto proteins = ProteinDbLoader::LoadProteinXML(final_Renamed, true, DecoyType::Reverse, std::vector<Modification*>(), false, std::vector<std::wstring>(), ok);
		//check length
		Assert::AreEqual(1, proteins[0].OneBasedPossibleLocalizedModifications->Count);
		//check location (key)
		Assert::AreEqual(true, proteins[0].OneBasedPossibleLocalizedModifications->ContainsKey(3));
		std::vector<Modification*> listOfMods = proteins[0].OneBasedPossibleLocalizedModifications[3];
		//check Type, count, ID
		Assert::AreEqual(listOfMods[0]->ModificationType, L"ConnorModType");
		Assert::AreEqual(listOfMods[0]->IdWithMotif, L"ConnorMod on P");
		Assert::AreEqual(listOfMods.size(), 1);
		Directory::Delete(outputFolder, true);
		File::Delete(xmlName);
		File::Delete(mzmlName);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete TestProteinWithMod' statement was not added since TestProteinWithMod was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete connorMod' statement was not added since connorMod was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task1' statement was not added since task1 was passed to a method or constructor. Handle memory management manually.
	}

	void GptmdPrunedDbTests::TestUserModSelectionInPrunedDB()
	{
		std::vector<(std::wstring, std::wstring)*> listOfModsFixed = {(L"Common Fixed", L"Carbamidomethyl of C"), (L"Common Fixed", L"Carbamidomethyl of U")};
		//Create Search Task
		SearchTask *task5 = new SearchTask();
		SearchParameters tempVar();
		task5->setSearchParameters(&tempVar);
		task5->getSearchParameters()->setWritePrunedDatabase(true);
		task5->getSearchParameters()->setSearchTarget(true);
		task5->getSearchParameters()->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
		CommonParameters tempVar2(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, = nullptr, = nullptr, listOfModsFixed);
		task5->setCommonParameters(&tempVar2);

		task5->getSearchParameters()->getModsToWriteSelection()[L"Mod"] = 0;
		task5->getSearchParameters()->getModsToWriteSelection()[L"Common Fixed"] = 1;
		task5->getSearchParameters()->getModsToWriteSelection()[L"Glycan"] = 2;
		task5->getSearchParameters()->getModsToWriteSelection()[L"missing"] = 3;

		//add task 1 to task list
		std::vector<(std::wstring, MetaMorpheusTask)*> taskList = {(L"task5", task5)};
		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"P", motif);
		ModificationMotif motif2;
		ModificationMotif::TryGetMotif(L"E", motif2);

		auto connorMod = new Modification(_originalId: L"ModToNotAppear", _modificationType: L"Mod", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);
		auto connorMod2 = new Modification(_originalId: L"Default(Mod in DB and Observed)", _modificationType: L"Common Fixed", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);
		auto connorMod3 = new Modification(_originalId: L"ModToAlwaysAppear", _modificationType: L"Glycan", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);
		auto connorMod4 = new Modification(_originalId: L"ModObservedNotinDB", _modificationType: L"missing", _target: motif2, _locationRestriction: L"Anywhere.", _monoisotopicMass: 5);

		GlobalVariables::AddMods({connorMod, connorMod2, connorMod3, connorMod4}, false);

		//create modification lists
		std::vector<Modification*> variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			task5->getCommonParameters()->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();
		std::vector<Modification*> fixedModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			task5->getCommonParameters()->ListOfModsFixed->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();

		//add modification to Protein object
		auto dictHere = std::unordered_map<int, std::vector<Modification*>>();
		Modification *modToAdd = connorMod;
		Modification *modToAdd2 = connorMod2;
		Modification *modToAdd3 = connorMod3;
		Modification *modToAdd4 = connorMod4;

		//add Fixed modifcation so can test if mod that is observed and not in DB
		fixedModifications.push_back(connorMod4);
		listOfModsFixed.push_back((connorMod4->ModificationType, connorMod4->IdWithMotif));

		dictHere.emplace(1, std::vector<Modification*> {modToAdd});
		dictHere.emplace(2, std::vector<Modification*> {modToAdd2}); //default
		dictHere.emplace(3, std::vector<Modification*> {modToAdd3}); //Alway Appear

		auto dictHere2 = std::unordered_map<int, std::vector<Modification*>>();
		dictHere2.emplace(1, std::vector<Modification*> {modToAdd});
		dictHere2.emplace(2, std::vector<Modification*> {modToAdd2}); //default
		dictHere2.emplace(3, std::vector<Modification*> {modToAdd3}); //Alway Appear
		dictHere2.emplace(4, std::vector<Modification*> {modToAdd4}); //observed
		//protein Creation (One with mod and one without)
		Protein *TestProteinWithModForDB = new Protein(L"PPPPPPPPPPE", L"accession1", L"organism", std::vector<std::tuple<std::wstring, std::wstring>>(), dictHere);
		Protein *TestProteinWithModObsevred = new Protein(L"PPPPPPPPPPE", L"accession1", L"organism", std::vector<std::tuple<std::wstring, std::wstring>>(), dictHere2);

		//First Write XML Database
		std::wstring xmlName = L"selectedMods.xml";
		std::wstring xmlName2 = L"selectedModsObvs.xml";

		//Add Mod to list and write XML input database
		std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> modList;
		auto Hash = std::unordered_set<std::tuple<int, Modification*>> {std::tuple<int, Modification*>(1, modToAdd), std::tuple<int, Modification*>(2, modToAdd2), std::tuple<int, Modification*>(3, modToAdd3), std::tuple<int, Modification*>(4, modToAdd4)};

		modList.emplace(L"test", Hash);
		ProteinDbWriter::WriteXmlDatabase(modList, {TestProteinWithModForDB}, xmlName);

		//Add Observed Only
		modList.emplace(L"test2", Hash);
		ProteinDbWriter::WriteXmlDatabase(modList, {TestProteinWithModObsevred}, xmlName2);

		//now create MZML data
		std::unordered_map<std::wstring, Modification*> ok;
		auto protein = ProteinDbLoader::LoadProteinXML(xmlName2, true, DecoyType::Reverse, std::vector<Modification*>(), false, std::vector<std::wstring>(), ok);
		auto digestedList = protein[0].Digest(task5->getCommonParameters()->getDigestionParams(), fixedModifications, variableModifications).ToList();

		//Set Peptide with 1 mod at position 3
		PeptideWithSetModifications *pepWithSetMods1 = digestedList[0];
		PeptideWithSetModifications *pepWithSetMods2 = digestedList[1];
		PeptideWithSetModifications *pepWithSetMods3 = digestedList[2];
		PeptideWithSetModifications *pepWithSetMods4 = digestedList[3];
		PeptideWithSetModifications *pepWithSetMods5 = digestedList[4];

		//CUSTOM PEP
		MsDataFile *myMsDataFile = new TestDataFile(std::vector<PeptideWithSetModifications*> {pepWithSetMods1, pepWithSetMods2, pepWithSetMods3, pepWithSetMods4, pepWithSetMods5});
		std::wstring mzmlName = LR"(newMzml.mzML)";
		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

		//make sure this runs correctly
		//run!
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestUserModSelectionInPrunedDB)");
		auto engine = new EverythingRunnerEngine(taskList, std::vector<std::wstring> {mzmlName}, std::vector<DbForTask*> {new DbForTask(xmlName, false)}, outputFolder);
		engine->Run();
		std::wstring final_Renamed = FileSystem::combine(MySetUpClass::outputFolder, L"task5", L"selectedModspruned.xml");
		auto proteins = ProteinDbLoader::LoadProteinXML(final_Renamed, true, DecoyType::Reverse, std::vector<Modification*>(), false, std::vector<std::wstring>(), ok);
		auto Dlist = proteins[0].GetVariantProteins().SelectMany([&] (std::any vp)
		{
			vp::Digest(task5->getCommonParameters()->getDigestionParams(), fixedModifications, variableModifications);
		}).ToList();
		Assert::AreEqual(Dlist[0].NumFixedMods, 1);

		//check length
		Assert::AreEqual(proteins[0].OneBasedPossibleLocalizedModifications->Count, 3);
		std::vector<Modification*> listOfLocalMods;
		listOfLocalMods.insert(listOfLocalMods.end(), (proteins[0].OneBasedPossibleLocalizedModifications[2]).begin(), (proteins[0].OneBasedPossibleLocalizedModifications[2]).end());
		listOfLocalMods.insert(listOfLocalMods.end(), (proteins[0].OneBasedPossibleLocalizedModifications[3]).begin(), (proteins[0].OneBasedPossibleLocalizedModifications[3]).end());
		listOfLocalMods.insert(listOfLocalMods.end(), (proteins[0].OneBasedPossibleLocalizedModifications[11]).begin(), (proteins[0].OneBasedPossibleLocalizedModifications[11]).end());

		//check Type, count, ID
		Assert::AreEqual(listOfLocalMods[0]->ModificationType, L"Common Fixed");
		Assert::AreEqual(listOfLocalMods[2]->ModificationType, L"missing");
		Assert::IsFalse(std::find(listOfLocalMods.begin(), listOfLocalMods.end(), connorMod) != listOfLocalMods.end()); //make sure that mod set not to show up is not in mod list

		Assert::AreEqual(listOfLocalMods[0]->IdWithMotif, L"Default(Mod in DB and Observed) on P");
		Assert::AreEqual(listOfLocalMods[1]->IdWithMotif, L"ModToAlwaysAppear on P");
		//Makes sure Mod that was not in the DB but was observed is in pruned DB
		Assert::AreEqual(listOfLocalMods[2]->IdWithMotif, L"ModObservedNotinDB on E");
		Assert::AreEqual(listOfLocalMods.size(), 3);
		Directory::Delete(outputFolder, true);
		File::Delete(mzmlName);
		File::Delete(xmlName);
		File::Delete(xmlName2);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete TestProteinWithModObsevred' statement was not added since TestProteinWithModObsevred was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete TestProteinWithModForDB' statement was not added since TestProteinWithModForDB was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete connorMod4' statement was not added since connorMod4 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete connorMod3' statement was not added since connorMod3 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete connorMod2' statement was not added since connorMod2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete connorMod' statement was not added since connorMod was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete task5' statement was not added since task5 was passed to a method or constructor. Handle memory management manually.
	}
}
