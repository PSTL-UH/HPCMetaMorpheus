#include "BinGenerationTest.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "TestDataFile.h"
#include "../TaskLayer/DbForTask.h"

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{

	void BinGenerationTest::TestBinGeneration()
	{
		SearchTask *st = new SearchTask();
		CommonParameters tempVar(, , = true, = true, = 3, = 12, = true, = false, = 1, 1, , , , , , , , , , , new DigestionParams(minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain));
		st->setCommonParameters(&tempVar);
		SearchParameters tempVar2();
		st->setSearchParameters(&tempVar2);
		st->getSearchParameters()->setDoHistogramAnalysis(true);
		st->getSearchParameters()->setMassDiffAcceptorType(MassDiffAcceptorType::Open);
		st->getSearchParameters()->setDecoyType(DecoyType::None);
		st->getSearchParameters()->setDoParsimony(true);
		st->getSearchParameters()->setDoQuantification(true);

		std::wstring proteinDbFilePath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"BinGenerationTest.xml");
		std::wstring mzmlFilePath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"BinGenerationTest.mzML");

		Protein *prot1 = new Protein(L"MEDEEK", L"prot1");
		Protein *prot2 = new Protein(L"MENEEK", L"prot2");

		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"D", motif);
		Modification *mod = new Modification(_target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);

		auto pep1_0 = prot1->Digest(st->getCommonParameters()->getDigestionParams(), std::vector<Modification*>(), std::vector<Modification*>()).First();
		auto pep1_10 = prot1->Digest(st->getCommonParameters()->getDigestionParams(), std::vector<Modification*>(), std::vector<Modification*>()).Last();

		Protein *prot3 = new Protein(L"MAAADAAAAAAAAAAAAAAA", L"prot3");

		auto pep2_0 = prot3->Digest(st->getCommonParameters()->getDigestionParams(), std::vector<Modification*>(), std::vector<Modification*>()).First();
		auto pep2_10 = prot3->Digest(st->getCommonParameters()->getDigestionParams(), std::vector<Modification*>(), std::vector<Modification*> {mod}).Last();

		Protein *prot4 = new Protein(L"MNNDNNNN", L"prot4");
		auto pep3_10 = prot4->Digest(st->getCommonParameters()->getDigestionParams(), std::vector<Modification*>(), std::vector<Modification*> {mod}).Last();

		std::vector<PeptideWithSetModifications*> pepsWithSetMods = {pep1_0, pep1_10, pep2_0, pep2_10, pep3_10};
		MsDataFile *myMsDataFile = new TestDataFile(pepsWithSetMods);

		std::vector<Protein*> proteinList = {prot1, prot2, prot3, prot4};

		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlFilePath, false);
		ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>(), proteinList, proteinDbFilePath);

		std::wstring output_folder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestBinGeneration");
		FileSystem::createDirectory(output_folder);
		st->RunTask(output_folder, {new DbForTask(proteinDbFilePath, false)}, {mzmlFilePath}, L"");

		Assert::AreEqual(3, File::ReadLines(FileSystem::combine(output_folder, LR"(MassDifferenceHistogram.tsv)")).size()());
		Directory::Delete(output_folder, true);
		File::Delete(proteinDbFilePath);
		File::Delete(mzmlFilePath);
		Directory::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(Task Settings)"), true);

//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
		delete prot4;
		delete prot3;
//C# TO C++ CONVERTER TODO TASK: A 'delete mod' statement was not added since mod was passed to a method or constructor. Handle memory management manually.
		delete prot2;
		delete prot1;
		delete st;
	}

	void BinGenerationTest::TestProteinSplitAcrossFiles()
	{
		SearchTask *st = new SearchTask();
		CommonParameters tempVar(, , = true, = true, = 3, = 12, = true, = false, = 1, 1, , , , , , , , , , , new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain));
		st->setCommonParameters(&tempVar);
		SearchParameters tempVar2();
		st->setSearchParameters(&tempVar2);
		st->getSearchParameters()->setDoHistogramAnalysis(true);
		st->getSearchParameters()->setMassDiffAcceptorType(MassDiffAcceptorType::Open);
		st->getSearchParameters()->setMatchBetweenRuns(true);
		st->getSearchParameters()->setDoQuantification(true);

		std::wstring proteinDbFilePath = FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestProteinSplitAcrossFiles.xml");
		std::wstring mzmlFilePath1 = FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestProteinSplitAcrossFiles1.mzML");
		std::wstring mzmlFilePath2 = FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestProteinSplitAcrossFiles2.mzML");

		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"D", motif);
		Modification *mod = new Modification(_originalId: L"mod1 on D", _modificationType: L"mt", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);

		std::unordered_map<int, std::vector<Modification*>> oneBasedModification =
		{
			{
				3, {mod}
			}
		};

		Protein *prot1 = new Protein(L"MEDEEK", L"prot1", oneBasedModifications: oneBasedModification);

		auto pep1 = prot1->Digest(st->getCommonParameters()->getDigestionParams(), std::vector<Modification*>(), std::vector<Modification*>()).First();
		auto pep2 = prot1->Digest(st->getCommonParameters()->getDigestionParams(), std::vector<Modification*>(), std::vector<Modification*>()).Last();

		std::vector<PeptideWithSetModifications*> listForFile1 = {pep1, pep2};
		std::vector<PeptideWithSetModifications*> listForFile2 = {pep2};
		MsDataFile *myMsDataFile1 = new TestDataFile(listForFile1);
		MsDataFile *myMsDataFile2 = new TestDataFile(listForFile2);

		std::vector<Protein*> proteinList = {prot1};

		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlFilePath1, false);
		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile2, mzmlFilePath2, false);
		ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>(), proteinList, proteinDbFilePath);

		std::wstring output_folder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestProteinSplitAcrossFiles");
		FileSystem::createDirectory(output_folder);

		st->RunTask(output_folder, {new DbForTask(proteinDbFilePath, false)}, {mzmlFilePath1, mzmlFilePath2}, L"");
		Directory::Delete(output_folder, true);
		File::Delete(proteinDbFilePath);
		File::Delete(mzmlFilePath1);
		File::Delete(mzmlFilePath2);
		Directory::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(Task Settings)"), true);

//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile2' statement was not added since myMsDataFile2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile1' statement was not added since myMsDataFile1 was passed to a method or constructor. Handle memory management manually.
		delete prot1;
		delete mod;
		delete st;
	}
}
