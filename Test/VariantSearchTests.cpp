#include "VariantSearchTests.h"
#include "../EngineLayer/CommonParameters.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "TestDataFile.h"
#include "../TaskLayer/DbForTask.h"

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

VariantSearchTests::CommonParameters *CommonParameters = new CommonParameters(, , = true, = true, = 3, = 12, = true, = false, = 1, 1, , , , , , , , , , , new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain, maxModsForPeptides: 1));

	void VariantSearchTests::SearchTests(int proteinIdx, int peptideIdx, bool containsVariant, const std::wstring &variantPsmShort, DecoyType *decoyType)
	{
		Stopwatch *stopwatch = new Stopwatch();
		stopwatch->Start();

		// Make sure can run the complete search task when multiple compact peptides may correspond to a single PWSM
		SearchTask *st = new SearchTask();
		SearchParameters tempVar();
		st->setSearchParameters(&tempVar);
		st->getSearchParameters()->setDoParsimony(true);
		st->getSearchParameters()->setDecoyType(decoyType);
		st->getSearchParameters()->setSearchTarget(decoyType == DecoyType::None);
		st->getSearchParameters()->setModPeptidesAreDifferent(false);
		CommonParameters tempVar2(, , , = true, = 3, = 12, = true, = false, = 1, 1, , , , , , , , new PpmTolerance(20), , , new DigestionParams(minPeptideLength: 2));
		st->setCommonParameters(&tempVar2);

		ModificationMotif motifV;
		ModificationMotif::TryGetMotif(L"V", motifV);
		Modification *mv = new Modification(L"mod", nullptr, L"type", nullptr, motifV, L"Anywhere.", nullptr, 42.01, std::unordered_map<std::wstring, std::vector<std::wstring>>(), nullptr, nullptr, nullptr, nullptr, nullptr);
		ModificationMotif motifP;
		ModificationMotif::TryGetMotif(L"P", motifP);
		Modification *mp = new Modification(L"mod", nullptr, L"type", nullptr, motifP, L"Anywhere.", nullptr, 42.01, std::unordered_map<std::wstring, std::vector<std::wstring>>(), nullptr, nullptr, nullptr, nullptr, nullptr);

		Protein(L"MPEPTIDE", L"protein1", sequenceVariations: std::vector<SequenceVariation*> tempVar3 = new Protein(L"MPEPTIDE", L"protein1", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar3::new SequenceVariation(4, 4, L"P", L"V", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", nullptr);
		Protein(L"MPEPTIDE", L"protein2", sequenceVariations: std::vector<SequenceVariation*> tempVar4 = new Protein(L"MPEPTIDE", L"protein2", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar4::new SequenceVariation(4, 5, L"PT", L"KT", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", nullptr);
		Protein(L"MPEPTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*> tempVar5 = new Protein(L"MPEPTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar5::new SequenceVariation(4, 4, L"P", L"PPP", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", nullptr);
		Protein(L"MPEPPPTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*> tempVar6 = new Protein(L"MPEPPPTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar6::new SequenceVariation(4, 6, L"PPP", L"P", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", nullptr);
		Protein(L"MPEPKPKTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*> tempVar7 = new Protein(L"MPEPKPKTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar7::new SequenceVariation(4, 7, L"PKPK", L"PK", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", nullptr);
		Protein(L"MPEPTAIDE", L"protein2", sequenceVariations: std::vector<SequenceVariation*> tempVar8 = new Protein(L"MPEPTAIDE", L"protein2", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar8::new SequenceVariation(4, 5, L"PTA", L"KT", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", nullptr);
		Protein(L"MPEKKAIDE", L"protein2", sequenceVariations: std::vector<SequenceVariation*> tempVar9 = new Protein(L"MPEKKAIDE", L"protein2", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar9::new SequenceVariation(4, 4, L"KKA", L"K", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", nullptr);
		Protein(L"MPEPTIDE", L"protein1", sequenceVariations: std::vector<SequenceVariation*> tempVar10 = new Protein(L"MPEPTIDE", L"protein1", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar10::new SequenceVariation(4, 4, L"P", L"V", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", {mv}::ToList()
		}
		});
		Protein(L"MPEPTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*> tempVar11 = new Protein(L"MPEPTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar11::new SequenceVariation(4, 4, L"P", L"PPP", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", {mp}::ToList()
		}
		});
		Protein(L"MPEPTIDEPEPTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*> tempVar12 = new Protein(L"MPEPTIDEPEPTIDE", L"protein3", sequenceVariations: std::vector<SequenceVariation*>();
		tempVar12::new SequenceVariation(4, 4, L"PTIDEPEPTIDE", L"PPP", LR"(1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30)", nullptr);
		std::vector<Protein*> proteins = {tempVar3), tempVar4), tempVar5), tempVar6), tempVar7), tempVar8), tempVar9), tempVar10), tempVar11), tempVar12)};
		PeptideWithSetModifications *pep = proteins[proteinIdx]->GetVariantProteins().SelectMany([&] (std::any p)
		{
			p::Digest(CommonParameters->getDigestionParams(), nullptr, nullptr);
		}).ToList()[peptideIdx];

		std::wstring xmlName = StringHelper::formatSimple(L"andguiaheov{0}.xml", std::to_wstring(proteinIdx));
		ProteinDbWriter::WriteXmlDatabase(nullptr, std::vector<std::vector<Protein*>>(proteinIdx) }, xmlName);

		std::wstring mzmlName = StringHelper::formatSimple(L"ajgdiv{0}.mzML", std::to_wstring(proteinIdx));
		MsDataFile *myMsDataFile = new TestDataFile({pep}, true);

		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, StringHelper::formatSimple(L"TestSearchWithVariants{0}", std::to_wstring(proteinIdx)));
		FileSystem::createDirectory(outputFolder);

		st->RunTask(outputFolder, {new DbForTask(xmlName, false)}, {mzmlName}, L"");
		auto psms = File::ReadAllLines(FileSystem::combine(outputFolder, L"AllPSMs.psmtsv"));

		Assert::IsTrue(psms.Any([&] (std::any line)
		{
			line->Contains(StringHelper::formatSimple(L"\t{0}\t", variantPsmShort) + (containsVariant ? variantPsmShort : L"\t"));
		}));

		Directory::Delete(outputFolder, true);
		File::Delete(mzmlName);
		File::Delete(xmlName);
		Directory::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(Task Settings)"), true);

		std::wcout << std::wstring::Format(L"Analysis time for VariantSearchTests.SearchTests({0},{1},{2},{3}): {4}h {5}m {6}s", std::to_wstring(proteinIdx), std::to_wstring(peptideIdx), StringHelper::toString(containsVariant), variantPsmShort, stopwatch->Elapsed.Hours, stopwatch->Elapsed.Minutes, stopwatch->Elapsed.Seconds) << std::endl;

//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete mp' statement was not added since mp was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete mv' statement was not added since mv was passed to a method or constructor. Handle memory management manually.
		delete st;
		delete stopwatch;
	}

	void VariantSearchTests::MoreTests(const std::wstring &filename, DecoyType *decoyType)
	{
		std::wstring xmlName = FileSystem::combine(TestContext::CurrentContext->TestDirectory, L"TestData", filename);
		std::any un;
		auto proteins = ProteinDbLoader::LoadProteinXML(xmlName, decoyType == DecoyType::None, decoyType, nullptr, false, nullptr, un);
		auto peps = proteins[1].Digest(CommonParameters->getDigestionParams(), nullptr, nullptr).ToList();
		PeptideWithSetModifications *pep = peps[peps.size() - 2];

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		std::wstring mzmlName = StringHelper::formatSimple(L"ajgdiv{0}{1}.mzML", filename, decoyType->ToString());
		MsDataFile *myMsDataFile = new TestDataFile({pep}, true);

		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, StringHelper::formatSimple(L"TestSearchWithVariants{0}{1}", filename, decoyType->ToString()));
		FileSystem::createDirectory(outputFolder);

		SearchTask *st = new SearchTask();
		SearchParameters tempVar();
		st->setSearchParameters(&tempVar);
		st->getSearchParameters()->setDoParsimony(true);
		st->getSearchParameters()->setDecoyType(decoyType);
		st->getSearchParameters()->setSearchTarget(decoyType == DecoyType::None);
		st->getSearchParameters()->setModPeptidesAreDifferent(false);
		CommonParameters tempVar2(, , , = true, = 3, = 12, = true, = false, = 1, 1, , , , , , , , new PpmTolerance(20), , , new DigestionParams(minPeptideLength: 2));
		st->setCommonParameters(&tempVar2);

		st->RunTask(outputFolder, {new DbForTask(xmlName, false)}, {mzmlName}, L"");
		auto psms = File::ReadAllLines(FileSystem::combine(outputFolder, L"AllPSMs.psmtsv"));

		//Assert.IsTrue(psms.Any(line => line.Contains($"\t{variantPsmShort}\t" + (containsVariant ? variantPsmShort : "\t"))));

		//Directory.Delete(outputFolder, true);
		//File.Delete(mzmlName);
		//File.Delete(xmlName);
		//Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);

		//Console.WriteLine($"Analysis time for VariantSearchTests.SearchTests({proteinIdx.ToString()},{peptideIdx.ToString()},{containsVariant.ToString()},{variantPsmShort}): {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");

		delete st;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
	}
}
