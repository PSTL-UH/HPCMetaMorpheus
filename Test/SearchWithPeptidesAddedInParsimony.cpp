#include "SearchWithPeptidesAddedInParsimony.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../EngineLayer/CommonParameters.h"
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

	void SearchWithPeptidesAddedInParsimony::SearchWithPeptidesAddedInParsimonyTest()
	{
		// Make sure can run the complete search task when multiple compact peptides may correspond to a single PWSM
		SearchTask *st = new SearchTask();
		SearchParameters tempVar();
		st->setSearchParameters(&tempVar);
		st->getSearchParameters()->setDoParsimony(true);
		st->getSearchParameters()->setDecoyType(DecoyType::None);
		st->getSearchParameters()->setModPeptidesAreDifferent(false);
		CommonParameters tempVar2(, , = true, = true, = 3, = 12, = true, = false, = 1, 1, , , , , , , , , , , new DigestionParams(minPeptideLength: 2));
		st->setCommonParameters(&tempVar2);

		std::wstring xmlName = L"andguiaheow.xml";

		DigestionParams tempVar3(maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: 2, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain, maxModsForPeptides: 1);
		CommonParameters *CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: &tempVar3);

		ModificationMotif motifA;
		ModificationMotif::TryGetMotif(L"A", motifA);
		Modification *alanineMod = new Modification(_originalId: L"111", _modificationType: L"mt", _target: motifA, _locationRestriction: L"Anywhere.", _monoisotopicMass: 111);

		auto variableModifications = std::vector<Modification*>();
		std::unordered_map<int, std::vector<Modification*>> oneBasedModifications1 =
		{
			{
				2, {alanineMod}
			}
		};
		Protein *protein1 = new Protein(L"MA", L"protein1", oneBasedModifications: oneBasedModifications1);
		// Alanine = Glycine + CH2

		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"G", motif1);

		Modification *glycineMod = new Modification(_originalId: L"CH2 on Glycine", _modificationType: L"mt", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: Chemistry::ChemicalFormula::ParseFormula(L"CH2").MonoisotopicMass);

		std::unordered_map<int, std::vector<Modification*>> oneBasedModifications2 =
		{
			{
				2, {glycineMod}
			}
		};
		Protein *protein2 = new Protein(L"MG", L"protein3", oneBasedModifications: oneBasedModifications2);

		PeptideWithSetModifications *pepMA = protein1->Digest(CommonParameters->getDigestionParams(), std::vector<Modification*>(), variableModifications).First();
		PeptideWithSetModifications *pepMA111 = protein1->Digest(CommonParameters->getDigestionParams(), std::vector<Modification*>(), variableModifications).Last();

		auto pepMG = protein2->Digest(CommonParameters->getDigestionParams(), std::vector<Modification*>(), variableModifications).First();

		ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>(), {protein1, protein2}, xmlName);

		std::wstring mzmlName = LR"(ajgdiu.mzML)";

		MsDataFile *myMsDataFile = new TestDataFile({pepMA, pepMG, pepMA111}, true);

		IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestSearchWithPeptidesAddedInParsimony)");
		FileSystem::createDirectory(outputFolder);

		st->RunTask(outputFolder, {new DbForTask(xmlName, false)}, {mzmlName}, L"");
		Directory::Delete(outputFolder, true);
		File::Delete(mzmlName);
		File::Delete(xmlName);
		Directory::Delete(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(Task Settings)"), true);

//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protein2' statement was not added since protein2 was passed to a method or constructor. Handle memory management manually.
		delete glycineMod;
//C# TO C++ CONVERTER TODO TASK: A 'delete protein1' statement was not added since protein1 was passed to a method or constructor. Handle memory management manually.
		delete alanineMod;
		delete CommonParameters;
		delete st;
	}
}
