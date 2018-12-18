#include "StefanParsimonyTest.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/ProteinParsimony/ProteinParsimonyEngine.h"
#include "../EngineLayer/ProteinParsimony/ProteinParsimonyResults.h"
#include "../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.h"
#include "../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrResults.h"

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace Test
{

	void StefanParsimonyTest::ParsimonyTreatModifiedFormsAsUnique()
	{
		bool modPeptidesAreUnique = true;

		// set up mods
		auto modDictionary = std::unordered_map<int, std::vector<Modification*>>();
		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"M", motif1);
		auto mod = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 15.99491461957);

		// modified version of protein
		auto protein1 = new Protein(L"PEPTIDEM", L"accession1");

		// unmodified version of protein
		auto protein2 = new Protein(L"YYYKPEPTIDEM", L"accession2");

		DigestionParams tempVar(protease: L"trypsin", minPeptideLength: 1);
		std::vector<PeptideWithSetModifications*> pwsmsFromProtein1 = protein1->Digest(&tempVar, {mod}, std::vector<Modification*>()).ToList(); //this is a fixed mod
		DigestionParams tempVar2(protease: L"trypsin", minPeptideLength: 1);
		std::vector<PeptideWithSetModifications*> pwsmsFromProtein2 = protein2->Digest(&tempVar2, std::vector<Modification*>(), std::vector<Modification*>()).ToList();

		// check to make sure mod is present
		PeptideWithSetModifications *modifiedPeptide = pwsmsFromProtein1[0];
		PeptideWithSetModifications *unmodifiedPeptide = pwsmsFromProtein2[1];

		Assert::That(!modifiedPeptide->FullSequence->Equals(unmodifiedPeptide->FullSequence)); // sequences should not be equal (one has a mod)
		Assert::That(modifiedPeptide->BaseSequence->Equals(unmodifiedPeptide->BaseSequence)); // base sequences should be equal
		Assert::That(modifiedPeptide->NumMods == 1); // methionine was oxidized on this protein
		Assert::That(unmodifiedPeptide->NumMods == 0); // there was no modification on this protein

		// build PSMs for parsimony
		std::vector<PeptideSpectralMatch*> psmsForParsimony;

		MzSpectrum tempVar3(new double[] {1}, new double[] {1}, false);
		MsDataScan *fakeScan = new MsDataScan(&tempVar3, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);

		CommonParameters tempVar4();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(fakeScan, 2, 0, L"File", &tempVar4);

		DigestionParams tempVar5();
		PeptideSpectralMatch *psm1 = new PeptideSpectralMatch(modifiedPeptide, 0, 10, 1, scan, &tempVar5, std::vector<MatchedFragmentIon*>());
		psm1->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
		psm1->ResolveAllAmbiguities();

		DigestionParams tempVar6();
		PeptideSpectralMatch *psm2 = new PeptideSpectralMatch(unmodifiedPeptide, 0, 10, 2, scan, &tempVar6, std::vector<MatchedFragmentIon*>());
		psm2->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
		psm2->ResolveAllAmbiguities();

		psmsForParsimony.push_back(psm1);
		psmsForParsimony.push_back(psm2);

		// apply parsimony
		CommonParameters tempVar7();
		ProteinParsimonyEngine *pae = new ProteinParsimonyEngine(psmsForParsimony, modPeptidesAreUnique, &tempVar7, std::vector<std::wstring>());

		// because the two chosen peptides are the same, we should end up with both protein accessions still in the list
		auto proteinParsimonyResult = static_cast<ProteinParsimonyResults*>(pae->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar8();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinParsimonyResult->getProteinGroups(), psmsForParsimony, false, true, true, &tempVar8, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		int countOfProteinGroups = results->SortedAndScoredProteinGroups.size();

		// because modified peptides were considered as unique, then there should be two protein groups after parsimony, and one protein accession for each peptide
		Assert::That(countOfProteinGroups == 2);
		Assert::That(results->SortedAndScoredProteinGroups.All([&] (std::any p)
		{
		delete proteinScoringEngine;
		delete pae;
//C# TO C++ CONVERTER TODO TASK: A 'delete psm2' statement was not added since psm2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete psm1' statement was not added since psm1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete fakeScan' statement was not added since fakeScan was passed to a method or constructor. Handle memory management manually.
		delete protein2;
		delete protein1;
//C# TO C++ CONVERTER TODO TASK: A 'delete mod' statement was not added since mod was passed to a method or constructor. Handle memory management manually.
			return p::Proteins->Count == 1;
		}));
		Assert::That(psm1->getProteinAccession() == L"accession1");
		Assert::That(psm2->getProteinAccession() == L"accession2");

		delete proteinScoringEngine;
		delete pae;
//C# TO C++ CONVERTER TODO TASK: A 'delete psm2' statement was not added since psm2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete psm1' statement was not added since psm1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete fakeScan' statement was not added since fakeScan was passed to a method or constructor. Handle memory management manually.
		delete protein2;
		delete protein1;
//C# TO C++ CONVERTER TODO TASK: A 'delete mod' statement was not added since mod was passed to a method or constructor. Handle memory management manually.
	}

	void StefanParsimonyTest::ParsimonyDontTreatModifiedFormsAsUnique()
	{
		bool modPeptidesAreUnique = false;

		// set up mods
		auto modDictionary = std::unordered_map<int, std::vector<Modification*>>();
		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"M", motif1);
		auto mod = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 15.99491461957);

		// modified version of protein
		auto protein1 = new Protein(L"PEPTIDEM", L"accession1");

		// unmodified version of protein
		auto protein2 = new Protein(L"YYYKPEPTIDEM", L"accession2");

		DigestionParams tempVar(protease: L"trypsin", minPeptideLength: 1);
		std::vector<PeptideWithSetModifications*> pwsmsFromProtein1 = protein1->Digest(&tempVar, {mod}, std::vector<Modification*>()).ToList(); //this is a fixed mod
		DigestionParams tempVar2(protease: L"trypsin", minPeptideLength: 1);
		std::vector<PeptideWithSetModifications*> pwsmsFromProtein2 = protein2->Digest(&tempVar2, std::vector<Modification*>(), std::vector<Modification*>()).ToList();

		// check to make sure mod is present
		PeptideWithSetModifications *modifiedPeptide = pwsmsFromProtein1[0];
		PeptideWithSetModifications *unmodifiedPeptide = pwsmsFromProtein2[1];

		Assert::That(!modifiedPeptide->FullSequence->Equals(unmodifiedPeptide->FullSequence)); // sequences should not be equal (one has a mod)
		Assert::That(modifiedPeptide->BaseSequence->Equals(unmodifiedPeptide->BaseSequence)); // base sequences should be equal
		Assert::That(modifiedPeptide->NumMods == 1); // methionine was oxidized on this protein
		Assert::That(unmodifiedPeptide->NumMods == 0); // there was no modification on this protein

		// build PSMs for parsimony
		std::vector<PeptideSpectralMatch*> psmsForParsimony;

		MzSpectrum tempVar3(new double[] {1}, new double[] {1}, false);
		MsDataScan *fakeScan = new MsDataScan(&tempVar3, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);

		CommonParameters tempVar4();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(fakeScan, 2, 0, L"File", &tempVar4);

		DigestionParams tempVar5();
		PeptideSpectralMatch *psm1 = new PeptideSpectralMatch(modifiedPeptide, 0, 10, 1, scan, &tempVar5, std::vector<MatchedFragmentIon*>());
		psm1->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
		psm1->ResolveAllAmbiguities();

		DigestionParams tempVar6();
		PeptideSpectralMatch *psm2 = new PeptideSpectralMatch(unmodifiedPeptide, 0, 10, 2, scan, &tempVar6, std::vector<MatchedFragmentIon*>());
		psm2->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
		psm2->ResolveAllAmbiguities();

		psmsForParsimony.push_back(psm1);
		psmsForParsimony.push_back(psm2);

		// apply parsimony
		CommonParameters tempVar7();
		ProteinParsimonyEngine *pae = new ProteinParsimonyEngine(psmsForParsimony, modPeptidesAreUnique, &tempVar7, std::vector<std::wstring>());

		// because the two chosen peptides are the same, we should end up with both protein accessions still in the list
		auto proteinParsimonyResult = static_cast<ProteinParsimonyResults*>(pae->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar8();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinParsimonyResult->getProteinGroups(), psmsForParsimony, false, true, true, &tempVar8, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		int countOfProteinGroups = results->SortedAndScoredProteinGroups.size();

		// because modified peptides were NOT considered as unique, 
		// then there should be one ambiguous protein group after parsimony, 
		// and two protein accessions for each peptide
		Assert::AreEqual(1, countOfProteinGroups);
		Assert::AreEqual(2, results->SortedAndScoredProteinGroups.front().Proteins->Count);
		Assert::IsNull(psm1->getProteinAccession());
		Assert::IsNull(psm2->getProteinAccession());

		delete proteinScoringEngine;
		delete pae;
//C# TO C++ CONVERTER TODO TASK: A 'delete psm2' statement was not added since psm2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete psm1' statement was not added since psm1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete fakeScan' statement was not added since fakeScan was passed to a method or constructor. Handle memory management manually.
		delete protein2;
		delete protein1;
//C# TO C++ CONVERTER TODO TASK: A 'delete mod' statement was not added since mod was passed to a method or constructor. Handle memory management manually.
	}

	void StefanParsimonyTest::ParsimonyWeirdCatch()
	{
		Protein *protein1 = new Protein(L"MATSIK", L"protein1", isDecoy: true);
		Protein *protein2 = new Protein(L"MATSIK", L"protein2");

		std::vector<Modification*> allKnownFixedModifications;
		DigestionParams *digestionParams = new DigestionParams(minPeptideLength: 5);
		std::vector<Modification*> variableModifications;
		auto pep1 = protein1->Digest(digestionParams, allKnownFixedModifications, variableModifications).First();
		auto pep2 = protein2->Digest(digestionParams, allKnownFixedModifications, variableModifications).First();

		// build the dictionary for input to parsimony
		MzSpectrum tempVar(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar2();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar2);

		DigestionParams tempVar3();
		std::vector<PeptideSpectralMatch*> psms = {new PeptideSpectralMatch(pep1,0,1,0, scan, &tempVar3, std::vector<MatchedFragmentIon*>())};

		// this PSM has a target and a decoy
		psms[0]->AddOrReplace(pep2, 1, 0, true, nullptr);

		std::for_each(psms.begin(), psms.end(), [&] (std::any p)
		{
			p::ResolveAllAmbiguities();
		});
		std::for_each(psms.begin(), psms.end(), [&] (std::any p)
		{
			p::SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
		});

		// apply parsimony
		CommonParameters tempVar4();
		ProteinParsimonyEngine *pae = new ProteinParsimonyEngine(psms, false, &tempVar4, std::vector<std::wstring>());

		// because the two chosen peptides are the same, we should end up with both protein accessions still in the list
		auto proteinParsimonyResult = static_cast<ProteinParsimonyResults*>(pae->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinParsimonyResult->getProteinGroups(), psms, false, true, true, &tempVar5, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		int countOfProteinGroups = results->SortedAndScoredProteinGroups.size();

		// only target protein gets generated
		Assert::That(countOfProteinGroups == 1);
		Assert::That(results->SortedAndScoredProteinGroups.front().Proteins->Count == 1);
		Assert::That(!results->SortedAndScoredProteinGroups.front().IsDecoy);

		delete proteinScoringEngine;
		delete pae;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
		delete protein2;
		delete protein1;
	}
}
