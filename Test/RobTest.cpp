#include "RobTest.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/ProteinParsimony/ProteinParsimonyEngine.h"
#include "../EngineLayer/ProteinParsimony/ProteinParsimonyResults.h"
#include "../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.h"
#include "../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrResults.h"
#include "../EngineLayer/ProteinParsimony/ProteinGroup.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace Test
{

	void RobTest::TestParsimony()
	{
		// creates some proteins to test parsimony with
		std::vector<std::wstring> proteinSequences = {L"AB--------", L"--C-------", L"---D---HHH--", L"-B-D---HHH--", L"-B--E-----", L"----EFG---", L"-----F----", L"--------I-", L"-B------I-", L"----EFG--J"};

		auto proteins = std::vector<Protein*>();

		for (int i = 0; i < proteinSequences.size(); i++)
		{
			Protein tempVar(proteinSequences[i], std::to_wstring(i + 1));
			proteins.push_back(&tempVar);
		}
		Protein tempVar2(L"-----F----*", L"D1", isDecoy: true);
		proteins.push_back(&tempVar2);
		Protein tempVar3(L"-----F----**", L"C1", isContaminant: true);
		proteins.push_back(&tempVar3);
		Protein tempVar4(L"----E----**", L"C2", isContaminant: true);
		proteins.push_back(&tempVar4);

		// create the protease
		std::vector<DigestionMotif*> digestionMotifs =
		{
			new DigestionMotif(L"A", nullptr, 1, nullptr),
			new DigestionMotif(L"B", nullptr, 1, nullptr),
			new DigestionMotif(L"C", nullptr, 1, nullptr),
			new DigestionMotif(L"D", nullptr, 1, nullptr),
			new DigestionMotif(L"E", nullptr, 1, nullptr),
			new DigestionMotif(L"F", nullptr, 1, nullptr),
			new DigestionMotif(L"G", nullptr, 1, nullptr),
			new DigestionMotif(L"H", nullptr, 1, nullptr),
			new DigestionMotif(L"I", nullptr, 1, nullptr),
			new DigestionMotif(L"J", nullptr, 1, nullptr),
			new DigestionMotif(L"-", nullptr, 1, nullptr)
		};

		auto protease = new Protease(L"test", CleavageSpecificity::Full, nullptr, nullptr, digestionMotifs);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 1);

		// digest the proteins
		auto peptides = std::unordered_set<PeptideWithSetModifications*>();
		for (auto protein : proteins)
		{
			for (PeptideWithSetModifications *peptide : protein->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()))
			{
//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//				switch (peptide.BaseSequence)
//ORIGINAL LINE: case "A":
				if (peptide->BaseSequence == L"A")
				{
						peptides.insert(peptide);
				}
//ORIGINAL LINE: case "B":
				else if (peptide->BaseSequence == L"B")
				{
						peptides.insert(peptide);
				}
//ORIGINAL LINE: case "C":
				else if (peptide->BaseSequence == L"C")
				{
						peptides.insert(peptide);
				}
//ORIGINAL LINE: case "D":
				else if (peptide->BaseSequence == L"D")
				{
						peptides.insert(peptide);
				}
//ORIGINAL LINE: case "E":
				else if (peptide->BaseSequence == L"E")
				{
						peptides.insert(peptide);
				}
//ORIGINAL LINE: case "F":
				else if (peptide->BaseSequence == L"F")
				{
						peptides.insert(peptide);
				}
//ORIGINAL LINE: case "G":
				else if (peptide->BaseSequence == L"G")
				{
						peptides.insert(peptide);
				}
//ORIGINAL LINE: case "H":
				else if (peptide->BaseSequence == L"H")
				{
						peptides.insert(peptide);
				}
//ORIGINAL LINE: case "I":
				else if (peptide->BaseSequence == L"I")
				{
						peptides.insert(peptide);
				}
			}
		}

		// create PSMs for the peptides
		std::unordered_map<std::wstring, PeptideSpectralMatch*> temp;

		MzSpectrum tempVar5(new double[] {1}, new double[] {1}, false);
		MsDataScan *fakeScan = new MsDataScan(&tempVar5, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);

		CommonParameters tempVar6();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(fakeScan, 2, 0, L"File", &tempVar6);

		for (auto peptide : peptides)
		{
			TValue psm;
			std::unordered_map<std::wstring, PeptideSpectralMatch*>::const_iterator temp_iterator = temp.find(peptide.BaseSequence);
			if (temp_iterator != temp.end())
			{
				psm = temp_iterator->second;
				psm::AddOrReplace(peptide, 1, 0, true, std::vector<MatchedFragmentIon*>());
			}
			else
			{
				psm = temp_iterator->second;
				PeptideSpectralMatch tempVar7(peptide, 0, 1, 0, scan, digestionParams, new std::vector<MatchedFragmentIon*>());
				temp.emplace(peptide->BaseSequence, &tempVar7);
			}
		}

		std::vector<PeptideSpectralMatch*> psms = temp.Values->ToList();

		for (auto psm : psms)
		{
			psm->ResolveAllAmbiguities();
			psm->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
		}

		// run parsimony
		CommonParameters tempVar8();
		ProteinParsimonyEngine *parsimonyEngine = new ProteinParsimonyEngine(psms, false, &tempVar8, std::vector<std::wstring>());
		auto parsimonyResults = static_cast<ProteinParsimonyResults*>(parsimonyEngine->Run());
		auto proteinGroups = parsimonyResults->getProteinGroups();

		CommonParameters tempVar9();
		ProteinScoringAndFdrEngine *proteinScoringAndFdrEngine = new ProteinScoringAndFdrEngine(proteinGroups, psms, true, false, true, &tempVar9, std::vector<std::wstring>());
		auto proteinScoringAndFdrResults = static_cast<ProteinScoringAndFdrResults*>(proteinScoringAndFdrEngine->Run());
		proteinGroups = proteinScoringAndFdrResults->SortedAndScoredProteinGroups;

		// select the PSMs' proteins
		std::vector<std::wstring> parsimonyProteinSequences = psms.SelectMany([&] (std::any p)
		{
			p::BestMatchingPeptides->Select([&] (std::any v)
			{
				v::Peptide::Protein;
			});
		})->Select([&] (std::any v)
		{
			v::BaseSequence;
		}).Distinct().ToList();

		// check that correct proteins are in parsimony list
		Assert->Contains(L"AB--------", parsimonyProteinSequences);
		Assert->Contains(L"--C-------", parsimonyProteinSequences);
		Assert->Contains(L"-B-D---HHH--", parsimonyProteinSequences);
		Assert->Contains(L"----E----**", parsimonyProteinSequences);
		Assert->Contains(L"-B------I-", parsimonyProteinSequences);
		Assert->Contains(L"----EFG---", parsimonyProteinSequences);
		Assert->Contains(L"----EFG--J", parsimonyProteinSequences);
		Assert::AreEqual(8, parsimonyProteinSequences.size());

		// sequence coverage test
		for (auto proteinGroup : proteinGroups)
		{
			for (auto coverage : proteinGroup->getSequenceCoveragePercent())
			{
				Assert::That(coverage <= 1.0);
			}
		}

		// test protein groups
		Assert::AreEqual(3, proteinGroups.size());
		Assert::AreEqual(1, proteinGroups.front().Proteins->Count);
		Assert::AreEqual(L"AB--------", proteinGroups.front().Proteins::First().BaseSequence);
		Assert::AreEqual(2, proteinGroups.front().AllPsmsBelowOnePercentFDR->Count);
		Assert::AreEqual(2, proteinGroups.front().ProteinGroupScore);

		delete proteinScoringAndFdrEngine;
		delete parsimonyEngine;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete fakeScan' statement was not added since fakeScan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}

	void RobTest::TestPTMOutput()
	{
		std::vector<Modification*> variableModifications;
		std::vector<Modification*> fixedModifications;

		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"S", motif);
		Modification tempVar(_originalId: L"resMod", _modificationType: L"HaHa", _target: motif, _locationRestriction: L"Anywhere.", _chemicalFormula: ChemicalFormula::ParseFormula(L"H"));
		variableModifications.push_back(&tempVar);

		auto proteinList = std::vector<Protein*> {new Protein(L"MNNNSKQQQ", L"accession")};
		auto protease = new Protease(L"CustomProtease", CleavageSpecificity::Full, nullptr, nullptr, std::vector<DigestionMotif*> {new DigestionMotif(L"K", nullptr, 1, nullptr)});
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);

		std::unordered_map<Modification*, unsigned short> modsDictionary =
		{
			{variableModifications.back(), 1}
		};

		DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, maxMissedCleavages: 0, minPeptideLength: 1);

		auto modPep = proteinList.front().Digest(digestionParams, fixedModifications, variableModifications).Last();
		std::unordered_set<PeptideWithSetModifications*> value = {modPep};
		PeptideWithSetModifications *compactPeptide1 = value.First();
		Assert::AreEqual(L"QQQ", value.First().FullSequence); //this might be base

		auto firstProtDigest = proteinList.front().Digest(digestionParams, fixedModifications, variableModifications).ToList();
		std::unordered_set<PeptideWithSetModifications*> value2 = std::vector<std::unordered_set<PeptideWithSetModifications*>>(0) };
		PeptideWithSetModifications *compactPeptide2 = value2.First();
		Assert::AreEqual(L"MNNNSK", value2.First().FullSequence); //this might be base

		std::unordered_set<PeptideWithSetModifications*> value2mod = std::vector<std::unordered_set<PeptideWithSetModifications*>>(1) };
		PeptideWithSetModifications *compactPeptide2mod = value2mod.Last();
		Assert::AreEqual(L"MNNNS[HaHa:resMod on S]K", value2mod.Last().FullSequence); //this might be base

		std::unordered_set<PeptideWithSetModifications*> value3 = std::vector<std::unordered_set<PeptideWithSetModifications*>>(2) };
		PeptideWithSetModifications *compactPeptide3 = value3.First();
		Assert::AreEqual(L"NNNSK", value3.First().FullSequence); //this might be base
		std::unordered_set<PeptideWithSetModifications*> value3mod = std::vector<std::unordered_set<PeptideWithSetModifications*>>(3) };

		PeptideWithSetModifications *compactPeptide3mod = value3mod.Last();
		Assert::AreEqual(L"NNNS[HaHa:resMod on S]K", value3mod.Last().FullSequence); //this might be base

		auto peptideList = std::unordered_set<PeptideWithSetModifications*>();
		for (auto protein : proteinList)
		{
			for (auto peptide : protein->Digest(digestionParams, std::vector<Modification*>(), variableModifications))
			{
				peptideList.insert(peptide);
			}
		}

		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *jdfk = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *ms2scan = new Ms2ScanWithSpecificMass(jdfk, 2, 0, L"File", &tempVar3);

		Tolerance *fragmentTolerance = new AbsoluteTolerance(0.01);

		auto match1 = new PeptideSpectralMatch(peptideList.ElementAt(0), 0, 10, 0, ms2scan, digestionParams, std::vector<MatchedFragmentIon*>())
		{
		};
		match1->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
		auto match2 = new PeptideSpectralMatch(peptideList.ElementAt(1), 0, 10, 0, ms2scan, digestionParams, std::vector<MatchedFragmentIon*>())
		{
		};
		match2->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
		auto match3 = new PeptideSpectralMatch(peptideList.ElementAt(1), 0, 10, 0, ms2scan, digestionParams, std::vector<MatchedFragmentIon*>())
		{
		};
		match3->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);

		std::vector<PeptideSpectralMatch*> psms = {match1, match2, match3};

		std::for_each(psms.begin(), psms.end(), [&] (std::any p)
		{
			p::ResolveAllAmbiguities();
		});

		CommonParameters tempVar4();
		ProteinParsimonyEngine *engine = new ProteinParsimonyEngine(psms, true, &tempVar4, std::vector<std::wstring> {L"ff"});
		auto cool = static_cast<ProteinParsimonyResults*>(engine->Run());
		auto proteinGroups = cool->getProteinGroups();

		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *f = new ProteinScoringAndFdrEngine(proteinGroups, psms, false, false, true, &tempVar5, std::vector<std::wstring>());
		f->Run();

		Assert::AreEqual(L"#aa5[resMod on S,info:occupancy=0.67(2/3)];", proteinGroups.front().ModsInfo[0]);

		delete f;
		delete engine;
		delete match3;
		delete match2;
		delete match1;
		delete fragmentTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete ms2scan' statement was not added since ms2scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete jdfk' statement was not added since jdfk was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}

	void RobTest::TestProteinGroupsAccessionOutputOrder()
	{
		auto p = std::unordered_set<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;

		// make protein B
		Protein tempVar(L"-----F----*", L"B", nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>(), isDecoy: true);
		p.insert(&tempVar);

		// make protein A
		Protein tempVar2(L"-----F----**", L"A", nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>(), isDecoy: true);
		p.insert(&tempVar2);

		// add protein B and A to the protein group
		ProteinGroup *testGroup = new ProteinGroup(p, nullptr, nullptr);

		// test order is AB and not BA
		Assert::That(testGroup->getProteinGroupName() == L"A|B");
		Assert::That(testGroup->getProteins().First().Accession->Equals(L"B"));

		delete testGroup;
	}
}
