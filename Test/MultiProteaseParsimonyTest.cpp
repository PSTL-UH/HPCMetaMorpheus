#include "MultiProteaseParsimonyTest.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/ProteinParsimony/ProteinParsimonyEngine.h"
#include "../EngineLayer/ProteinParsimony/ProteinParsimonyResults.h"
#include "../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.h"
#include "../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrResults.h"
#include "../EngineLayer/ProteinParsimony/ProteinGroup.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "SetUpTests.h"

using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;

namespace Test
{

	void MultiProteaseParsimonyTest::MultiProteaseTest()
	{
		std::vector<std::wstring> sequences = {L"ABCKPEPR", L"BRPEPR", L"ARPEPR"};

		auto p = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			p.push_back(&tempVar);
		}

		DigestionParams *digestionParams_Tryp = new DigestionParams(protease: L"trypsin", minPeptideLength: 1);
		DigestionParams *digestionParams_ArgC = new DigestionParams(protease: L"Arg-C", minPeptideLength: 1);

		PeptideWithSetModifications *pepA_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams_Tryp, oneBasedStartResidueInProtein: 5, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"PEPR", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_2T = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams_Tryp, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"PEPR", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_3T = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: digestionParams_Tryp, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"PEPR", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepB_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams_Tryp, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABCK", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_2A = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams_ArgC, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"PEPR", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_3A = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: digestionParams_ArgC, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"PEPR", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);

		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		PeptideSpectralMatch *psmPEPR_T = new PeptideSpectralMatch(pepA_1T, 0, 10, 0, scan, digestionParams_Tryp, std::vector<MatchedFragmentIon*>());
		psmPEPR_T->AddOrReplace(pepA_2T, 10, 0, true, std::vector<MatchedFragmentIon*>());
		psmPEPR_T->AddOrReplace(pepA_3T, 10, 0, true, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmPEPR_A = new PeptideSpectralMatch(pepA_2A, 0, 10, 0, scan, digestionParams_ArgC, std::vector<MatchedFragmentIon*>());
		psmPEPR_A->AddOrReplace(pepA_3A, 10, 0, true, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmABCK_T = new PeptideSpectralMatch(pepB_1T, 0, 10, 0, scan, digestionParams_Tryp, std::vector<MatchedFragmentIon*>());

		std::vector<PeptideSpectralMatch*> psms = {psmPEPR_T, psmPEPR_A, psmABCK_T};
		std::for_each(psms.begin(), psms.end(), [&] (std::any j)
		{
			j::ResolveAllAmbiguities();
		});
		std::for_each(psms.begin(), psms.end(), [&] (std::any j)
		{
			j::SetFdrValues(1, 0, 0, 1, 0, 0, NAN, NAN, NAN, false);
		});

		std::unordered_set<DigestionParams*> digestionParamsList;
		digestionParamsList.insert(digestionParams_Tryp);
		digestionParamsList.insert(digestionParams_ArgC);
		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"M", motif1);
		Modification *mod = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modVarList = {mod};

		ModificationMotif motif2;
		ModificationMotif::TryGetMotif(L"M", motif2);
		Modification *mod2 = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif2, _locationRestriction: L"Anyhwere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modFixedList = {mod};

		CommonParameters tempVar4();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar4, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar5, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;
		Assert::AreEqual(2, proteinGroups.size());
		auto proteinGroup1 = proteinGroups.Where([&] (std::any h)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmABCK_T;
		delete psmPEPR_A;
		delete psmPEPR_T;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_3A' statement was not added since pepA_3A was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2A' statement was not added since pepA_2A was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1T' statement was not added since pepB_1T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_3T' statement was not added since pepA_3T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2T' statement was not added since pepA_2T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1T' statement was not added since pepA_1T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams_ArgC' statement was not added since digestionParams_ArgC was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams_Tryp' statement was not added since digestionParams_Tryp was passed to a method or constructor. Handle memory management manually.
			return h->ProteinGroupName == L"1";
		}).First();
		Assert::AreEqual(1, proteinGroup1->UniquePeptides->Count);
		Assert::AreEqual(2, proteinGroup1->AllPeptides->Count);
		auto proteinGroup2 = proteinGroups.Where([&] (std::any h)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmABCK_T;
		delete psmPEPR_A;
		delete psmPEPR_T;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_3A' statement was not added since pepA_3A was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2A' statement was not added since pepA_2A was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1T' statement was not added since pepB_1T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_3T' statement was not added since pepA_3T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2T' statement was not added since pepA_2T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1T' statement was not added since pepA_1T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams_ArgC' statement was not added since digestionParams_ArgC was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams_Tryp' statement was not added since digestionParams_Tryp was passed to a method or constructor. Handle memory management manually.
			return h->ProteinGroupName == L"2|3";
		}).First();
		Assert::AreEqual(0, proteinGroup2->UniquePeptides->Count);
		Assert::AreEqual(4, proteinGroup2->AllPeptides->Count);

		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmABCK_T;
		delete psmPEPR_A;
		delete psmPEPR_T;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_3A' statement was not added since pepA_3A was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2A' statement was not added since pepA_2A was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1T' statement was not added since pepB_1T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_3T' statement was not added since pepA_3T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2T' statement was not added since pepA_2T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1T' statement was not added since pepA_1T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams_ArgC' statement was not added since digestionParams_ArgC was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams_Tryp' statement was not added since digestionParams_Tryp was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseSamePeptideSameProteinsDifferentProteases()
	{
		std::vector<std::wstring> sequences = {L"-XYZ--ABC", L"-XYZ-EFG-ABC"};
		//both proteases are cleaving at the same spots to simulate trypsin and argC producing the same peptides

		std::vector<DigestionMotif*> motifs =
		{
			new DigestionMotif(L"-", nullptr, 1, nullptr),
			new DigestionMotif(L"Z", nullptr, 1, nullptr)
		};
		auto protease1 = new Protease(L"proteaseA1", CleavageSpecificity::Full, nullptr, nullptr, motifs);
		ProteaseDictionary::Dictionary->Add(protease1->Name, protease1);
		auto protease2 = new Protease(L"proteaseB1", CleavageSpecificity::Full, nullptr, nullptr, motifs);
		ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);

		//a hashset of peptideWithSetModifications from all proteases
		auto pwsmList = std::unordered_set<PeptideWithSetModifications*>();

		auto proteinList = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			proteinList.push_back(&tempVar);
		}

		DigestionParams *digestionParams = new DigestionParams(protease: protease1->Name, minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: protease2->Name, minPeptideLength: 1);

		PeptideWithSetModifications *pepA_1Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 12, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepB_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"EFG", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_1Dp2 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_2Dp2 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 12, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);

		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		PeptideSpectralMatch *psmABC_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		psmABC_Dp1->AddOrReplace(pepA_2Dp1, 10, 0, true, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmABC_Dp2 = new PeptideSpectralMatch(pepA_1Dp2, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		psmABC_Dp2->AddOrReplace(pepA_2Dp2, 10, 0, true, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmEFG_Dp1 = new PeptideSpectralMatch(pepB_2Dp1, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());

		// builds psm list to match to peptides
		std::vector<PeptideSpectralMatch*> psms = {psmABC_Dp1, psmABC_Dp2, psmEFG_Dp1};

		std::for_each(psms.begin(), psms.end(), [&] (std::any p)
		{
			p::ResolveAllAmbiguities();
		});
		std::for_each(psms.begin(), psms.end(), [&] (std::any p)
		{
			p::SetFdrValues(1, 0, 0, 1, 0, 0, NAN, NAN, NAN, false);
		});

		std::unordered_set<DigestionParams*> digestionParamsList;
		digestionParamsList.insert(digestionParams);
		digestionParamsList.insert(digestionParams2);
		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"M", motif1);
		Modification *mod = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modVarList = {mod};

		ModificationMotif motif2;
		ModificationMotif::TryGetMotif(L"M", motif2);
		Modification *mod2 = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif2, _locationRestriction: L"Anyhwere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modFixedList = {mod};

		CommonParameters tempVar4();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar4, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar5, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;
		// should result in 1 protein group (protein2)
		Assert::AreEqual(1, proteinGroups.size());
		Assert::AreEqual(L"2", proteinGroups.ElementAt(0).ProteinGroupName);

		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmEFG_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp2' statement was not added since pepA_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp2' statement was not added since pepA_1Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp1' statement was not added since pepB_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp1' statement was not added since pepA_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease1' statement was not added since protease1 was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseParsimony_SharedSequenceCanBeUniquePeptide()
	{
		std::vector<std::wstring> sequences = {L"-XYZ--ABC", L"-XYZ-EFGABC"};

		std::vector<DigestionMotif*> motifs1 =
		{
			new DigestionMotif(L"-", nullptr, 1, nullptr),
			new DigestionMotif(L"G", nullptr, 1, nullptr)
		};
		std::vector<DigestionMotif*> motifs2 = {new DigestionMotif(L"G", nullptr, 1, nullptr)};

		auto protease1 = new Protease(L"proteaseA", CleavageSpecificity::Full, nullptr, nullptr, motifs1);
		ProteaseDictionary::Dictionary->Add(protease1->Name, protease1);
		auto protease2 = new Protease(L"proteaseB", CleavageSpecificity::Full, nullptr, nullptr, motifs2);
		ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);

		//a hashset of peptideWithSetModifications from all proteases
		auto pwsmList = std::unordered_set<PeptideWithSetModifications*>();

		auto proteinList = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			proteinList.push_back(&tempVar);
		}

		DigestionParams *digestionParams = new DigestionParams(protease: protease1->Name, minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: protease2->Name, minPeptideLength: 1);

		PeptideWithSetModifications *pepA_1Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"XYZ", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"XYZ", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepB_1Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepC_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"EFGABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepB_2Dp2 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);

		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		PeptideSpectralMatch *psmABC_Dp1 = new PeptideSpectralMatch(pepB_1Dp1, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmABC_Dp2 = new PeptideSpectralMatch(pepB_2Dp2, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmEFGABC_Dp1 = new PeptideSpectralMatch(pepC_2Dp1, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmXYZ_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		psmXYZ_Dp1->AddOrReplace(pepA_2Dp1, 10, 0, true, std::vector<MatchedFragmentIon*>());

		// builds psm list to match to peptides
		std::vector<PeptideSpectralMatch*> psms = {psmABC_Dp1, psmABC_Dp2, psmEFGABC_Dp1, psmXYZ_Dp1};

		std::for_each(psms.begin(), psms.end(), [&] (std::any p)
		{
			p::ResolveAllAmbiguities();
		});
		std::for_each(psms.begin(), psms.end(), [&] (std::any p)
		{
			p::SetFdrValues(1, 0, 0, 1, 0, 0, NAN, NAN, NAN, false);
		});

		std::unordered_set<DigestionParams*> digestionParamsList;
		digestionParamsList.insert(digestionParams);
		digestionParamsList.insert(digestionParams2);
		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"M", motif1);
		Modification *mod = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modVarList = {mod};

		ModificationMotif motif2;
		ModificationMotif::TryGetMotif(L"M", motif2);
		Modification *mod2 = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif2, _locationRestriction: L"Anyhwere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modFixedList = {mod};

		CommonParameters tempVar4();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar4, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar5, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;
		Assert::AreEqual(2, proteinGroups.size());

		auto proteinGroup1 = proteinGroups.Where([&] (std::any p)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmXYZ_Dp1;
		delete psmEFGABC_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepC_2Dp1' statement was not added since pepC_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp1' statement was not added since pepA_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease1' statement was not added since protease1 was passed to a method or constructor. Handle memory management manually.
			return p->ProteinGroupName == L"1";
		}).First();
		Assert::AreEqual(2, proteinGroup1->AllPeptides->Count);
		Assert::AreEqual(1, proteinGroup1->UniquePeptides->Count);
		auto pg1pep1 = proteinGroup1->AllPeptides.Where([&] (std::any p)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmXYZ_Dp1;
		delete psmEFGABC_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepC_2Dp1' statement was not added since pepC_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp1' statement was not added since pepA_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease1' statement was not added since protease1 was passed to a method or constructor. Handle memory management manually.
			return p->BaseSequence == L"XYZ";
		}).First();
		Assert::That(pg1pep1->DigestionParams.Protease->Name == L"proteaseA");
		auto pg1pep2 = proteinGroup1->AllPeptides.Where([&] (std::any p)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmXYZ_Dp1;
		delete psmEFGABC_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepC_2Dp1' statement was not added since pepC_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp1' statement was not added since pepA_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease1' statement was not added since protease1 was passed to a method or constructor. Handle memory management manually.
			return p->BaseSequence == L"ABC";
		}).First();
		Assert::That(pg1pep2->DigestionParams.Protease->Name == L"proteaseA");
		Assert::That(proteinGroup1->UniquePeptides.First().BaseSequence->Equals(L"ABC"));

		auto proteinGroup2 = proteinGroups.Where([&] (std::any p)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmXYZ_Dp1;
		delete psmEFGABC_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepC_2Dp1' statement was not added since pepC_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp1' statement was not added since pepA_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease1' statement was not added since protease1 was passed to a method or constructor. Handle memory management manually.
			return p->ProteinGroupName == L"2";
		}).First();
		Assert::AreEqual(3, proteinGroup2->AllPeptides->Count);
		Assert::AreEqual(2, proteinGroup2->UniquePeptides->Count);
		auto pg2pep1 = proteinGroup2->AllPeptides.Where([&] (std::any p)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmXYZ_Dp1;
		delete psmEFGABC_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepC_2Dp1' statement was not added since pepC_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp1' statement was not added since pepA_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease1' statement was not added since protease1 was passed to a method or constructor. Handle memory management manually.
			return p->BaseSequence == L"XYZ";
		}).First();
		Assert::That(pg2pep1->DigestionParams.Protease->Name == L"proteaseA");
		auto pg2pep2 = proteinGroup2->AllPeptides.Where([&] (std::any p)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmXYZ_Dp1;
		delete psmEFGABC_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepC_2Dp1' statement was not added since pepC_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp1' statement was not added since pepA_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease1' statement was not added since protease1 was passed to a method or constructor. Handle memory management manually.
			return p->BaseSequence == L"ABC";
		}).First();
		Assert::That(pg2pep2->DigestionParams.Protease->Name == L"proteaseB");
		auto pg2pep3 = proteinGroup2->AllPeptides.Where([&] (std::any p)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmXYZ_Dp1;
		delete psmEFGABC_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepC_2Dp1' statement was not added since pepC_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp1' statement was not added since pepA_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease1' statement was not added since protease1 was passed to a method or constructor. Handle memory management manually.
			return p->BaseSequence == L"EFGABC";
		}).First();
		Assert::That(pg2pep3->DigestionParams.Protease->Name == L"proteaseA");
		auto uniquePeptideSequences = proteinGroup2->UniquePeptides->Select([&] (std::any p)
		{
			p::BaseSequence;
		}).ToList();
		Assert::That(std::find(uniquePeptideSequences.begin(), uniquePeptideSequences.end(), L"ABC") != uniquePeptideSequences.end());
		Assert::That(std::find(uniquePeptideSequences.begin(), uniquePeptideSequences.end(), L"EFGABC") != uniquePeptideSequences.end());

		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmXYZ_Dp1;
		delete psmEFGABC_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepC_2Dp1' statement was not added since pepC_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp1' statement was not added since pepA_2Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease1' statement was not added since protease1 was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseParsimony_IndistringuishableProteinsNowDistinguishable()
	{
		std::vector<std::wstring> sequences = {L"ABCEFG", L"EFGABC"};

		std::vector<std::tuple<std::wstring, FragmentationTerminus*>> sequencesInducingCleavage = {std::tuple<std::wstring, FragmentationTerminus*>(L"C", FragmentationTerminus::C)};
		std::vector<std::tuple<std::wstring, FragmentationTerminus*>> sequencesInducingCleavage2 = {std::tuple<std::wstring, FragmentationTerminus*>(L"G", FragmentationTerminus::C)};

		std::vector<DigestionMotif*> motifs1 = {new DigestionMotif(L"C", nullptr, 1, nullptr)};
		std::vector<DigestionMotif*> motifs2 = {new DigestionMotif(L"G", nullptr, 1, nullptr)};

		auto protease = new Protease(L"testA", CleavageSpecificity::Full, nullptr, nullptr, motifs1);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		auto protease2 = new Protease(L"testB", CleavageSpecificity::Full, nullptr, nullptr, motifs2);
		ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);
		auto peptideList = std::unordered_set<PeptideWithSetModifications*>();

		auto p = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			p.push_back(&tempVar);
		}

		DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: protease2->Name, minPeptideLength: 1);

		PeptideWithSetModifications *pepA_1Dp1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_2Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 4, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepB_1Dp1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 4, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"EFG", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepB_2Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"EFG", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);

		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		PeptideSpectralMatch *psmABC_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmABC_Dp2 = new PeptideSpectralMatch(pepA_2Dp2, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmEFG_Dp1 = new PeptideSpectralMatch(pepB_1Dp1, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmEFG_Dp2 = new PeptideSpectralMatch(pepB_2Dp2, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());

		// builds psm list to match to peptides
		std::vector<PeptideSpectralMatch*> psms = {psmABC_Dp1, psmABC_Dp2, psmEFG_Dp1, psmEFG_Dp2};

		std::for_each(psms.begin(), psms.end(), [&] (std::any h)
		{
			h::ResolveAllAmbiguities();
		});
		std::for_each(psms.begin(), psms.end(), [&] (std::any h)
		{
			h::SetFdrValues(1, 0, 0, 1, 0, 0, NAN, NAN, NAN, false);
		});

		CommonParameters tempVar4();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar4, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar5, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;

		Assert::AreEqual(2, proteinGroups.size());

		// check first protein group
		ProteinGroup *pg1 = proteinGroups.Where([&] (std::any v)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete psmEFG_Dp2;
		delete psmEFG_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp2' statement was not added since pepA_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
			return v->ProteinGroupName == L"1";
		}).First();
		PeptideWithSetModifications *pg1pep1 = pg1->getAllPeptides().Where([&] (std::any v)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete psmEFG_Dp2;
		delete psmEFG_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp2' statement was not added since pepA_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
			return v->BaseSequence == L"ABC";
		}).First();
		PeptideWithSetModifications *pg1pep2 = pg1->getAllPeptides().Where([&] (std::any v)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete psmEFG_Dp2;
		delete psmEFG_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp2' statement was not added since pepA_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
			return v->BaseSequence == L"EFG";
		}).First();
		Assert::That(std::find(pg1->UniquePeptides.begin(), pg1->UniquePeptides.end(), pg1pep1) != pg1->UniquePeptides.end());
		Assert::That(pg1pep1->DigestionParams.Protease->Name == L"testA");
		Assert::That(std::find(pg1->UniquePeptides.begin(), pg1->UniquePeptides.end(), pg1pep2) != pg1->UniquePeptides.end());
		Assert::That(pg1pep2->DigestionParams.Protease->Name == L"testA");
		Assert::That(pg1->getAllPeptides().size() == 2);
		Assert::That(pg1->getUniquePeptides().size() == 2);

		// check second protein group
		ProteinGroup *pg2 = proteinGroups.Where([&] (std::any v)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete psmEFG_Dp2;
		delete psmEFG_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp2' statement was not added since pepA_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
			return v->ProteinGroupName == L"2";
		}).First();
		PeptideWithSetModifications *pg2pep1 = pg2->getAllPeptides().Where([&] (std::any v)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete psmEFG_Dp2;
		delete psmEFG_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp2' statement was not added since pepA_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
			return v->BaseSequence == L"ABC";
		}).First();
		PeptideWithSetModifications *pg2pep2 = pg2->getAllPeptides().Where([&] (std::any v)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete psmEFG_Dp2;
		delete psmEFG_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp2' statement was not added since pepA_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
			return v->BaseSequence == L"EFG";
		}).First();
		Assert::That(std::find(pg2->UniquePeptides.begin(), pg2->UniquePeptides.end(), pg2pep1) != pg2->UniquePeptides.end());
		Assert::That(pg2pep1->DigestionParams.Protease->Name == L"testB");
		Assert::That(std::find(pg2->UniquePeptides.begin(), pg2->UniquePeptides.end(), pg2pep2) != pg2->UniquePeptides.end());
		Assert::That(pg2pep2->DigestionParams.Protease->Name == L"testB");
		Assert::That(pg2->getAllPeptides().size() == 2);
		Assert::That(pg2->getUniquePeptides().size() == 2);

		delete proteinScoringEngine;
		delete ppe;
		delete psmEFG_Dp2;
		delete psmEFG_Dp1;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_2Dp2' statement was not added since pepB_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepB_1Dp1' statement was not added since pepB_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp2' statement was not added since pepA_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseParsimony_SameAminoAcidsResultInTwoUniquePeptidesForOneProtein()
	{
		std::vector<std::wstring> sequences = {L"-XYZ-EFGABC"};

		std::vector<DigestionMotif*> motifs1 =
		{
			new DigestionMotif(L"A", nullptr, 0, nullptr),
			new DigestionMotif(L"Z", nullptr, 1, nullptr)
		};
		std::vector<DigestionMotif*> motifs2 = {new DigestionMotif(L"G", nullptr, 1, nullptr)};

		auto protease = new Protease(L"test1", CleavageSpecificity::Full, nullptr, nullptr, motifs1);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		auto protease2 = new Protease(L"test2", CleavageSpecificity::Full, nullptr, nullptr, motifs2);
		ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);

		auto p = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			p.push_back(&tempVar);
		}

		DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: protease2->Name, minPeptideLength: 1);

		PeptideWithSetModifications *pepA = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepB = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);

		auto peptideList = std::unordered_set<PeptideWithSetModifications*> {pepA, pepB};

		// builds psm list to match to peptides
		std::vector<PeptideSpectralMatch*> psms;
		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		// creates psms for specific PWSM
		for (auto pwsm : peptideList)
		{
			if (pwsm->DigestionParams == digestionParams)
			{
				PeptideSpectralMatch tempVar4(pwsm, 0, 10, 0, scan, digestionParams, new std::vector<MatchedFragmentIon*>());
				psms.push_back(&tempVar4);
			}
			if (pwsm->DigestionParams == digestionParams2)
			{
				PeptideSpectralMatch tempVar5(pwsm, 0, 10, 0, scan, digestionParams2, new std::vector<MatchedFragmentIon*>());
				psms.push_back(&tempVar5);
			}
			psms.back().SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
			psms.back().ResolveAllAmbiguities();
		}

		CommonParameters tempVar6();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar6, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar7();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar7, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;

		Assert::AreEqual(1, proteinGroups.size());
		Assert::AreEqual(L"1", proteinGroups.ElementAt(0).ProteinGroupName);
		Assert::AreEqual(2, proteinGroups.ElementAt(0).UniquePeptides->Count);

		delete proteinScoringEngine;
		delete ppe;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
		delete pepB;
		delete pepA;
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseParsimony_TestingPeptideBaseSequenceCanBeBothSharedAndUnique()
	{
		std::vector<std::wstring> sequences = {L"-XYZ--ABC", L"-XYZ-EFGABC", L"-XYZ-GABC"};

		std::vector<DigestionMotif*> motifs1 =
		{
			new DigestionMotif(L"-", nullptr, 1, nullptr),
			new DigestionMotif(L"Z", nullptr, 1, nullptr)
		};
		std::vector<DigestionMotif*> motifs2 = {new DigestionMotif(L"G", nullptr, 1, nullptr)};

		auto protease = new Protease(L"test5", CleavageSpecificity::Full, nullptr, nullptr, motifs1);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		auto protease2 = new Protease(L"test6", CleavageSpecificity::Full, nullptr, nullptr, motifs2);
		ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);
		auto peptideList = std::unordered_set<PeptideWithSetModifications*>();

		auto p = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			p.push_back(&tempVar);
		}

		DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: protease2->Name, minPeptideLength: 1);

		PeptideWithSetModifications *pepA_1Dp1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_2Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepA_3Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);

		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		PeptideSpectralMatch *psmABC_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmABC_Dp2 = new PeptideSpectralMatch(pepA_2Dp2, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		psmABC_Dp2->AddOrReplace(pepA_3Dp2, 10, 0, true, std::vector<MatchedFragmentIon*>());

		// builds psm list to match to peptides
		std::vector<PeptideSpectralMatch*> psms = {psmABC_Dp1, psmABC_Dp2};

		std::for_each(psms.begin(), psms.end(), [&] (std::any h)
		{
			h::ResolveAllAmbiguities();
		});
		std::for_each(psms.begin(), psms.end(), [&] (std::any h)
		{
			h::SetFdrValues(1, 0, 0, 1, 0, 0, NAN, NAN, NAN, false);
		});

		CommonParameters tempVar4();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar4, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar5, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;

		Assert::AreEqual(2, proteinGroups.size());
		if (proteinGroups.ElementAt(0)->ProteinGroupName == L"1")
		{
			Assert::AreEqual(L"1", proteinGroups.ElementAt(0).ProteinGroupName);
			Assert::AreEqual(1, proteinGroups.ElementAt(0).UniquePeptides->Count);
			Assert::AreEqual(L"ABC", proteinGroups.ElementAt(0).UniquePeptides::ElementAt(0).FullSequence);
			Assert::AreEqual(L"2|3", proteinGroups.ElementAt(1).ProteinGroupName);
			Assert::AreEqual(0, proteinGroups.ElementAt(1).UniquePeptides->Count);
			Assert::AreEqual(2, proteinGroups.ElementAt(1).AllPeptides->Count);
		}
		else
		{
			Assert::AreEqual(L"1", proteinGroups.ElementAt(1).ProteinGroupName);
			Assert::AreEqual(1, proteinGroups.ElementAt(1).UniquePeptides->Count);
			Assert::AreEqual(L"ABC", proteinGroups.ElementAt(1).UniquePeptides::ElementAt(0).FullSequence);
			Assert::AreEqual(L"2|3", proteinGroups.ElementAt(0).ProteinGroupName);
			Assert::AreEqual(0, proteinGroups.ElementAt(0).UniquePeptides->Count);
			Assert::AreEqual(2, proteinGroups.ElementAt(0).AllPeptides->Count);
		}

		delete proteinScoringEngine;
		delete ppe;
		delete psmABC_Dp2;
		delete psmABC_Dp1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_3Dp2' statement was not added since pepA_3Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_2Dp2' statement was not added since pepA_2Dp2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepA_1Dp1' statement was not added since pepA_1Dp1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseParsimony_BaseSequenceCanBeSharedOrUniqueButOnlyUnqiuePSMSeen()
	{
		std::vector<std::wstring> sequences = {L"-XYZ--ABC", L"-XYZ-EFGABC", L"-XYZ-GABC"};

		std::vector<DigestionMotif*> motifs1 =
		{
			new DigestionMotif(L"-", nullptr, 1, nullptr),
			new DigestionMotif(L"Z", nullptr, 1, nullptr)
		};
		std::vector<DigestionMotif*> motifs2 = {new DigestionMotif(L"G", nullptr, 1, nullptr)};


		auto protease = new Protease(L"test3", CleavageSpecificity::Full, nullptr, nullptr, motifs1);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		auto protease2 = new Protease(L"test4", CleavageSpecificity::Full, nullptr, nullptr, motifs2);
		ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);

		auto peptideList = std::vector<PeptideWithSetModifications*>();
		auto p = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			p.push_back(&tempVar);
		}

		DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: protease2->Name, minPeptideLength: 1);

		for (auto protein : p)
		{
			for (auto peptide : protein->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()))
			{
//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//				switch (peptide.BaseSequence)
//ORIGINAL LINE: case "ABC":
				if (peptide->BaseSequence == L"ABC")
				{
						peptideList.push_back(peptide);
				}
			}
		}

		// builds psm list to match to peptides
		std::vector<PeptideSpectralMatch*> psms;
		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		// creates psms for specific PWSM
		for (auto pwsm : peptideList)
		{
			if (pwsm->DigestionParams == digestionParams)
			{
				PeptideSpectralMatch tempVar4(pwsm, 0, 10, 0, scan, digestionParams, new std::vector<MatchedFragmentIon*>());
				psms.push_back(&tempVar4);
			}
			psms.back().SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
			psms.back().ResolveAllAmbiguities();
		}

		CommonParameters tempVar5();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar5, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar6();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar6, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;

		Assert::AreEqual(1, proteinGroups.size());
		Assert::AreEqual(L"1", proteinGroups.ElementAt(0).ProteinGroupName);
		Assert::AreEqual(1, proteinGroups.ElementAt(0).UniquePeptides->Count);
		Assert::AreEqual(L"ABC", proteinGroups.ElementAt(0).UniquePeptides::ElementAt(0).FullSequence);

		delete proteinScoringEngine;
		delete ppe;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
		delete digestionParams2;
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::TestPSMFdrFiltering_Simulated()
	{
		std::vector<std::wstring> sequences = {L"-XYZ--ABC", L"-XYZ-EFGABC"};

		std::vector<DigestionMotif*> motifs1 =
		{
			new DigestionMotif(L"-", nullptr, 1, nullptr),
			new DigestionMotif(L"Z", nullptr, 1, nullptr)
		};
		std::vector<DigestionMotif*> motifs2 = {new DigestionMotif(L"G", nullptr, 1, nullptr)};

		auto protease = new Protease(L"testC", CleavageSpecificity::Full, nullptr, nullptr, motifs1);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		auto protease2 = new Protease(L"testD", CleavageSpecificity::Full, nullptr, nullptr, motifs2);
		ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);
		auto peptideList = std::vector<PeptideWithSetModifications*>();

		auto p = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			p.push_back(&tempVar);
		}

		DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: protease2->Name, minPeptideLength: 1);

		for (auto protein : p)
		{
			for (auto peptide : protein->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()))
			{
//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//				switch (peptide.BaseSequence)
//ORIGINAL LINE: case "ABC":
				if (peptide->BaseSequence == L"ABC")
				{
						peptideList.push_back(peptide);
				}
//ORIGINAL LINE: case "EFGABC":
				else if (peptide->BaseSequence == L"EFGABC")
				{
						peptideList.push_back(peptide);
				}
//ORIGINAL LINE: case "XYZ":
				else if (peptide->BaseSequence == L"XYZ")
				{
						peptideList.push_back(peptide);
				}
			}
			for (auto peptide : protein->Digest(digestionParams2, std::vector<Modification*>(), std::vector<Modification*>()))
			{
//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//				switch (peptide.BaseSequence)
//ORIGINAL LINE: case "ABC":
				if (peptide->BaseSequence == L"ABC")
				{
						peptideList.push_back(peptide);
				}
//ORIGINAL LINE: case "-XYZ-EFG":
				else if (peptide->BaseSequence == L"-XYZ-EFG")
				{
						peptideList.push_back(peptide);
				}
			}
		}

		// builds psm list to match to peptides
		std::vector<PeptideSpectralMatch*> psms;
		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		double goodFdr = 0.00100;
		double badFdr = 0.0200;

		// creates psms for specific PWSM
		for (int i = 0; i < peptideList.size()(); i++)
		{
			if (peptideList[i]->DigestionParams == digestionParams)
			{
				PeptideSpectralMatch tempVar4(peptideList[i], 0, 10, 0, scan, digestionParams, new std::vector<MatchedFragmentIon*>());
				psms.push_back(&tempVar4);
			}
			if (peptideList[i]->DigestionParams == digestionParams2)
			{
				PeptideSpectralMatch tempVar5(peptideList[i], 0, 10, 0, scan, digestionParams2, new std::vector<MatchedFragmentIon*>());
				psms.push_back(&tempVar5);
			}

			switch (i)
			{
				case 0:
					psms.back().SetFdrValues(0, 0, goodFdr, 0, 0, badFdr, 0, 0, 0, false);
					break;

				case 1:
					psms.back().SetFdrValues(0, 0, badFdr, 0, 0, goodFdr, 0, 0, 0, false);
					break;

				case 2:
					psms.back().SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);
					break;

				case 3:
					psms.back().SetFdrValues(0, 0, badFdr, 0, 0, badFdr, 0, 0, 0, false);
					break;

				case 4:
					psms.back().SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);
					break;

				case 5:
					psms.back().SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);
					break;
			}
			psms.back().ResolveAllAmbiguities();
		}

		CommonParameters tempVar6();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar6, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar7();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar7, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;

		psms.ElementAt(5).SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);

		//this iscopy of code that filteres psms in PostSearch Analysis Task
		auto fdrFilteredPsms = std::vector<PeptideSpectralMatch*>();
		for (auto psm : psms)
		{
			if (psm != nullptr && psm->getFdrInfo()->getQValue() <= 0.0100 && psm->getFdrInfo()->getQValueNotch() <= 0.0100)
			{
				fdrFilteredPsms.push_back(psm);
			}
		}

		Assert::AreEqual(3, fdrFilteredPsms.size());

		auto test1 = std::find(fdrFilteredPsms.begin(), fdrFilteredPsms.end(), psms.ElementAt(2)) != fdrFilteredPsms.end());
		auto test2 = std::find(fdrFilteredPsms.begin(), fdrFilteredPsms.end(), psms.ElementAt(4)) != fdrFilteredPsms.end());
		auto test3 = std::find(fdrFilteredPsms.begin(), fdrFilteredPsms.end(), psms.ElementAt(5)) != fdrFilteredPsms.end());
		auto test4 = std::find(fdrFilteredPsms.begin(), fdrFilteredPsms.end(), psms.ElementAt(0)) != fdrFilteredPsms.end());
		auto test5 = std::find(fdrFilteredPsms.begin(), fdrFilteredPsms.end(), psms.ElementAt(1)) != fdrFilteredPsms.end());
		auto test6 = std::find(fdrFilteredPsms.begin(), fdrFilteredPsms.end(), psms.ElementAt(3)) != fdrFilteredPsms.end());
		Assert::AreEqual(true, test1);
		Assert::AreEqual(true, test2);
		Assert::AreEqual(true, test3);
		Assert::AreEqual(false, test4);
		Assert::AreEqual(false, test5);
		Assert::AreEqual(false, test6);

		delete proteinScoringEngine;
		delete ppe;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::TestPSMFdrFiltering_RealFile()
	{
		SearchTask *Task1 = new SearchTask();
		CommonParameters tempVar(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, = nullptr, = nullptr, = nullptr, 1);
		Task1->setCommonParameters(&tempVar);
		SearchParameters tempVar2();
		Task1->setSearchParameters(&tempVar2);
		Task1->getSearchParameters()->setDoParsimony(true);
		Task1->getSearchParameters()->setSearchTarget(true);
		Task1->getSearchParameters()->setWritePrunedDatabase(true);
		Task1->getSearchParameters()->setSearchType(SearchType::Classic);
		std::wstring mzmlName = LR"(TestData\PrunedDbSpectra.mzml)";
		std::wstring fastaName = LR"(TestData\DbForPrunedDb.fasta)";
		std::wstring outputFolder = FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TestPSMFdrFiltering_RealFileTest)");

		auto engine = new EverythingRunnerEngine(std::vector<(std::wstring, MetaMorpheusTask)*> {(L"TestPSMFdrFiltering_RealFile", Task1)}, std::vector<std::wstring> {mzmlName}, std::vector<DbForTask*> {new DbForTask(fastaName, false)}, outputFolder);
		engine->Run();

		auto thisTaskOutputFolder = FileSystem::combine(MySetUpClass::outputFolder, LR"(TestPSMFdrFiltering_RealFile)");

		auto psms = FileSystem::combine(thisTaskOutputFolder, L"AllPSMs.psmtsv");

		Assert::AreEqual(12, File::ReadLines(psms).size()());
		auto protGroups = FileSystem::combine(thisTaskOutputFolder, L"AllProteinGroups.tsv");

		Assert::AreEqual(7, File::ReadLines(protGroups).size()());
		Directory::Delete(outputFolder, true);

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete Task1' statement was not added since Task1 was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseParsimony_TestingSameBaseSequenceSharedandUniqueMoreComplexSample()
	{
		std::vector<std::wstring> sequences = {L"-XYZ--ABC", L"-XYZ-XYGABC", L"-ABGEFG-XYZ", L"-XYZ-GEFGABC"};
		auto p = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			p.push_back(&tempVar);
		}

		std::vector<DigestionMotif*> motifs1 =
		{
			new DigestionMotif(L"-", nullptr, 1, nullptr),
			new DigestionMotif(L"Z", nullptr, 1, nullptr)
		};
		std::vector<DigestionMotif*> motifs2 = {new DigestionMotif(L"G", nullptr, 1, nullptr)};

		auto protease = new Protease(L"proteaseAlpha", CleavageSpecificity::Full, nullptr, nullptr, motifs1);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		auto protease2 = new Protease(L"proteaseBeta", CleavageSpecificity::Full, nullptr, nullptr, motifs2);
		ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);

		DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: protease2->Name, minPeptideLength: 1);

		PeptideWithSetModifications *pepABC_1Alpha = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepABC_2Beta = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepABC_4Beta = new PeptideWithSetModifications(protein: p.ElementAt(3), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 12, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABC", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepEFG_4Beta = new PeptideWithSetModifications(protein: p.ElementAt(3), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"EFG", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepEFG_3Beta = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 5, oneBasedEndResidueInProtein: 7, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"EFG", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);

		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		PeptideSpectralMatch *psmABC_Alpha = new PeptideSpectralMatch(pepABC_1Alpha, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmABC_Beta = new PeptideSpectralMatch(pepABC_2Beta, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		psmABC_Beta->AddOrReplace(pepABC_4Beta, 10, 0, true, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmEFG_Beta = new PeptideSpectralMatch(pepEFG_3Beta, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		psmEFG_Beta->AddOrReplace(pepEFG_4Beta, 10, 0, true, std::vector<MatchedFragmentIon*>());

		std::vector<PeptideSpectralMatch*> psms = {psmABC_Alpha, psmABC_Beta, psmEFG_Beta};
		std::for_each(psms.begin(), psms.end(), [&] (std::any j)
		{
			j::ResolveAllAmbiguities();
		});
		std::for_each(psms.begin(), psms.end(), [&] (std::any j)
		{
			j::SetFdrValues(1, 0, 0, 1, 0, 0, NAN, NAN, NAN, false);
		});

		std::unordered_set<DigestionParams*> digestionParamsList;
		digestionParamsList.insert(digestionParams);
		digestionParamsList.insert(digestionParams2);
		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"M", motif1);
		Modification *mod = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modVarList = {mod};

		ModificationMotif motif2;
		ModificationMotif::TryGetMotif(L"M", motif2);
		Modification *mod2 = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif2, _locationRestriction: L"Anyhwere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modFixedList = {mod};

		CommonParameters tempVar4();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar4, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar5, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;
		Assert::AreEqual(2, proteinGroups.size());
		auto proteinGroup1 = proteinGroups.Where([&] (std::any h)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmEFG_Beta;
		delete psmABC_Beta;
		delete psmABC_Alpha;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepEFG_3Beta' statement was not added since pepEFG_3Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepEFG_4Beta' statement was not added since pepEFG_4Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_4Beta' statement was not added since pepABC_4Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_2Beta' statement was not added since pepABC_2Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_1Alpha' statement was not added since pepABC_1Alpha was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
			return h->ProteinGroupName == L"1";
		}).First();
		Assert::AreEqual(1, proteinGroup1->UniquePeptides->Count);
		Assert::AreEqual(1, proteinGroup1->AllPeptides->Count);
		auto proteinGroup2 = proteinGroups.Where([&] (std::any h)
		{
		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmEFG_Beta;
		delete psmABC_Beta;
		delete psmABC_Alpha;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepEFG_3Beta' statement was not added since pepEFG_3Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepEFG_4Beta' statement was not added since pepEFG_4Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_4Beta' statement was not added since pepABC_4Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_2Beta' statement was not added since pepABC_2Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_1Alpha' statement was not added since pepABC_1Alpha was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
			return h->ProteinGroupName == L"4";
		}).First();
		Assert::AreEqual(0, proteinGroup2->UniquePeptides->Count);
		Assert::AreEqual(2, proteinGroup2->AllPeptides->Count);

		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmEFG_Beta;
		delete psmABC_Beta;
		delete psmABC_Alpha;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepEFG_3Beta' statement was not added since pepEFG_3Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepEFG_4Beta' statement was not added since pepEFG_4Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_4Beta' statement was not added since pepABC_4Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_2Beta' statement was not added since pepABC_2Beta was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_1Alpha' statement was not added since pepABC_1Alpha was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseParsimony_TestingActuallyIndistinguisableProteins()
	{
		std::vector<std::wstring> sequences = {L"KABCKXYZK", L"KXYZKABCK"};
		auto p = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			p.push_back(&tempVar);
		}

		DigestionParams *digestionParams = new DigestionParams(protease: L"trypsin", minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: L"Lys-C (don't cleave before proline)", minPeptideLength: 1);

		PeptideWithSetModifications *pepABCK_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABCK", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepABCK_2T = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABCK", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepABCK_1L = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABCK", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepABCK_2L = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"ABCK", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepXYZK_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"XYZK", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepXYZK_2T = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"XYZK", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepXYZK_1L = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"XYZK", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);
		PeptideWithSetModifications *pepXYZK_2L = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity::Full, peptideDescription: L"XYZK", missedCleavages: 0, allModsOneIsNterminus: std::unordered_map<int, Modification*>(), numFixedMods: 0);

		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		PeptideSpectralMatch *psmABCK_T = new PeptideSpectralMatch(pepABCK_1T, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		psmABCK_T->AddOrReplace(pepABCK_2T, 10, 0, true, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmABCK_L = new PeptideSpectralMatch(pepABCK_1L, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		psmABCK_L->AddOrReplace(pepABCK_2L, 10, 0, true, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmXYZK_T = new PeptideSpectralMatch(pepXYZK_1T, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		psmXYZK_T->AddOrReplace(pepXYZK_2T, 10, 0, true, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmXYZK_L = new PeptideSpectralMatch(pepXYZK_1L, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		psmXYZK_L->AddOrReplace(pepXYZK_2L, 10, 0, true, std::vector<MatchedFragmentIon*>());

		std::vector<PeptideSpectralMatch*> psms = {psmABCK_T, psmABCK_L, psmXYZK_T, psmXYZK_L};
		std::for_each(psms.begin(), psms.end(), [&] (std::any j)
		{
			j::ResolveAllAmbiguities();
		});
		std::for_each(psms.begin(), psms.end(), [&] (std::any j)
		{
			j::SetFdrValues(1, 0, 0, 1, 0, 0, NAN, NAN, NAN, false);
		});

		std::unordered_set<DigestionParams*> digestionParamsList;
		digestionParamsList.insert(digestionParams);
		digestionParamsList.insert(digestionParams2);
		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"M", motif1);
		Modification *mod = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modVarList = {mod};

		ModificationMotif motif2;
		ModificationMotif::TryGetMotif(L"M", motif2);
		Modification *mod2 = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif2, _locationRestriction: L"Anyhwere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modFixedList = {mod};

		CommonParameters tempVar4();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar4, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar5, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;
		Assert::AreEqual(1, proteinGroups.size());
		Assert::AreEqual(L"1|2", proteinGroups.ElementAt(0).ProteinGroupName);
		Assert::AreEqual(8, proteinGroups.ElementAt(0).AllPeptides->Count);

		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmXYZK_L;
		delete psmXYZK_T;
		delete psmABCK_L;
		delete psmABCK_T;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepXYZK_2L' statement was not added since pepXYZK_2L was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepXYZK_1L' statement was not added since pepXYZK_1L was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepXYZK_2T' statement was not added since pepXYZK_2T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepXYZK_1T' statement was not added since pepXYZK_1T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABCK_2L' statement was not added since pepABCK_2L was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABCK_1L' statement was not added since pepABCK_1L was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABCK_2T' statement was not added since pepABCK_2T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABCK_1T' statement was not added since pepABCK_1T was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseParsimony_TestingGreedyAlgorithm()
	{
		std::vector<std::wstring> sequences = {L"-ABC-XYZ", L"-ABC-GABC", L"-XYZ-GABC"};
		auto p = std::vector<Protein*>();
		std::vector<std::tuple<std::wstring, std::wstring>> gn;
		for (int i = 0; i < sequences.size(); i++)
		{
			Protein tempVar(sequences[i], std::to_wstring(i + 1), nullptr, gn, new std::unordered_map<int, std::vector<Modification*>>());
			p.push_back(&tempVar);
		}

		std::vector<DigestionMotif*> motifs1 =
		{
			new DigestionMotif(L"-", nullptr, 0, nullptr),
			new DigestionMotif(L"-", nullptr, 1, nullptr)
		};
		std::vector<DigestionMotif*> motifs2 = {new DigestionMotif(L"G", nullptr, 1, nullptr)};

		auto protease = new Protease(L"proteaseDash", CleavageSpecificity::Full, nullptr, nullptr, motifs1);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		auto protease2 = new Protease(L"proteaseG", CleavageSpecificity::Full, nullptr, nullptr, motifs2);
		ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);

		DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 1);
		DigestionParams *digestionParams2 = new DigestionParams(protease: protease2->Name, minPeptideLength: 1);

		PeptideWithSetModifications *pepABC_1Dash = new PeptideWithSetModifications(p.ElementAt(0), digestionParams, 2, 4, CleavageSpecificity::Unknown, L"ABCK", 0, std::unordered_map<int, Modification*>(), 0);
		PeptideWithSetModifications *pepABC_2Dash = new PeptideWithSetModifications(p.ElementAt(1), digestionParams, 7, 9, CleavageSpecificity::Unknown, L"ABCK", 0, std::unordered_map<int, Modification*>(), 0);
		PeptideWithSetModifications *pepABC_2G = new PeptideWithSetModifications(p.ElementAt(1), digestionParams2, 2, 4, CleavageSpecificity::Unknown, L"ABCK", 0, std::unordered_map<int, Modification*>(), 0);
		PeptideWithSetModifications *pepABC_3G = new PeptideWithSetModifications(p.ElementAt(2), digestionParams2, 7, 9, CleavageSpecificity::Unknown, L"ABCK", 0, std::unordered_map<int, Modification*>(), 0);

		MzSpectrum tempVar2(new double[] {1}, new double[] {1}, false);
		MsDataScan *dfb = new MsDataScan(&tempVar2, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar3);

		PeptideSpectralMatch *psmABC_Dash = new PeptideSpectralMatch(pepABC_1Dash, 0, 10, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		psmABC_Dash->AddOrReplace(pepABC_2Dash, 10, 0, true, std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *psmABC_G = new PeptideSpectralMatch(pepABC_2G, 0, 10, 0, scan, digestionParams2, std::vector<MatchedFragmentIon*>());
		psmABC_G->AddOrReplace(pepABC_3G, 10, 0, true, std::vector<MatchedFragmentIon*>());

		std::vector<PeptideSpectralMatch*> psms = {psmABC_Dash, psmABC_G};
		std::for_each(psms.begin(), psms.end(), [&] (std::any j)
		{
			j::ResolveAllAmbiguities();
		});
		std::for_each(psms.begin(), psms.end(), [&] (std::any j)
		{
			j::SetFdrValues(1, 0, 0, 1, 0, 0, NAN, NAN, NAN, false);
		});

		std::unordered_set<DigestionParams*> digestionParamsList;
		digestionParamsList.insert(digestionParams);
		digestionParamsList.insert(digestionParams2);
		ModificationMotif motif1;
		ModificationMotif::TryGetMotif(L"M", motif1);
		Modification *mod = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif1, _locationRestriction: L"Anywhere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modVarList = {mod};

		ModificationMotif motif2;
		ModificationMotif::TryGetMotif(L"M", motif2);
		Modification *mod2 = new Modification(_originalId: L"Oxidation of M", _modificationType: L"Common Variable", _target: motif2, _locationRestriction: L"Anyhwere.", _monoisotopicMass: 15.99491461957);
		std::vector<Modification*> modFixedList = {mod};

		CommonParameters tempVar4();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(psms, false, &tempVar4, nullptr);
		auto proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(ppe->Run());

		// score protein groups and merge indistinguishable ones
		CommonParameters tempVar5();
		ProteinScoringAndFdrEngine *proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults->getProteinGroups(), psms, false, true, true, &tempVar5, std::vector<std::wstring>());
		auto results = static_cast<ProteinScoringAndFdrResults*>(proteinScoringEngine->Run());

		std::vector<ProteinGroup*> &proteinGroups = results->SortedAndScoredProteinGroups;
		Assert::AreEqual(1, proteinGroups.size());
		Assert::AreEqual(L"2", proteinGroups.ElementAt(0).ProteinGroupName);
		Assert::AreEqual(2, proteinGroups.ElementAt(0).AllPeptides->Count);

		delete proteinScoringEngine;
		delete ppe;
		delete mod2;
		delete mod;
		delete psmABC_G;
		delete psmABC_Dash;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_3G' statement was not added since pepABC_3G was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_2G' statement was not added since pepABC_2G was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_2Dash' statement was not added since pepABC_2Dash was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pepABC_1Dash' statement was not added since pepABC_1Dash was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams2' statement was not added since digestionParams2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}

	void MultiProteaseParsimonyTest::MultiProteaseParsimony_TestingProteaseSpecificFDRCalculations()
	{
		// two protease options
		DigestionParams *tryp = new DigestionParams(protease: L"trypsin");
		DigestionParams *gluC = new DigestionParams(protease: L"Glu-C");

		// target or decoy protein
		Protein *t = new Protein(L"P", L"1");
		Protein *d = new Protein(L"P", L"2", isDecoy: true);

		MzSpectrum tempVar(new double[1], new double[1], false);
		MsDataScan *dfb = new MsDataScan(&tempVar, 0, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 0, nullptr);

		CommonParameters tempVar2();
		Ms2ScanWithSpecificMass *scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, L"File", &tempVar2);
		std::vector<MatchedFragmentIon*> f;

		PeptideWithSetModifications tempVar3(L"P", nullptr, digestionParams: tryp, p: t);
		PeptideWithSetModifications tempVar12(L"P", nullptr, digestionParams: tryp, p: t);
		PeptideWithSetModifications tempVar11(L"P", nullptr, digestionParams: tryp, p: d);
		PeptideWithSetModifications tempVar10(L"P", nullptr, digestionParams: tryp, p: d);
		PeptideWithSetModifications tempVar9(L"P", nullptr, digestionParams: tryp, p: t);
		PeptideWithSetModifications tempVar8(L"P", nullptr, digestionParams: gluC, p: t);
		PeptideWithSetModifications tempVar7(L"P", nullptr, digestionParams: gluC, p: d);
		PeptideWithSetModifications tempVar6(L"P", nullptr, digestionParams: gluC, p: t);
		PeptideWithSetModifications tempVar5(L"P", nullptr, digestionParams: tryp, p: t);
		PeptideWithSetModifications tempVar4(L"P", nullptr, digestionParams: gluC, p: t);
		std::vector<PeptideSpectralMatch*> psms =
		{
			new PeptideSpectralMatch(&tempVar3, 0, 20, 1, scan, tryp, f),
			new PeptideSpectralMatch(&tempVar4, 0, 19, 1, scan, gluC, f),
			new PeptideSpectralMatch(&tempVar5, 0, 18, 1, scan, tryp, f),
			new PeptideSpectralMatch(&tempVar6, 0, 17, 1, scan, gluC, f),
			new PeptideSpectralMatch(&tempVar7, 0, 16, 1, scan, gluC, f),
			new PeptideSpectralMatch(&tempVar8, 0, 15, 1, scan, gluC, f),
			new PeptideSpectralMatch(&tempVar9, 0, 14, 1, scan, tryp, f),
			new PeptideSpectralMatch(&tempVar10, 0, 13, 1, scan, tryp, f),
			new PeptideSpectralMatch(&tempVar11, 0, 12, 1, scan, tryp, f),
			new PeptideSpectralMatch(&tempVar12, 0, 11, 1, scan, tryp, f)
		};

		std::for_each(psms.begin(), psms.end(), [&] (std::any p)
		{
			p::ResolveAllAmbiguities();
		});

		FdrAnalysisEngine tempVar13(psms, 0, new CommonParameters(), new std::vector<std::wstring>());
		(&tempVar13)->Run();
		psms = psms.OrderByDescending([&] (std::any p)
		{
			p::Score;
		}).ToList();

		Assert::AreEqual(0.00, std::round(psms[0]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));
		Assert::AreEqual(0.00, std::round(psms[1]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));
		Assert::AreEqual(0.00, std::round(psms[2]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));
		Assert::AreEqual(0.00, std::round(psms[3]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));
		Assert::AreEqual(0.33, std::round(psms[4]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));
		Assert::AreEqual(0.33, std::round(psms[5]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));
		Assert::AreEqual(0.00, std::round(psms[6]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));
		Assert::AreEqual(0.33, std::round(psms[7]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));
		Assert::AreEqual(0.50, std::round(psms[8]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));
		Assert::AreEqual(0.50, std::round(psms[9]->getFdrInfo().getQValue() * std::pow(10, 2)) / std::pow(10, 2));

//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete dfb' statement was not added since dfb was passed to a method or constructor. Handle memory management manually.
		delete d;
		delete t;
//C# TO C++ CONVERTER TODO TASK: A 'delete gluC' statement was not added since gluC was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete tryp' statement was not added since tryp was passed to a method or constructor. Handle memory management manually.
	}
}
