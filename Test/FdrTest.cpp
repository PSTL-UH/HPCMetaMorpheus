#include "FdrTest.h"
#include "../EngineLayer/PrecursorSearchModes/DotMassDiffAcceptor.h"
#include "../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "TestDataFile.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../TaskLayer/SearchTask/SearchParameters.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/SearchTask/SearchTask.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::Indexing;
using namespace EngineLayer::ModernSearch;
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

	void FdrTest::FdrTestMethod()
	{
		PpmTolerance tempVar(5);
		MassDiffAcceptor *searchModes = new DotMassDiffAcceptor(L"", {0, 1.0029}, &tempVar);
		std::vector<std::wstring> nestedIds;

		Protein *p = new Protein(L"MNKNNKNNNKNNNNK", nullptr);
		DigestionParams *digestionParams = new DigestionParams();
		auto digested = p->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()).ToList();

		PeptideWithSetModifications *pep1 = digested[0];
		PeptideWithSetModifications *pep2 = digested[1];
		PeptideWithSetModifications *pep3 = digested[2];
		PeptideWithSetModifications *pep4 = digested[3];

		TestDataFile *t = new TestDataFile(std::vector<PeptideWithSetModifications*> {pep1, pep2, pep3});

		MsDataScan *mzLibScan1 = t->GetOneBasedScan(2);
		CommonParameters tempVar2();
		Ms2ScanWithSpecificMass *scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, pep1->MonoisotopicMass.ToMz(1), 1, L"", &tempVar2);
		PeptideSpectralMatch *psm1 = new PeptideSpectralMatch(pep1, 0, 3, 0, scan1, digestionParams, std::vector<MatchedFragmentIon*>());

		MsDataScan *mzLibScan2 = t->GetOneBasedScan(4);
		CommonParameters tempVar3();
		Ms2ScanWithSpecificMass *scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, pep2->MonoisotopicMass.ToMz(1), 1, L"", &tempVar3);
		PeptideSpectralMatch *psm2 = new PeptideSpectralMatch(pep2, 1, 2, 1, scan2, digestionParams, std::vector<MatchedFragmentIon*>());

		MsDataScan *mzLibScan3 = t->GetOneBasedScan(6);
		CommonParameters tempVar4();
		Ms2ScanWithSpecificMass *scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, pep3->MonoisotopicMass.ToMz(1), 1, L"", &tempVar4);
		PeptideSpectralMatch *psm3 = new PeptideSpectralMatch(pep3, 0, 1, 2, scan3, digestionParams, std::vector<MatchedFragmentIon*>());

		psm3->AddOrReplace(pep4, 1, 1, true, std::vector<MatchedFragmentIon*>());

		auto newPsms = std::vector<PeptideSpectralMatch*> {psm1, psm2, psm3};
		for (auto psm : newPsms)
		{
			psm->ResolveAllAmbiguities();
		}


		CommonParameters *cp = new CommonParameters(, = DissociationType::HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, true);

		FdrAnalysisEngine *fdr = new FdrAnalysisEngine(newPsms, searchModes->NumNotches, cp, nestedIds);

		fdr->Run();

		Assert::AreEqual(2, searchModes->NumNotches);
		Assert::AreEqual(0, newPsms[0]->getFdrInfo().getCumulativeDecoyNotch());
		Assert::AreEqual(1, newPsms[0]->getFdrInfo().getCumulativeTargetNotch());
		Assert::AreEqual(0, newPsms[1]->getFdrInfo().getCumulativeDecoyNotch());
		Assert::AreEqual(1, newPsms[1]->getFdrInfo().getCumulativeTargetNotch());
		Assert::AreEqual(0, newPsms[2]->getFdrInfo().getCumulativeDecoyNotch());
		Assert::AreEqual(1, newPsms[2]->getFdrInfo().getCumulativeTargetNotch());

		Assert::AreEqual(0, newPsms[0]->getFdrInfo().getCumulativeDecoy());
		Assert::AreEqual(1, newPsms[0]->getFdrInfo().getCumulativeTarget());
		Assert::AreEqual(0, newPsms[1]->getFdrInfo().getCumulativeDecoy());
		Assert::AreEqual(2, newPsms[1]->getFdrInfo().getCumulativeTarget());
		Assert::AreEqual(0, newPsms[2]->getFdrInfo().getCumulativeDecoy());
		Assert::AreEqual(3, newPsms[2]->getFdrInfo().getCumulativeTarget());

		delete fdr;
//C# TO C++ CONVERTER TODO TASK: A 'delete cp' statement was not added since cp was passed to a method or constructor. Handle memory management manually.
		delete psm3;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan3' statement was not added since scan3 was passed to a method or constructor. Handle memory management manually.
		delete psm2;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan2' statement was not added since scan2 was passed to a method or constructor. Handle memory management manually.
		delete psm1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan1' statement was not added since scan1 was passed to a method or constructor. Handle memory management manually.
		delete t;
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
		delete p;
		delete searchModes;
	}

	void FdrTest::TestDeltaValues()
	{
		DigestionParams tempVar(minPeptideLength: 5);
		CommonParameters *CommonParameters = new CommonParameters(scoreCutoff: 1, useDeltaScore: true, digestionParams: &tempVar);

		SearchParameters *SearchParameters = new SearchParameters();
		SearchParameters->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
		std::vector<Modification*> variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			CommonParameters->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();
		std::vector<Modification*> fixedModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			CommonParameters->ListOfModsFixed->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();

		// Generate data for files
		Protein *TargetProtein1 = new Protein(L"TIDEANTHE", L"accession1");
		Protein *TargetProtein2 = new Protein(L"TIDELVE", L"accession2");
		Protein *TargetProtein3 = new Protein(L"TIDENIE", L"accession3");
		Protein *TargetProteinLost = new Protein(L"PEPTIDEANTHE", L"accession4");
		Protein *DecoyProteinFound = new Protein(L"PETPLEDQGTHE", L"accessiond", isDecoy: true);

		MsDataFile *myMsDataFile = new TestDataFile(std::vector<PeptideWithSetModifications*> {TargetProtein1->Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).ToList()[0], TargetProtein2->Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).ToList()[0], TargetProtein3->Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).ToList()[0], DecoyProteinFound->Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).ToList()[0]});

		auto proteinList = std::vector<Protein*> {TargetProtein1, TargetProtein2, TargetProtein3, TargetProteinLost, DecoyProteinFound};

		auto searchModes = new SinglePpmAroundZeroSearchMode(5);

		bool DoPrecursorDeconvolution = true;
		bool UseProvidedPrecursorInfo = true;
		double DeconvolutionIntensityRatio = 4;
		int DeconvolutionMaxAssumedChargeState = 10;
		Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);

		CommonParameters tempVar2();
		auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, L"", &tempVar2).OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		//check better when using delta
		std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
		ClassicSearchEngine tempVar3(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new std::vector<std::wstring>());
		(&tempVar3)->Run();

		auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::None, CommonParameters, 30000, false, std::vector<FileInfo*>(), std::vector<std::wstring>());
		auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
		MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(), SearchParameters->getMassDiffAcceptorType(), SearchParameters->getCustomMdac());

		std::vector<PeptideSpectralMatch*> allPsmsArrayModern(listOfSortedms2Scans.size());
		ModernSearchEngine tempVar4(allPsmsArrayModern, listOfSortedms2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(), 0, CommonParameters, massDiffAcceptor, 0, new std::vector<std::wstring>());
		(&tempVar4)->Run();

		FdrAnalysisEngine tempVar5(allPsmsArray.ToList(), 1, CommonParameters, new std::vector<std::wstring>());
		FdrAnalysisResults *fdrResultsClassicDelta = static_cast<FdrAnalysisResults*>((&tempVar5)->Run());
		FdrAnalysisEngine tempVar6(allPsmsArrayModern.ToList(), 1, CommonParameters, new std::vector<std::wstring>());
		FdrAnalysisResults *fdrResultsModernDelta = static_cast<FdrAnalysisResults*>((&tempVar6)->Run());
		Assert::IsTrue(fdrResultsClassicDelta->getPsmsWithin1PercentFdr() == 3);
		Assert::IsTrue(fdrResultsModernDelta->getPsmsWithin1PercentFdr() == 3);

		DigestionParams tempVar7(minPeptideLength: 5);
		CommonParameters = new CommonParameters(digestionParams: &tempVar7);

		//check worse when using score
		FdrAnalysisEngine tempVar8(allPsmsArray.ToList(), 1, CommonParameters, new std::vector<std::wstring>());
		FdrAnalysisResults *fdrResultsClassic = static_cast<FdrAnalysisResults*>((&tempVar8)->Run());
		FdrAnalysisEngine tempVar9(allPsmsArray.ToList(), 1, CommonParameters, new std::vector<std::wstring>());
		FdrAnalysisResults *fdrResultsModern = static_cast<FdrAnalysisResults*>((&tempVar9)->Run());
		Assert::IsTrue(fdrResultsClassic->getPsmsWithin1PercentFdr() == 0);
		Assert::IsTrue(fdrResultsModern->getPsmsWithin1PercentFdr() == 0);

		//check that when delta is bad, we used the score
		// Generate data for files
		Protein *DecoyProtein1 = new Protein(L"TLEDAGGTHE", L"accession1d", isDecoy: true);
		Protein *DecoyProtein2 = new Protein(L"TLEDLVE", L"accession2d", isDecoy: true);
		Protein *DecoyProtein3 = new Protein(L"TLEDNIE", L"accession3d", isDecoy: true);
		Protein *DecoyProteinShiny = new Protein(L"GGGGGG", L"accessionShinyd", isDecoy: true);

		myMsDataFile = new TestDataFile(std::vector<PeptideWithSetModifications*> {TargetProtein1->Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).ToList()[0], TargetProtein2->Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).ToList()[0], TargetProtein3->Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).ToList()[0], DecoyProteinShiny->Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).ToList()[0]});

		proteinList = {TargetProtein1, DecoyProtein1, TargetProtein2, DecoyProtein2, TargetProtein3, DecoyProtein3, DecoyProteinShiny};

		CommonParameters tempVar10();
		listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, L"", &tempVar10).OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		//check no change when using delta
		allPsmsArray = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
		ClassicSearchEngine tempVar11(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new std::vector<std::wstring>());
		(&tempVar11)->Run();

		DigestionParams tempVar12(minPeptideLength: 5);
		CommonParameters = new CommonParameters(useDeltaScore: true, digestionParams: &tempVar12);

		indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::None, CommonParameters, 30000, false, std::vector<FileInfo*>(), std::vector<std::wstring>());
		indexResults = static_cast<IndexingResults*>(indexEngine->Run());
		massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(), SearchParameters->getMassDiffAcceptorType(), SearchParameters->getCustomMdac());
		allPsmsArrayModern = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
		ModernSearchEngine tempVar13(allPsmsArrayModern, listOfSortedms2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(), 0, CommonParameters, massDiffAcceptor, 0, new std::vector<std::wstring>());
		(&tempVar13)->Run();

		FdrAnalysisEngine tempVar14(allPsmsArray.ToList(), 1, CommonParameters, new std::vector<std::wstring>());
		fdrResultsClassicDelta = static_cast<FdrAnalysisResults*>((&tempVar14)->Run());
		FdrAnalysisEngine tempVar15(allPsmsArrayModern.ToList(), 1, CommonParameters, new std::vector<std::wstring>());
		fdrResultsModernDelta = static_cast<FdrAnalysisResults*>((&tempVar15)->Run());
		Assert::IsTrue(fdrResultsClassicDelta->getPsmsWithin1PercentFdr() == 3);
		Assert::IsTrue(fdrResultsModernDelta->getPsmsWithin1PercentFdr() == 3);

		DigestionParams tempVar16(minPeptideLength: 5);
		CommonParameters = new CommonParameters(digestionParams: &tempVar16);

		//check no change when using score
		FdrAnalysisEngine tempVar17(allPsmsArray.ToList(), 1, CommonParameters, new std::vector<std::wstring>());
		fdrResultsClassic = static_cast<FdrAnalysisResults*>((&tempVar17)->Run());
		FdrAnalysisEngine tempVar18(allPsmsArrayModern.ToList(), 1, CommonParameters, new std::vector<std::wstring>());
		fdrResultsModern = static_cast<FdrAnalysisResults*>((&tempVar18)->Run());
		Assert::IsTrue(fdrResultsClassic->getPsmsWithin1PercentFdr() == 3);
		Assert::IsTrue(fdrResultsModern->getPsmsWithin1PercentFdr() == 3);

		delete DecoyProteinShiny;
		delete DecoyProtein3;
		delete DecoyProtein2;
		delete DecoyProtein1;
		delete indexEngine;
		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
		delete DecoyProteinFound;
		delete TargetProteinLost;
		delete TargetProtein3;
		delete TargetProtein2;
		delete TargetProtein1;
		delete SearchParameters;
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
	}
}
