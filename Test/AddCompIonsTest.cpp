#include "AddCompIonsTest.h"
#include "TestDataFile.h"
#include "../EngineLayer/PrecursorSearchModes/OpenMassDiffAcceptor.h"
#include "../EngineLayer/CommonParameters.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../TaskLayer/SearchTask/SearchParameters.h"
#include "../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/MetaMorpheusEngine.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
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

	void AddCompIonsTest::TestAddCompIonsClassic()
	{
		auto myMsDataFile = new TestDataFile();
		auto variableModifications = std::vector<Modification*>();
		auto fixedModifications = std::vector<Modification*>();
		auto proteinList = std::vector<Protein*> {new Protein(L"QXQ", nullptr)};

		auto productMassTolerance = new AbsoluteTolerance(0.01);
		auto searchModes = new OpenSearchMode();

		bool DoPrecursorDeconvolution = true;
		bool UseProvidedPrecursorInfo = true;
		double DeconvolutionIntensityRatio = 4;
		int DeconvolutionMaxAssumedChargeState = 10;
		Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);

		CommonParameters tempVar();
		auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, L"", &tempVar).OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		std::vector<DigestionMotif*> motifs = {new DigestionMotif(L"K", nullptr, 1, nullptr)};
		Protease *protease = new Protease(L"Custom Protease3", CleavageSpecificity::Full, nullptr, nullptr, motifs);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);

		DigestionParams tempVar2(protease: protease->Name, maxMissedCleavages: 0, minPeptideLength: 1);
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar2, scoreCutoff: 1, addCompIons: false);
		std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
		ClassicSearchEngine tempVar3(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new std::vector<std::wstring>());
		(&tempVar3)->Run();

		DigestionParams tempVar4(protease: protease->Name, maxMissedCleavages: 0, minPeptideLength: 1);
		CommonParameters *CommonParameters2 = new CommonParameters(digestionParams: &tempVar4, scoreCutoff: 1, addCompIons: true);

		std::vector<PeptideSpectralMatch*> allPsmsArray2(listOfSortedms2Scans.size());
		ClassicSearchEngine tempVar5(allPsmsArray2, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters2, new std::vector<std::wstring>());
		(&tempVar5)->Run();

		double scoreT = allPsmsArray2[0]->getScore();
		double scoreF = allPsmsArray[0]->getScore();

		// Single search mode
		Assert::AreEqual(allPsmsArray.size(), allPsmsArray2.size());

		// Single ms2 scan
		Assert::AreEqual(allPsmsArray.size(), allPsmsArray2.size());

		Assert::IsTrue(scoreT > 1);

		Assert::AreEqual(allPsmsArray[0]->getScanNumber(), allPsmsArray2[0]->getScanNumber());

		Assert::IsTrue(scoreT == scoreF * 3 && scoreT > scoreF + 2);

//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters2' statement was not added since CommonParameters2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes was passed to a method or constructor. Handle memory management manually.
		delete productMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
	}

	void AddCompIonsTest::TestCompIons_ModernSearch()
	{
		auto myMsDataFile = new TestDataFile();
		auto variableModifications = std::vector<Modification*>();
		auto fixedModifications = std::vector<Modification*>();
		auto localizeableModifications = std::vector<Modification*>();
		std::unordered_map<Modification*, unsigned short> modsDictionary;
		for (auto mod : fixedModifications)
		{
			modsDictionary.emplace(mod, 0);
		}
		int ii = 1;
		for (auto mod : variableModifications)
		{
			modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
			ii++;
		}
		for (auto mod : localizeableModifications)
		{
			modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
			ii++;
		}

		auto proteinList = std::vector<Protein*> {new Protein(L"MNNNKQQQ", nullptr)};

		SearchParameters *SearchParameters = new SearchParameters();
		SearchParameters->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
		SearchParameters->setSearchTarget(true);
		std::vector<DigestionMotif*> motifs = {new DigestionMotif(L"K", nullptr, 1, nullptr)};
		Protease *protease = new Protease(L"singleN4", CleavageSpecificity::Full, nullptr, nullptr, motifs);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		DigestionParams tempVar(protease: protease->Name, minPeptideLength: 1);
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1);
		DigestionParams tempVar2(protease: protease->Name, minPeptideLength: 1);
		CommonParameters *withCompIons = new CommonParameters(digestionParams: &tempVar2, scoreCutoff: 1, addCompIons: true);

		auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters, SearchParameters->getMaxFragmentSize(), false, std::vector<FileInfo*>(), std::vector<std::wstring>());

		auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());

		bool DoPrecursorDeconvolution = true;
		bool UseProvidedPrecursorInfo = true;
		double DeconvolutionIntensityRatio = 4;
		int DeconvolutionMaxAssumedChargeState = 10;
		Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);

		CommonParameters tempVar3();
		auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, L"", &tempVar3).OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(), SearchParameters->getMassDiffAcceptorType(), SearchParameters->getCustomMdac());

		// without complementary ions
		std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
		ModernSearchEngine tempVar4(allPsmsArray, listOfSortedms2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(), 0, CommonParameters, massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(), new std::vector<std::wstring>());
		(&tempVar4)->Run();

		// with complementary ions
		std::vector<PeptideSpectralMatch*> allPsmsArray2(listOfSortedms2Scans.size());
		ModernSearchEngine tempVar5(allPsmsArray2, listOfSortedms2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(), 0, withCompIons, massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(), new std::vector<std::wstring>());
		(&tempVar5)->Run();

		// Single search mode
		Assert::AreEqual(allPsmsArray.size(), allPsmsArray2.size());

		// Single ms2 scan
		Assert::AreEqual(allPsmsArray.size(), allPsmsArray2.size());
		Assert::That(allPsmsArray[0] != nullptr);
		Assert::That(allPsmsArray2[0] != nullptr);

		Assert::IsTrue(allPsmsArray2[0]->getScore() > 1);

		Assert::AreEqual(allPsmsArray[0]->getScanNumber(), allPsmsArray2[0]->getScanNumber());

		Assert::IsTrue(allPsmsArray2[0]->Score <= allPsmsArray[0]->getScore() * 2 && allPsmsArray2[0]->getScore() > allPsmsArray[0]->getScore() + 3);

		delete DeconvolutionMassTolerance;
		delete indexEngine;
//C# TO C++ CONVERTER TODO TASK: A 'delete withCompIons' statement was not added since withCompIons was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
		delete SearchParameters;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
	}

	void AddCompIonsTest::TestCompIons_MatchIonsScore()
	{
		TestDataFile *t = new TestDataFile();
		Tolerance *productMassTolerance = new AbsoluteTolerance(0.01);
		double precursorMass = 300;
		//The below theoretical does not accurately represent B-Y ions
		std::vector<double> sorted_theoretical_product_masses_for_this_peptide = {precursorMass + (2 * Constants::ProtonMass) - 275.1350, precursorMass + (2 * Constants::ProtonMass) - 258.127, precursorMass + (2 * Constants::ProtonMass) - 257.1244, 50, 60, 70, 147.0764, precursorMass + (2 * Constants::ProtonMass) - 147.0764, precursorMass + (2 * Constants::ProtonMass) - 70, precursorMass + (2 * Constants::ProtonMass) - 60, precursorMass + (2 * Constants::ProtonMass) - 50, 257.1244, 258.127, 275.1350}; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
		std::vector<Product*> productsWithLocalizedMassDiff;
		for (auto d : sorted_theoretical_product_masses_for_this_peptide)
		{
			NeutralTerminusFragment *frag = new NeutralTerminusFragment(FragmentationTerminus::Both, d, 1, 1);
			Product tempVar(ProductType::b, frag, 0);
			productsWithLocalizedMassDiff.push_back(&tempVar);

//C# TO C++ CONVERTER TODO TASK: A 'delete frag' statement was not added since frag was passed to a method or constructor. Handle memory management manually.
		}
		CommonParameters *commonParametersNoComp = new CommonParameters();
		AbsoluteTolerance tempVar2(0.01);
		commonParametersNoComp->setProductMassTolerance(&tempVar2);
		AbsoluteTolerance tempVar3(0.01);
		CommonParameters *commonParametersWithComp = new CommonParameters(, , = true, = true, = 3, = 12, = true, true, , , , , , , , , &tempVar3);

		MsDataScan *scan = t->GetOneBasedScan(2);
		CommonParameters tempVar4();
		auto scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(1), 1, L"", &tempVar4);
		std::vector<MatchedFragmentIon*> matchedIons = MetaMorpheusEngine::MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParametersNoComp);

		std::vector<MatchedFragmentIon*> matchedCompIons = MetaMorpheusEngine::MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParametersWithComp);
		matchedCompIons.insert(matchedCompIons.end(), matchedIons.begin(), matchedIons.end());

		// score when the mass-diff is on this residue
		double localizedScore = MetaMorpheusEngine::CalculatePeptideScore(scan, matchedIons, 0);
		double scoreNormal = MetaMorpheusEngine::CalculatePeptideScore(scan, matchedIons, 0);
		double scoreComp = MetaMorpheusEngine::CalculatePeptideScore(scan, matchedCompIons, 0);
		Assert::IsTrue(scoreNormal * 2 == scoreComp && scoreComp > scoreNormal + 1);

//C# TO C++ CONVERTER TODO TASK: A 'delete scanWithMass' statement was not added since scanWithMass was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete commonParametersWithComp' statement was not added since commonParametersWithComp was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete commonParametersNoComp' statement was not added since commonParametersNoComp was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete productMassTolerance' statement was not added since productMassTolerance was passed to a method or constructor. Handle memory management manually.
		delete t;
	}

	void AddCompIonsTest::TestCompIons_MatchIons()
	{
		TestDataFile *t = new TestDataFile(0.0001);
		Tolerance *productMassTolerance = new AbsoluteTolerance(0.01);
		double precursorMass = 402.18629720155;
		//The below theoretical does not accurately represent B-Y ions
		std::vector<double> sorted_theoretical_product_masses_for_this_peptide = {50, 60, 70, 147.0764 - Constants::ProtonMass, 200, 215, 230, 245, precursorMass + Constants::ProtonMass - 147.0764, 258.127, 275.1350, precursorMass + (2 * Constants::ProtonMass) - 70, precursorMass + (2 * Constants::ProtonMass) - 60, precursorMass + (2 * Constants::ProtonMass) - 50}; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
		std::vector<double> matchedIonMassesT;
		std::vector<double> matchedDaErrorT;
		std::vector<double> matchedPpmErrorT;
		std::vector<double> matchedIonIntensityT;
		std::vector<double> matchedIonMassesF;
		std::vector<double> matchedDaErrorF;
		std::vector<double> matchedPpmErrorF;
		std::vector<double> matchedIonIntensityF;

		std::vector<int> matchedIonSeriesT;
		std::vector<int> matchedIonSeriesF;

		//MetaMorpheusEngine.MatchIons(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonSeriesT, matchedIonMassesT, matchedDaErrorT, matchedPpmErrorT, matchedIonIntensityT, precursorMass, ProductType.B, true);
		//MetaMorpheusEngine.MatchIons(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonSeriesF, matchedIonMassesF, matchedDaErrorF, matchedPpmErrorF, matchedIonIntensityF, precursorMass, ProductType.B, false);

		//Test the number of series is doubled
		Assert::IsTrue(matchedIonSeriesT.size() == matchedIonSeriesF.size() * 2);
		//Test the number of ions is doubled
		Assert::IsTrue(matchedIonMassesT.size() == matchedIonMassesF.size() * 2);
		//Test the number of da errors is doubled
		Assert::IsTrue(matchedDaErrorT.size() == matchedDaErrorF.size() * 2);
		//test the number of ppm errors is doubled
		Assert::IsTrue(matchedPpmErrorT.size() == matchedPpmErrorF.size() * 2);
		//test the number of the intensity values is doubled
		Assert::IsTrue(matchedIonIntensityT.size() == matchedIonIntensityF.size() * 2);
		for (auto d : matchedDaErrorF)
		{
			Assert::IsTrue(d <= 0.01);
		}
		for (auto d : matchedDaErrorT)
		{
			Assert::IsTrue(d <= 0.01);
		}

		delete productMassTolerance;
		delete t;
	}
}
