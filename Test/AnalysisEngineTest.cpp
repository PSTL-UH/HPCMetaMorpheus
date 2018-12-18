#include "AnalysisEngineTest.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "TestDataFile.h"
#include "../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../TaskLayer/MetaMorpheusTask.h"

using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::HistogramAnalysis;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;

namespace Test
{

	void AnalysisEngineTests::TestAnalysisEngineTests()
	{
		std::vector<DigestionMotif*> motifs = {new DigestionMotif(L"K", nullptr, 1, nullptr)};
		Protease *protease = new Protease(L"Custom Protease5", CleavageSpecificity::Full, nullptr, nullptr, motifs);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		DigestionParams tempVar(protease: protease->Name, maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: 1042);
		PpmTolerance tempVar2(10);
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1, productMassTolerance: &tempVar2);

		std::vector<Modification*> localizeableModifications;
		std::vector<Modification*> variableModifications;
		std::vector<Modification*> fixedModifications;

		std::unordered_map<Modification*, unsigned short> modsDictionary;
		for (auto mod : fixedModifications)
		{
			modsDictionary.emplace(mod, 0);
		}
		int i = 1;
		for (auto mod : variableModifications)
		{
			modsDictionary.emplace(mod, static_cast<unsigned short>(i));
			i++;
		}
		for (auto mod : localizeableModifications)
		{
			modsDictionary.emplace(mod, static_cast<unsigned short>(i));
			i++;
		}

		auto proteinList = std::vector<Protein*> {new Protein(L"MNNNKQQQ", L"accession")};
		auto modPep = proteinList.front().Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).Last();
		std::unordered_set<PeptideWithSetModifications*> value1 = {modPep};
		PeptideWithSetModifications *compactPeptide1 = value1.First();

		Assert::AreEqual(L"QQQ", value1.First().BaseSequence);
		auto modPep2 = proteinList.front().Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).First();
		std::unordered_set<PeptideWithSetModifications*> value2 = {modPep2};
		PeptideWithSetModifications *compactPeptide2 = value2.First();

		Assert::AreEqual(L"MNNNK", value2.First().BaseSequence);

		auto modPep3 = proteinList.front().Digest(CommonParameters->getDigestionParams(), fixedModifications, variableModifications).ToList()[1];
		std::unordered_set<PeptideWithSetModifications*> value3 = {modPep3};
		PeptideWithSetModifications *compactPeptide3 = value3.First();
		Assert::AreEqual(L"NNNK", value3.First().BaseSequence);

		//newPsms[0] = new List<PsmParent>[] { new List<PsmParent>{ new PsmModern(compactPeptide1, null, 1,  1, 2, 2, 1,1, 1, 1, 3,0) },
		//                                     new List<PsmParent>{  new PsmModern(compactPeptide2, null, 2,2+132.040,3,3,2,2,2,2,2,0) },
		//                                     new List<PsmParent>{ new PsmModern(compactPeptide3, null, 3, 3, 4, 3, 3, 3, 3, 3, 3, 0)} };

		MsDataScan tempVar3(new MzSpectrum(new double[] {1}, new double[] {1}, false), 2, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=1", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 1, nullptr);
		CommonParameters tempVar4();
		Ms2ScanWithSpecificMass *scanA = new Ms2ScanWithSpecificMass(&tempVar3, 1, 1, L"", &tempVar4);
		MsDataScan tempVar5(new MzSpectrum(new double[] {1}, new double[] {1}, false), 3, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=2", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 1, nullptr);
		CommonParameters tempVar6();
		Ms2ScanWithSpecificMass *scanB = new Ms2ScanWithSpecificMass(&tempVar5, 2 + 132.040, 1, L"", &tempVar6);
		MsDataScan tempVar7(new MzSpectrum(new double[] {1}, new double[] {1}, false), 4, 1, true, Polarity::Positive, NAN, nullptr, nullptr, MZAnalyzerType::Orbitrap, NAN, nullptr, nullptr, L"scan=3", NAN, nullptr, nullptr, NAN, nullptr, DissociationType::AnyActivationType, 1, nullptr);
		CommonParameters tempVar8();
		Ms2ScanWithSpecificMass *scanC = new Ms2ScanWithSpecificMass(&tempVar7, 3, 1, L"", &tempVar8);

		PeptideSpectralMatch *matchA = new PeptideSpectralMatch(compactPeptide1, 0, 0, 0, scanA, CommonParameters->getDigestionParams(), std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *matchB = new PeptideSpectralMatch(compactPeptide2, 0, 0, 0, scanB, CommonParameters->getDigestionParams(), std::vector<MatchedFragmentIon*>());
		PeptideSpectralMatch *matchC = new PeptideSpectralMatch(compactPeptide3, 0, 0, 0, scanC, CommonParameters->getDigestionParams(), std::vector<MatchedFragmentIon*>());

		auto newPsms = std::vector<PeptideSpectralMatch*> {matchA, matchB, matchC};

		MsDataFile *myMsDataFile = new TestDataFile(std::vector<PeptideWithSetModifications*> {value1.First(), value2.First(), value3.First()});

		auto searchMode = new SinglePpmAroundZeroSearchMode(5);
		std::function<void(const std::vector<PeptideSpectralMatch*>&, const std::wstring&, const std::vector<std::wstring>&)> action2 = [&] (std::vector<PeptideSpectralMatch*> &l, const std::wstring &s, std::vector<std::wstring> &sdf)
		{
			;
		};

		bool DoPrecursorDeconvolution = true;
		bool UseProvidedPrecursorInfo = true;
		double DeconvolutionIntensityRatio = 4;
		int DeconvolutionMaxAssumedChargeState = 10;
		Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);

		CommonParameters tempVar9();
		auto arrayOfMs2ScansSortedByMass = MetaMorpheusTask::GetMs2Scans(myMsDataFile, L"", &tempVar9).OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		std::function<void(BinTreeStructure*, const std::wstring&)> action1 = [&] (BinTreeStructure *l, const std::wstring &s)
		{
			Assert::AreEqual(1, l.FinalBins->Count);
		};

		FdrAnalysisEngine *engine = new FdrAnalysisEngine(newPsms, searchMode->NumNotches, CommonParameters, std::vector<std::wstring> {L"ff"});

		engine->Run();

		delete engine;
		delete DeconvolutionMassTolerance;
		delete searchMode;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
		delete matchC;
		delete matchB;
		delete matchA;
//C# TO C++ CONVERTER TODO TASK: A 'delete scanC' statement was not added since scanC was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scanB' statement was not added since scanB was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scanA' statement was not added since scanA was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}
}
