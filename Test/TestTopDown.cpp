#include "TestTopDown.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../EngineLayer/PeptideSpectralMatch.h"

using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::Indexing;
using namespace EngineLayer::ModernSearch;
using namespace IO::MzML;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{

	void TestTopDown::TestClassicSearchEngineTopDown()
	{
		DigestionParams tempVar(protease: L"top-down");
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1, assumeOrphanPeaksAreZ1Fragments: false);

		auto variableModifications = std::vector<Modification*>();
		auto fixedModifications = std::vector<Modification*>();
		auto proteinList = std::vector<Protein*> {new Protein(L"MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE", L"P40312")};

		auto myMsDataFile = Mzml::LoadAllStaticData(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TopDownTestData\slicedTDYeast.mzML)"));

		auto searchMode = new SinglePpmAroundZeroSearchMode(5);

		bool DoPrecursorDeconvolution = true;
		bool UseProvidedPrecursorInfo = false;
		double DeconvolutionIntensityRatio = 4;
		int DeconvolutionMaxAssumedChargeState = 60;
		Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);

		CommonParameters tempVar2();
		auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, L"", &tempVar2).OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
		ClassicSearchEngine tempVar3(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchMode, CommonParameters, new std::vector<std::wstring>());
		(&tempVar3)->Run();

		auto psm = allPsmsArray.Where([&] (std::any p)
		{
		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchMode' statement was not added since searchMode was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
			return p != nullptr;
		}).FirstOrDefault();
		Assert::That(psm->MatchedFragmentIons->Count > 50);

		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchMode' statement was not added since searchMode was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
	}

	void TestTopDown::TestModernSearchEngineTopDown()
	{
		DigestionParams tempVar(protease: L"top-down");
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1, assumeOrphanPeaksAreZ1Fragments: false);

		auto variableModifications = std::vector<Modification*>();
		auto fixedModifications = std::vector<Modification*>();
		auto proteinList = std::vector<Protein*> {new Protein(L"MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE", L"P40312")};

		auto myMsDataFile = Mzml::LoadAllStaticData(FileSystem::combine(TestContext::CurrentContext->TestDirectory, LR"(TopDownTestData\slicedTDYeast.mzML)"));

		auto searchMode = new SinglePpmAroundZeroSearchMode(5);

		bool DoPrecursorDeconvolution = true;
		bool UseProvidedPrecursorInfo = false;
		double DeconvolutionIntensityRatio = 4;
		int DeconvolutionMaxAssumedChargeState = 60;
		Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);

		CommonParameters tempVar2();
		auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, L"", &tempVar2).OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());

		auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters, 30000, false, std::vector<FileInfo*>(), std::vector<std::wstring>());
		auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());

		ModernSearchEngine tempVar3(allPsmsArray, listOfSortedms2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(), 0, CommonParameters, searchMode, 0, new std::vector<std::wstring>());
		(&tempVar3)->Run();

		auto psm = allPsmsArray.Where([&] (std::any p)
		{
		delete indexEngine;
		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchMode' statement was not added since searchMode was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
			return p != nullptr;
		}).FirstOrDefault();
		Assert::That(psm->MatchedFragmentIons->Count > 50);

		delete indexEngine;
		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchMode' statement was not added since searchMode was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
	}
}
