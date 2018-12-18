#include "MyPeptideTest.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/PrecursorSearchModes/OpenMassDiffAcceptor.h"

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
using namespace UsefulProteomicsDatabases;

namespace Test
{

	void MyPeptideTest::TestIdenticalPeaks()
	{
		std::unordered_map<int, std::vector<Modification*>> mods;
		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"M", motif);
		mods.emplace(1, std::vector<Modification*> {new Modification(_originalId: L"Hehe", _target: motif, _locationRestriction: L"Anywhere.", _monoisotopicMass: 18.010565)});
		auto prot = new Protein(L"MMMM", nullptr, nullptr, nullptr, mods);
		DigestionParams *digestionParams = new DigestionParams(minPeptideLength: 1);
		auto ye = prot->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()).First();

		auto massArray = ye->Fragment(DissociationType::HCD, FragmentationTerminus::Both)->Select([&] (std::any p)
		{
			p::NeutralMass;
		})->ToArray();
		std::sort(massArray.begin(), massArray.end());
		std::vector<double> intensities = {1, 1, 1, 1};
		std::vector<double> mz = {massArray[0].ToMz(1), massArray[2].ToMz(1), massArray[4].ToMz(1), 10000};
		MzSpectrum *massSpectrum = new MzSpectrum(mz, intensities, false);
		MzRange tempVar(300, 2000);
		MsDataScan *scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity::Positive, 1, &tempVar, L"", MZAnalyzerType::Unknown, massSpectrum->SumOfAllY, nullptr, nullptr, L"scan=1", 0, nullptr, nullptr, 0, nullptr, DissociationType::Unknown, 1, nullptr);

		std::vector<PeptideSpectralMatch*> globalPsms(1);
		CommonParameters tempVar2();
		std::vector<Ms2ScanWithSpecificMass*> arrayOfSortedMS2Scans = {new Ms2ScanWithSpecificMass(scan, 0, 1, L"", &tempVar2)};
		PpmTolerance tempVar3(5);
		DigestionParams tempVar4(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain);
		CommonParameters *CommonParameters = new CommonParameters(productMassTolerance: &tempVar3, scoreCutoff: 1, digestionParams: &tempVar4);
		OpenSearchMode tempVar5();
		ClassicSearchEngine *cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, std::vector<Modification*>(), std::vector<Modification*>(), {prot}, &tempVar5, CommonParameters, std::vector<std::wstring>());

		cse->Run();
		Assert::AreEqual(3, globalPsms[0]->MatchedFragmentIons->Count);

		delete cse;
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete massSpectrum' statement was not added since massSpectrum was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete prot' statement was not added since prot was passed to a method or constructor. Handle memory management manually.
	}

	void MyPeptideTest::TestLastPeaks()
	{
		std::unordered_map<int, std::vector<Modification*>> mods;
		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"M", motif);
		auto prot = new Protein(L"MMMM", nullptr, nullptr, nullptr, mods);
		DigestionParams *digestionParams = new DigestionParams(minPeptideLength: 1);
		PeptideWithSetModifications *thePep = prot->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()).First();

		auto massArray = thePep->Fragment(DissociationType::HCD, FragmentationTerminus::Both)->Select([&] (std::any p)
		{
			p::NeutralMass;
		})->ToArray();
		std::sort(massArray.begin(), massArray.end());
		std::vector<double> intensities = {1, 1, 1};
		std::vector<double> mz = {1, 2, massArray[4].ToMz(1)};
		MzSpectrum *massSpectrum = new MzSpectrum(mz, intensities, false);
		MzRange tempVar(300, 2000);
		MsDataScan *scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity::Positive, 1, &tempVar, L"", MZAnalyzerType::Unknown, massSpectrum->SumOfAllY, nullptr, nullptr, L"scan=1", 0, nullptr, nullptr, 0, nullptr, DissociationType::Unknown, 1, nullptr);

		std::vector<PeptideSpectralMatch*> globalPsms(1);
		CommonParameters tempVar2();
		std::vector<Ms2ScanWithSpecificMass*> arrayOfSortedMS2Scans = {new Ms2ScanWithSpecificMass(scan, 0, 1, L"", &tempVar2)};
		PpmTolerance tempVar3(5);
		DigestionParams tempVar4(maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: std::numeric_limits<int>::max(), initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain);
		CommonParameters *CommonParameters = new CommonParameters(scoreCutoff: 1, productMassTolerance: &tempVar3, digestionParams: &tempVar4);
		OpenSearchMode tempVar5();
		ClassicSearchEngine *cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, std::vector<Modification*>(), std::vector<Modification*>(), {prot}, &tempVar5, CommonParameters, std::vector<std::wstring>());

		cse->Run();
		Assert::Less(globalPsms[0]->getScore(), 2);
		Assert::Greater(globalPsms[0]->getScore(), 1);

		delete cse;
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete massSpectrum' statement was not added since massSpectrum was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete prot' statement was not added since prot was passed to a method or constructor. Handle memory management manually.
	}

	void MyPeptideTest::TestVeryCloseExperimentalsClassic()
	{
		std::unordered_map<int, std::vector<Modification*>> mods;
		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"M", motif);
		auto prot = new Protein(L"MMMM", nullptr, nullptr, nullptr, mods);
		DigestionParams *digestionParams = new DigestionParams(minPeptideLength: 1);
		auto thePep = prot->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()).First();

		auto massArray = thePep->Fragment(DissociationType::HCD, FragmentationTerminus::Both)->Select([&] (std::any p)
		{
			p::NeutralMass;
		})->ToArray();
		std::sort(massArray.begin(), massArray.end());
		std::vector<double> intensities = {1, 1, 1, 1};
		std::vector<double> mz = {1, 2, massArray[4].ToMz(1), massArray[4].ToMz(1) + 1e-9};
		MzSpectrum *massSpectrum = new MzSpectrum(mz, intensities, false);
		MzRange tempVar(300, 2000);
		MsDataScan *scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity::Positive, 1, &tempVar, L"", MZAnalyzerType::Unknown, massSpectrum->SumOfAllY, nullptr, nullptr, L"scan=1", 0, nullptr, nullptr, 0, nullptr, DissociationType::Unknown, 1, nullptr);

		std::vector<PeptideSpectralMatch*> globalPsms(1);
		CommonParameters tempVar2();
		std::vector<Ms2ScanWithSpecificMass*> arrayOfSortedMS2Scans = {new Ms2ScanWithSpecificMass(scan, 0, 1, L"", &tempVar2)};
		PpmTolerance tempVar3(5);
		DigestionParams tempVar4(maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: std::numeric_limits<int>::max(), initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain);
		CommonParameters *CommonParameters = new CommonParameters(productMassTolerance: &tempVar3, scoreCutoff: 1, digestionParams: &tempVar4);
		OpenSearchMode tempVar5();
		ClassicSearchEngine *cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, std::vector<Modification*>(), std::vector<Modification*>(), {prot}, &tempVar5, CommonParameters, std::vector<std::wstring>());

		cse->Run();
		Assert::Less(globalPsms[0]->getScore(), 2);
		Assert::Greater(globalPsms[0]->getScore(), 1);

		delete cse;
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete massSpectrum' statement was not added since massSpectrum was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete prot' statement was not added since prot was passed to a method or constructor. Handle memory management manually.
	}

	void MyPeptideTest::TestVeryCloseExperimentalsModern()
	{
		std::unordered_map<int, std::vector<Modification*>> mods;
		ModificationMotif motif;
		ModificationMotif::TryGetMotif(L"M", motif);
		auto prot = new Protein(L"MMMM", nullptr, nullptr, nullptr, mods);
		DigestionParams *digestionParams = new DigestionParams(minPeptideLength: 1);
		auto thePep = prot->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()).First();

		auto massArray = thePep->Fragment(DissociationType::HCD, FragmentationTerminus::Both)->Select([&] (std::any p)
		{
			p::NeutralMass;
		})->ToArray();
		std::sort(massArray.begin(), massArray.end());
		std::vector<double> intensities = {1, 1, 1, 1};
		std::vector<double> mz = {1, 2, massArray[4].ToMz(1), massArray[4].ToMz(1) + 1e-9};
		MzSpectrum *massSpectrum = new MzSpectrum(mz, intensities, false);
		MzRange tempVar(300, 2000);
		MsDataScan *scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity::Positive, 1, &tempVar, L"", MZAnalyzerType::Unknown, massSpectrum->SumOfAllY, nullptr, nullptr, L"scan=1", 0, nullptr, nullptr, 0, nullptr, DissociationType::Unknown, 1, nullptr);

		std::vector<PeptideSpectralMatch*> globalPsms(1);
		CommonParameters tempVar2();
		std::vector<Ms2ScanWithSpecificMass*> arrayOfSortedMS2Scans = {new Ms2ScanWithSpecificMass(scan, 600, 1, L"", &tempVar2)};
		PpmTolerance tempVar3(5);
		DigestionParams tempVar4(maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: std::numeric_limits<int>::max(), initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain);
		CommonParameters *CommonParameters = new CommonParameters(productMassTolerance: &tempVar3, scoreCutoff: 1, digestionParams: &tempVar4);

		auto indexEngine = new IndexingEngine(std::vector<Protein*> {prot}, std::vector<Modification*>(), std::vector<Modification*>(), 1, DecoyType::Reverse, CommonParameters, 30000, false, std::vector<FileInfo*>(), std::vector<std::wstring>());
		auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
		OpenSearchMode tempVar5();
		auto cse = new ModernSearchEngine(globalPsms, arrayOfSortedMS2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(), 0, CommonParameters, &tempVar5, 0, std::vector<std::wstring>());

		cse->Run();
		Assert::Less(globalPsms[0]->getScore(), 2);
		Assert::Greater(globalPsms[0]->getScore(), 1);

		delete cse;
		delete indexEngine;
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete massSpectrum' statement was not added since massSpectrum was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete prot' statement was not added since prot was passed to a method or constructor. Handle memory management manually.
	}

	void MyPeptideTest::TestAllNaN()
	{
		std::unordered_map<int, std::vector<Modification*>> mods;
		auto prot = new Protein(L"XMMM", nullptr, nullptr, nullptr, mods);
		DigestionParams *digestionParams = new DigestionParams(minPeptideLength: 1);
		auto thePep = prot->Digest(digestionParams, std::vector<Modification*>(), std::vector<Modification*>()).First();

		auto massArray = thePep->Fragment(DissociationType::HCD, FragmentationTerminus::Both)->Select([&] (std::any p)
		{
			p::NeutralMass;
		})->ToArray();
		std::sort(massArray.begin(), massArray.end());
		std::vector<double> intensities = {1, 1, 1, 1};
		std::vector<double> mz = {1, 2, 3, 4};
		MzSpectrum *massSpectrum = new MzSpectrum(mz, intensities, false);
		MzRange tempVar(300, 2000);
		MsDataScan *scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity::Positive, 1, &tempVar, L"", MZAnalyzerType::Unknown, massSpectrum->SumOfAllY, nullptr, nullptr, L"scan=1", 0, nullptr, nullptr, 0, nullptr, DissociationType::Unknown, 1, nullptr);

		std::vector<PeptideSpectralMatch*> globalPsms(1);
		CommonParameters tempVar2();
		std::vector<Ms2ScanWithSpecificMass*> arrayOfSortedMS2Scans = {new Ms2ScanWithSpecificMass(scan, 0, 0, L"", &tempVar2)};
		PpmTolerance tempVar3(5);
		DigestionParams tempVar4(maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: std::numeric_limits<int>::max(), initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain);
		CommonParameters *CommonParameters = new CommonParameters(productMassTolerance: &tempVar3, scoreCutoff: 1, digestionParams: &tempVar4);

		OpenSearchMode tempVar5();
		ClassicSearchEngine *cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, std::vector<Modification*>(), std::vector<Modification*>(), {prot}, &tempVar5, CommonParameters, std::vector<std::wstring>());

		cse->Run();
		Assert::IsNull(globalPsms[0]);

		delete cse;
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete massSpectrum' statement was not added since massSpectrum was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete prot' statement was not added since prot was passed to a method or constructor. Handle memory management manually.
	}
}
