#include "CoIsolationTests.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../EngineLayer/PeptideSpectralMatch.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;

namespace Test
{

	void CoIsolationTests::TestCoIsolation()
	{
		std::vector<DigestionMotif*> motifs = {new DigestionMotif(L"K", nullptr, 1, nullptr)};
		Protease *protease = new Protease(L"CustProtease", CleavageSpecificity::Full, nullptr, nullptr, motifs);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		DigestionParams tempVar(protease->Name, minPeptideLength: 1);
		CommonParameters *CommonParameters = new CommonParameters(scoreCutoff: 1, deconvolutionIntensityRatio: 50, digestionParams: &tempVar);

		auto variableModifications = std::vector<Modification*>();
		auto fixedModifications = std::vector<Modification*>();
		auto proteinList = std::vector<Protein*> {new Protein(L"MNNNKNDNK", nullptr)};

		auto searchModes = new SinglePpmAroundZeroSearchMode(5);

		Proteomics::AminoAcidPolymer::Peptide *pep1 = new Proteomics::AminoAcidPolymer::Peptide(L"NNNK");
		Proteomics::AminoAcidPolymer::Peptide *pep2 = new Proteomics::AminoAcidPolymer::Peptide(L"NDNK");

		auto dist1 = IsotopicDistribution::GetDistribution(pep1->GetChemicalFormula(), 0.1, 0.01);

		auto dist2 = IsotopicDistribution::GetDistribution(pep2->GetChemicalFormula(), 0.1, 0.01);

		std::vector<MsDataScan*> Scans(2);
		std::vector<double> ms1intensities = {0.8, 0.8, 0.2, 0.02, 0.2, 0.02};
		std::vector<double> ms1mzs = dist1->Masses.Concat(dist2->Masses)->OrderBy([&] (std::any b)
		{
		delete pep2;
		delete pep1;
		delete searchModes;
		delete CommonParameters;
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
			return b;
		})->Select([&] (std::any b)
		{
			b::ToMz(1);
		})->ToArray();

		double selectedIonMz = ms1mzs[1];

		MzSpectrum *MS1 = new MzSpectrum(ms1mzs, ms1intensities, false);

		MzRange tempVar2(300, 2000);
		Scans[0] = new MsDataScan(MS1, 1, 1, false, Polarity::Positive, 1.0, &tempVar2, L"first spectrum", MZAnalyzerType::Unknown, MS1->SumOfAllY, nullptr, nullptr, L"scan=1");

		std::vector<double> ms2intensities = {1, 1, 1, 1, 1};
		std::vector<double> ms2mzs = {146.106.ToMz(1), 228.086.ToMz(1), 229.07.ToMz(1), 260.148.ToMz(1), 342.129.ToMz(1)};
		MzSpectrum *MS2 = new MzSpectrum(ms2mzs, ms2intensities, false);
		double isolationMZ = selectedIonMz;
		MzRange tempVar3(100, 1500);
		Scans[1] = new MsDataScan(MS2, 2, 2, false, Polarity::Positive, 2.0, &tempVar3, L"second spectrum", MZAnalyzerType::Unknown, MS2->SumOfAllY, nullptr, nullptr, L"scan=2", selectedIonMz, nullptr, nullptr, isolationMZ, 2.5, DissociationType::HCD, 1, nullptr);

		auto myMsDataFile = new MsDataFile(Scans, nullptr);

		CommonParameters tempVar4(deconvolutionIntensityRatio: 50);
		auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, L"", &tempVar4).OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
		;
		ClassicSearchEngine tempVar5(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new std::vector<std::wstring>());
		(&tempVar5)->Run();

		// Two matches for this single scan! Corresponding to two co-isolated masses
		Assert::AreEqual(2, allPsmsArray.size());

		Assert::IsTrue(allPsmsArray[0]->getScore() > 1);
		Assert::AreEqual(2, allPsmsArray[0]->getScanNumber());

		Assert::AreEqual(L"NNNK", allPsmsArray[0]->getBaseSequence());
		Assert::AreEqual(L"NDNK", allPsmsArray[1]->getBaseSequence());

//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MS2' statement was not added since MS2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MS1' statement was not added since MS1 was passed to a method or constructor. Handle memory management manually.
		delete pep2;
		delete pep1;
//C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}
}
