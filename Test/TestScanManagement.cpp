#include "TestScanManagement.h"
#include "TestDataFile.h"
#include "../EngineLayer/CommonParameters.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace TaskLayer;

namespace Test
{

	void TestScanManagement::TestGetCombinedMs2Scans()
	{
		auto myMsDataFile = new TestDataFile(5);

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

		//Write prime code to combine MS2MS3
		std::unordered_map<int, double> listOfScanPrecusor;

		std::vector<MsDataScan*> ListOfSortedMsScans;
		std::vector<Ms2ScanWithSpecificMass*> test;

		for (auto ms2scan : myMsDataFile->GetAllScansList().Where([&] (std::any x)
		{
		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
			return x::MsnOrder != 1;
		}))
		{
			if (ms2scan->MsnOrder == 2 && !listOfScanPrecusor.Contains(KeyValuePair<int, double>(ms2scan::OneBasedPrecursorScanNumber->Value, ms2scan::SelectedIonMZ->Value)))
			{
				if (ms2scan::OneBasedPrecursorScanNumber.HasValue)
				{
					listOfScanPrecusor.emplace(ms2scan::OneBasedPrecursorScanNumber->Value, ms2scan::SelectedIonMZ->Value);
					std::vector<int> currentScanMS2OneBasedScanNumber;
					currentScanMS2OneBasedScanNumber.push_back(ms2scan::OneBasedScanNumber);
					auto mz2 = ms2scan::MassSpectrum::XArray::ToList();
					auto intensities2 = ms2scan::MassSpectrum::YArray::ToList();
					for (int i = 1; i < 7; i++)
					{
						if (ms2scan::OneBasedScanNumber + i <= myMsDataFile->NumSpectra)
						{
							auto x = myMsDataFile->GetOneBasedScan(ms2scan::OneBasedScanNumber + i);
							//var x = myMsDataFile.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>().ElementAt(i);

							if (x->MsnOrder == 2 && x->SelectedIonMZ == ms2scan::SelectedIonMZ)
							{
								currentScanMS2OneBasedScanNumber.push_back(x->OneBasedScanNumber);
								mz2.AddRange(x->MassSpectrum.XArray.ToList());
								intensities2.AddRange(x->MassSpectrum.YArray.ToList());
							}
							if (x->MsnOrder == 3 && std::find(currentScanMS2OneBasedScanNumber.begin(), currentScanMS2OneBasedScanNumber.end(), x->OneBasedPrecursorScanNumber->Value) != currentScanMS2OneBasedScanNumber.end())
							{
								mz2.AddRange(x->MassSpectrum.XArray.ToList());
								intensities2.AddRange(x->MassSpectrum.YArray.ToList());
							}
						}
					}
					auto MassSpectrum2 = new MzSpectrum(mz2.ToArray(), intensities2.ToArray(), false);
					MsDataScan tempVar2(MassSpectrum2, ms2scan::OneBasedScanNumber, ms2scan::MsnOrder, ms2scan::IsCentroid, Polarity::Positive, ms2scan::RetentionTime, ms2scan::ScanWindowRange, ms2scan::ScanFilter, ms2scan::MzAnalyzer, ms2scan::TotalIonCurrent, ms2scan::InjectionTime, nullptr, L"", ms2scan::SelectedIonMZ, ms2scan::SelectedIonChargeStateGuess, ms2scan::SelectedIonIntensity, ms2scan::IsolationMz, nullptr, ms2scan::DissociationType, ms2scan::OneBasedPrecursorScanNumber, ms2scan::SelectedIonMonoisotopicGuessMz);
					ListOfSortedMsScans.push_back(&tempVar2);

//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
				}
			}
		}
		for (auto ms2scan : ListOfSortedMsScans.Where([&] (std::any x)
		{
		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
			return x::MsnOrder != 1;
		}))
		{
			Ms2ScanWithSpecificMass tempVar3(ms2scan, ms2scan::SelectedIonMonoisotopicGuessMz->Value, ms2scan::SelectedIonChargeStateGuess->Value, L"", new CommonParameters());
			test.push_back(&tempVar3);
		}
		auto testToArray = test.OrderBy([&] (std::any b)
		{
			b::PrecursorMass;
		})->ToArray();

		//Using function to combine MS2MS3
		//var listOfSortedms2Scans2 = MetaMorpheusTask.GetCombinedMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

		//Assert.AreEqual(5, myMsDataFile.NumSpectra);
		//Assert.AreEqual(1, listOfSortedms2Scans2.Count());

		delete DeconvolutionMassTolerance;
//C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
	}
}
