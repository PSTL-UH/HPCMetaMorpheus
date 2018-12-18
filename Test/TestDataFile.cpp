#include "TestDataFile.h"

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace Test
{

	TestDataFile::TestDataFile() : MsDataFile(nullptr, nullptr, nullptr, nullptr, nullptr)
	{
		auto mz1 = std::vector<double> {50, 60, 70, 80, 90, 402.18629720155.ToMz(2)};
		auto intensities1 = std::vector<double> {1, 1, 1, 1, 1, 1};
		auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
		auto ScansHere = std::vector<MsDataScan*>();
		MzLibUtil::MzRange tempVar2(0, 10000);
		ScansHere.new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=1");

		auto mz2 = std::vector<double> {50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350};
		auto intensities2 = std::vector<double> {1, 1, 1, 1, 1, 1, 1};
		auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
		MsDataScan tempVar3(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), L"f", MZAnalyzerType::Unknown, 100000, 1, nullptr, L"scan=2", 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType::HCD, 1, 402.18629720155.ToMz(2));
		ScansHere.push_back(&tempVar3);

		Scans = ScansHere.ToArray();

//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
	}

	TestDataFile::TestDataFile(double closeMassDifference) : MsDataFile(nullptr, nullptr, nullptr, nullptr, nullptr)
	{
		auto mz1 = std::vector<double> {50, 60, 70, 80, 90, 402.18629720155.ToMz(2)};
		auto intensities1 = std::vector<double> {1, 1, 1, 1, 1, 1};
		auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

		auto ScansHere = std::vector<MsDataScan*>();
		MzLibUtil::MzRange tempVar2(0, 10000);
		ScansHere.new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=1");
		auto mz2 = std::vector<double> {50, 60, 70, 147.0764, 258.132 - closeMassDifference - Constants::ProtonMass, 258.132 - Constants::ProtonMass, 275.1350};
		auto intensities2 = std::vector<double> {1, 1, 1, 1, 1, 1, 1};
		auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
		MsDataScan tempVar3(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), L"f", MZAnalyzerType::Unknown, 100000, 1, nullptr, L"scan=2", 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType::HCD, 1, 402.18629720155.ToMz(2));
		ScansHere.push_back(&tempVar3);

		Scans = ScansHere.ToArray();

//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
	}

	TestDataFile::TestDataFile(bool emptyScan) : MsDataFile(nullptr, nullptr, nullptr, nullptr, nullptr)
	{
		auto mz1 = std::vector<double> {50};
		auto intensities1 = std::vector<double> {1};
		auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

		auto ScansHere = std::vector<MsDataScan*>();
		MzLibUtil::MzRange tempVar2(0, 10000);
		ScansHere.new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=1");
		auto mz2 = std::vector<double> {1};
		auto intensities2 = std::vector<double> {1};
		auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
		MsDataScan tempVar3(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), L"f", MZAnalyzerType::Unknown, 100000, 1, nullptr, L"scan=2", 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType::HCD, 1, 402.18629720155.ToMz(2));
		ScansHere.push_back(&tempVar3);

		Scans = ScansHere.ToArray();

//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
	}

	TestDataFile::TestDataFile(const std::wstring &slightlyLargerDataFile) : MsDataFile(nullptr, nullptr, nullptr, nullptr, nullptr)
	{
		auto mz1 = std::vector<double> {50, 60, 70, 80, 90, 630.27216.ToMz(2)};
		auto intensities1 = std::vector<double> {1, 1, 1, 1, 1, 1};
		auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

		auto ScansHere = std::vector<MsDataScan*>();
		MzLibUtil::MzRange tempVar2(0, 10000);
		ScansHere.new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=1");
		auto mz2 = std::vector<double> {50, 60, 70, 76.0393, 133.0608, 147.0764, 190.0822, 247.1037, 257.1244, 258.127, 275.1350, 385.1830, 442.2045, 630.27216};
		auto intensities2 = std::vector<double> {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
		MsDataScan tempVar3(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), L"f", MZAnalyzerType::Unknown, 100000, 1, nullptr, L"scan=2", 630.27216.ToMz(2), 2, 1, 630.27216.ToMz(2), 2, DissociationType::HCD, 1, 630.27216.ToMz(2));
		ScansHere.push_back(&tempVar3);

		Scans = ScansHere.ToArray();

//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
	}

	TestDataFile::TestDataFile(double precursor, std::vector<double> &products) : MsDataFile(nullptr, nullptr, nullptr, nullptr, nullptr)
	{
		auto mz1 = std::vector<double> {precursor.ToMz(2)};
		auto intensities1 = std::vector<double> {1};
		auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

		auto ScansHere = std::vector<MsDataScan*>();
		MzLibUtil::MzRange tempVar2(0, 10000);
		ScansHere.new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=1");
		auto mz2 = products;
		auto intensities2 = std::vector<double>(products.size());
		for (int i = 0; i < intensities2.size(); i++)
		{
			intensities2[i] = 1;
		}
		auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
		MsDataScan tempVar3(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), L"f", MZAnalyzerType::Unknown, 100000, 1, nullptr, L"scan=2", precursor.ToMz(2), 2, 1, precursor.ToMz(2), 2, DissociationType::HCD, 1, precursor.ToMz(2));
		ScansHere.push_back(&tempVar3);

		Scans = ScansHere.ToArray();

//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
	}

	TestDataFile::TestDataFile(std::vector<PeptideWithSetModifications*> &pepWithSetModss, bool additionalMasses) : MsDataFile(LR"(no nativeID format)", L"mzML format", nullptr, L"SHA-1", LR"(C:\fake.mzML)", nullptr)
	{
		auto ScansHere = std::vector<MsDataScan*>();
		for (int i = 0; i < pepWithSetModss.size(); i++)
		{
			auto pepWithSetMods = pepWithSetModss[i];
			auto mz1 = std::vector<double> {pepWithSetMods->MonoisotopicMass.ToMz(3), (pepWithSetMods->MonoisotopicMass + Constants::C13MinusC12).ToMz(3), (pepWithSetMods->MonoisotopicMass + 2 * Constants::C13MinusC12).ToMz(3), pepWithSetMods->MonoisotopicMass.ToMz(2), (pepWithSetMods->MonoisotopicMass + Constants::C13MinusC12).ToMz(2), (pepWithSetMods->MonoisotopicMass + 2 * Constants::C13MinusC12).ToMz(2)};
			auto intensities1 = std::vector<double> {1, 0.5, 0.25, 1, 0.5, 0.25};
			auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

			MsDataScan tempVar2(MassSpectrum1, 2 * i + 1, 1, true, Polarity::Positive, 2 * i, new MzLibUtil::MzRange(0, 10000), L"gg", MZAnalyzerType::Orbitrap, 1000, 1, nullptr, L"scan=1");
			ScansHere.push_back(&tempVar2);

			std::vector<double> mz2;
			std::vector<double> intensities2;
			std::vector<double> additionalMassesArray;
			if (additionalMasses)
			{
				additionalMassesArray = {260.08307817722, 397.14199003569, 498.18966850487, 612.23259594625, 683.2697097314, 146.10552769922, 217.14264148437};
			}
			else
			{
				additionalMassesArray = std::vector<double>();
			}
			for (auto aok : pepWithSetMods->Fragment(DissociationType::HCD, FragmentationTerminus::Both))
			{
				mz2.push_back(aok->NeutralMass.ToMz(1));
				mz2.push_back((aok->NeutralMass + Constants::C13MinusC12).ToMz(1));
				intensities2.push_back(1);
				intensities2.push_back(1);
			}
			auto MassSpectrum2 = new MzSpectrum(mz2.OrderBy([&] (std::any b)
			{
			delete MassSpectrum2;
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
				return b;
			})->ToArray(), intensities2.ToArray(), false);
			MsDataScan tempVar3(MassSpectrum2, 2 * i + 2, 2, true, Polarity::Positive, 2 * i + 1, new MzLibUtil::MzRange(0, 10000), L"gg", MZAnalyzerType::Orbitrap, 234734, 1, nullptr, L"scan=2", pepWithSetMods->MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods->MonoisotopicMass.ToMz(2), 2, DissociationType::HCD, 2 * i + 1, pepWithSetMods->MonoisotopicMass.ToMz(2));
			ScansHere.push_back(&tempVar3);

//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
		}
		Scans = ScansHere.ToArray();
	}

	TestDataFile::TestDataFile(PeptideWithSetModifications *pepWithSetMods) : MsDataFile(LR"(no nativeID format)", L"mzML format", nullptr, L"SHA-1", LR"(C:\fake.mzML)", nullptr)
	{
		auto mz1 = std::vector<double> {pepWithSetMods->MonoisotopicMass.ToMz(2), (pepWithSetMods->MonoisotopicMass + Constants::C13MinusC12).ToMz(2), (pepWithSetMods->MonoisotopicMass + 2 * Constants::C13MinusC12).ToMz(2)};
		auto intensities1 = std::vector<double> {1, 1, 1};
		auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

		auto ScansHere = std::vector<MsDataScan*>();
		MzLibUtil::MzRange tempVar2(0, 10000);
		ScansHere.new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=1");

		std::vector<double> mz2;
		std::vector<double> intensities2;
		for (auto aok : pepWithSetMods->Fragment(DissociationType::HCD, FragmentationTerminus::Both))
		{
			mz2.push_back(aok->NeutralMass.ToMz(1));
			mz2.push_back((aok->NeutralMass + Constants::C13MinusC12).ToMz(1));
			intensities2.push_back(1);
			intensities2.push_back(1);
		}
		auto MassSpectrum2 = new MzSpectrum(mz2.OrderBy([&] (std::any b)
		{
		delete MassSpectrum2;
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
			return b;
		})->ToArray(), intensities2.ToArray(), false);
		MzLibUtil::MzRange tempVar3(0, 10000);
		auto scan2 = new MsDataScan(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, &tempVar3, L"df", MZAnalyzerType::Orbitrap, 234734, 1, nullptr, L"scan=2", pepWithSetMods->MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods->MonoisotopicMass.ToMz(2), 2, DissociationType::HCD, 1, pepWithSetMods->MonoisotopicMass.ToMz(2));
		scan2->ComputeSelectedPeakIntensity(MassSpectrum1);
		scan2->ComputeMonoisotopicPeakIntensity(MassSpectrum1);
		ScansHere.push_back(scan2);
		Scans = ScansHere.ToArray();

//C# TO C++ CONVERTER TODO TASK: A 'delete scan2' statement was not added since scan2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
	}

	TestDataFile::TestDataFile(PeptideWithSetModifications *pepWithSetMods, const std::wstring &v) : MsDataFile(nullptr, nullptr, nullptr, nullptr, nullptr)
	{
		if (v == L"quadratic")
		{
			// Add three ms1 peaks with charge 2, exact
			auto MassSpectrum1 = new MzSpectrum(std::vector<double> {pepWithSetMods->MonoisotopicMass.ToMz(2), (pepWithSetMods->MonoisotopicMass + Constants::C13MinusC12).ToMz(2), (pepWithSetMods->MonoisotopicMass + 2 * Constants::C13MinusC12).ToMz(2)}, std::vector<double> {1, 1, 1}, false);

			std::vector<double> mz2;
			std::vector<double> intensities2;
			for (auto aok : pepWithSetMods->Fragment(DissociationType::HCD, FragmentationTerminus::Both))
			{
				auto t1 = aok->NeutralMass.ToMz(1);
				auto c = 0.0000001;
				mz2.push_back(t1 + c * std::pow(t1, 2));
				auto t2 = (aok->NeutralMass + Constants::C13MinusC12).ToMz(1);
				mz2.push_back(t2 + c * std::pow(t2, 2));
				intensities2.push_back(1);
				intensities2.push_back(1);
			}
			auto MassSpectrum2 = new MzSpectrum(mz2.OrderBy([&] (std::any b)
			{
			delete MassSpectrum2;
			delete MassSpectrum1;
				return b;
			})->ToArray(), intensities2.ToArray(), false);

			MzLibUtil::MzRange tempVar2(0, 10000);
			auto scan2 = new MsDataScan(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, &tempVar2, L"df", MZAnalyzerType::Orbitrap, 234734, 1, nullptr, L"scan=2", pepWithSetMods->MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods->MonoisotopicMass.ToMz(2), 2, DissociationType::HCD, 1, pepWithSetMods->MonoisotopicMass.ToMz(2));
			scan2->ComputeSelectedPeakIntensity(MassSpectrum1);
			scan2->ComputeMonoisotopicPeakIntensity(MassSpectrum1);
			MzLibUtil::MzRange tempVar3(0, 10000);
			auto ScansHere = std::vector<MsDataScan*>
			{
				new MsDataScan(MassSpectrum1,1, 1, true, Polarity::Positive, 1, &tempVar3, L"ff", MZAnalyzerType::Unknown, 1000,1, nullptr, L"scan=1"),
				scan2
			};
			Scans = ScansHere.ToArray();

			delete scan2;
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
		}
	}

	TestDataFile::TestDataFile(PeptideWithSetModifications *pepWithSetMods, int charge, double intensity, double rt) : MsDataFile(nullptr, nullptr, nullptr, nullptr, nullptr)
	{
		auto mz1 = std::vector<double> {pepWithSetMods->MonoisotopicMass.ToMz(charge), (pepWithSetMods->MonoisotopicMass + 1.003).ToMz(charge), (pepWithSetMods->MonoisotopicMass + 2.005).ToMz(charge)};
		auto intensities1 = std::vector<double> {intensity, intensity * 10, intensity / 10};
		auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

		auto ScansHere = std::vector<MsDataScan*>();
		MzLibUtil::MzRange tempVar2(0, 10000);
		ScansHere.new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, rt, &tempVar2, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=1");

		std::vector<double> mz2;
		std::vector<double> intensities2;
		for (auto aok : pepWithSetMods->Fragment(DissociationType::HCD,FragmentationTerminus::Both))
		{
			mz2.push_back(aok->NeutralMass.ToMz(1));
			mz2.push_back((aok->NeutralMass + 1.003).ToMz(1));
			intensities2.push_back(intensity);
			intensities2.push_back(intensity);
		}
		auto MassSpectrum2 = new MzSpectrum(mz2.OrderBy([&] (std::any b)
		{
		delete MassSpectrum2;
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
			return b;
		})->ToArray(), intensities2.ToArray(), false);
		MzLibUtil::MzRange tempVar3(0, 10000);
		auto scan2 = new MsDataScan(MassSpectrum2, 2, 2, true, Polarity::Positive, rt + 0.01, &tempVar3, L"df", MZAnalyzerType::Orbitrap, 234734, 1, nullptr, L"scan=2", pepWithSetMods->MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods->MonoisotopicMass.ToMz(2), 2, DissociationType::HCD, 1, pepWithSetMods->MonoisotopicMass.ToMz(2));
		scan2->ComputeSelectedPeakIntensity(MassSpectrum1);
		scan2->ComputeMonoisotopicPeakIntensity(MassSpectrum1);
		ScansHere.push_back(scan2);
		Scans = ScansHere.ToArray();

//C# TO C++ CONVERTER TODO TASK: A 'delete scan2' statement was not added since scan2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
	}

	TestDataFile::TestDataFile(int MS3) : MsDataFile(nullptr, nullptr, nullptr, nullptr, nullptr)
	{
		auto mz1 = std::vector<double> {50, 60, 70, 80, 90, 764.1376.ToMz(2)};
		auto intensities1 = std::vector<double> {1, 1, 1, 1, 1, 1};
		auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
		auto ScansHere = std::vector<MsDataScan*>();
		MzLibUtil::MzRange tempVar2(0, 10000);
		ScansHere.new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=1");

		auto mz2 = std::vector<double> {52, 62, 72, 147.0764, 257.1244, 258.127, 275.1350, 502};
		auto intensities2 = std::vector<double> {1, 1, 1, 1, 1, 1, 1, 1};
		auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
		MsDataScan tempVar3(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), L"f", MZAnalyzerType::Unknown, 100000, 1, nullptr, L"scan=2", 764.1376.ToMz(2), 2, 1, 764.1376.ToMz(2), 2, DissociationType::CID, 1, 764.1376.ToMz(1));
		ScansHere.push_back(&tempVar3);

		auto mz3 = std::vector<double> {53, 63, 73, 148.0764, 258.1244, 259.127, 276.1350, 503};
		auto intensities3 = std::vector<double> {1, 1, 1, 1, 1, 1, 1, 1};
		auto MassSpectrum3 = new MzSpectrum(mz3, intensities3, false);
		MsDataScan tempVar4(MassSpectrum3, 3, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), L"f", MZAnalyzerType::Unknown, 100000, 1, nullptr, L"scan=3", 764.1376.ToMz(2), 2, 1, 764.1376.ToMz(2), 2, DissociationType::ETD, 1, 764.1376.ToMz(1));
		ScansHere.push_back(&tempVar4);

		auto mz4 = std::vector<double> {54, 64, 74, 149.0764, 259.1244, 260.127, 277.1350, 504};
		auto intensities4 = std::vector<double> {1, 1, 1, 1, 1, 1, 1, 1};
		auto MassSpectrum4 = new MzSpectrum(mz4, intensities4, false);
		MsDataScan tempVar5(MassSpectrum4, 4, 3, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), L"f", MZAnalyzerType::Unknown, 100000, 1, nullptr, L"scan=4", 275.1350.ToMz(1), 1, 1, 275.1350.ToMz(1), 1, DissociationType::HCD, 2, 275.1350.ToMz(1));
		ScansHere.push_back(&tempVar5);

		auto mz5 = std::vector<double> {55, 65, 75, 150.0764, 260.1244, 261.127, 278.1350, 505};
		auto intensities5 = std::vector<double> {1, 1, 1, 1, 1, 1, 1, 1};
		auto MassSpectrum5 = new MzSpectrum(mz5, intensities5, false);
		MsDataScan tempVar6(MassSpectrum5, 5, 3, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), L"f", MZAnalyzerType::Unknown, 100000, 1, nullptr, L"scan=5", 257.1244.ToMz(1), 1, 1, 257.1244.ToMz(1), 1, DissociationType::HCD, 2, 257.1244.ToMz(1));
		ScansHere.push_back(&tempVar6);

		Scans = ScansHere.ToArray();

//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum5' statement was not added since MassSpectrum5 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum4' statement was not added since MassSpectrum4 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum3' statement was not added since MassSpectrum3 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
	}

	TestDataFile::TestDataFile(std::vector<double> &ms2Mz, std::vector<double> &ms2Intensities, double precursorMass, int precursorZ, double rt) : MsDataFile(nullptr, nullptr, nullptr, nullptr, nullptr)
	{
		auto ms1 = new MzSpectrum(std::vector<double> {precursorMass.ToMz(precursorZ), (precursorMass + 1.003).ToMz(precursorZ)}, std::vector<double> {1, 1}, false);
		auto ms2 = new MzSpectrum(ms2Mz, ms2Intensities, false);

		auto ScansHere = std::vector<MsDataScan*>();
		MzLibUtil::MzRange tempVar2(0, 10000);
		ScansHere.new MsDataScan(ms1, 1, 1, true, Polarity::Positive, rt, &tempVar2, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=1");
		MzLibUtil::MzRange tempVar3(0, 10000);
		ScansHere.new MsDataScan(ms2, 1, 2, true, Polarity::Positive, rt + 0.01, &tempVar3, L"ff", MZAnalyzerType::Unknown, 1000, 1, nullptr, L"scan=2", precursorMass.ToMz(precursorZ), precursorZ, 1, precursorMass.ToMz(precursorZ), 1.0, DissociationType::HCD, 1, precursorMass.ToMz(precursorZ));

		Scans = ScansHere.ToArray();

//C# TO C++ CONVERTER TODO TASK: A 'delete ms2' statement was not added since ms2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete ms1' statement was not added since ms1 was passed to a method or constructor. Handle memory management manually.
	}

	std::wstring TestDataFile::getFilePath() const
	{
		return L"TestDataFile";
	}

	std::wstring TestDataFile::getName() const
	{
		return L"TestDataFile";
	}

	void TestDataFile::ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities)
	{
		MzSpectrum *massSpectrum = new MzSpectrum(mz, intensities, false);
		Scans[0] = new MsDataScan(massSpectrum, Scans[0].OneBasedScanNumber, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum->SumOfAllY, Scans[0].InjectionTime, nullptr, Scans[0].NativeId);

//C# TO C++ CONVERTER TODO TASK: A 'delete massSpectrum' statement was not added since massSpectrum was passed to a method or constructor. Handle memory management manually.
	}

	MsDataScan *TestDataFile::GetOneBasedScan(int scanNumber)
	{
		return Scans[scanNumber - 1];
	}

	std::vector<MsDataScan*> TestDataFile::GetMS1Scans()
	{
		throw NotImplementedException();
	}
}
