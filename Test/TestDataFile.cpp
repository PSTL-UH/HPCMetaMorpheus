#include "TestDataFile.h"

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace Test
{
    
    TestDataFile::TestDataFile()  : MsDataFile(2,  new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        
        //auto mz1 = std::vector<double> {50, 60, 70, 80, 90, 402.18629720155.ToMz(2)};
        auto mz1 = std::vector<double> {50, 60, 70, 80, 90, Chemistry::ClassExtensions::ToMz(402.18629720155,2)};
        auto intensities1 = std::vector<double> {1, 1, 1, 1, 1, 1};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

        MzLibUtil::MzRange tempVar2(0, 10000);
        auto ScansHere = new std::vector<MsDataScan*>();
        auto s = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, "ff",
                                MZAnalyzerType::Unknown, 1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(s);
        
        auto mz2 = std::vector<double> {50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350};
        auto intensities2 = std::vector<double> {1, 1, 1, 1, 1, 1, 1};
        auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
        auto tempVar3 = new MsDataScan (MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), "f",
                                        MZAnalyzerType::Unknown, 100000, 1, std::vector<std::vector<double>>(), "scan=2",
                                        Chemistry::ClassExtensions::ToMz(402.18629720155,2),
                                        2, 1,
                                        Chemistry::ClassExtensions::ToMz(402.18629720155, 2), 2, DissociationType::HCD, 1,
                                        Chemistry::ClassExtensions::ToMz(402.18629720155, 2));
        ScansHere->push_back(tempVar3);
        
        Scans = *ScansHere; //.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    TestDataFile::TestDataFile(double closeMassDifference) :
        MsDataFile(2,  new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        auto mz1 = std::vector<double> {50, 60, 70, 80, 90, Chemistry::ClassExtensions::ToMz(402.18629720155,2)};
        auto intensities1 = std::vector<double> {1, 1, 1, 1, 1, 1};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
        
        MzLibUtil::MzRange tempVar2(0, 10000);
        auto ScansHere = new std::vector<MsDataScan*>();
        auto S = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, "ff",
                               MZAnalyzerType::Unknown, 1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(S);
        std::vector<double> mz2 = {50, 60, 70, 147.0764, 258.132 - closeMassDifference - Constants::protonMass,
                                   258.132 - Constants::protonMass, 275.1350};
        std::vector<double> intensities2 =  {1, 1, 1, 1, 1, 1, 1};
        auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
        auto tempVar3 =  new MsDataScan (MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), "f",
                                         MZAnalyzerType::Unknown, 100000, 1, std::vector<std::vector<double>>(), "scan=2",
                                         Chemistry::ClassExtensions::ToMz(402.18629720155, 2), 2, 1,
                                         Chemistry::ClassExtensions::ToMz(402.18629720155,2), 2, DissociationType::HCD, 1,
                                         Chemistry::ClassExtensions::ToMz(402.18629720155, 2));
        ScansHere->push_back(tempVar3);
        
        Scans = *ScansHere; //.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    TestDataFile::TestDataFile(bool emptyScan) :
        MsDataFile( 2,  new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        std::vector<double> mz1 = {50};
        std::vector<double> intensities1 =  {1};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
        
        MzLibUtil::MzRange tempVar2(0, 10000);
        auto ScansHere = new std::vector<MsDataScan*>();
        auto S  = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, "ff",
                                 MZAnalyzerType::Unknown, 1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(S);
        
        std::vector<double> mz2 = {1};
        std::vector<double> intensities2 =  {1};
        auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
        auto tempVar3 = new MsDataScan (MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), "f",
                                        MZAnalyzerType::Unknown, 100000, 1, std::vector<std::vector<double>>(), "scan=2",
                                        Chemistry::ClassExtensions::ToMz(402.18629720155, 2),
                                        2, 1,
                                        Chemistry::ClassExtensions::ToMz(402.18629720155, 2), 2, DissociationType::HCD, 1,
                                        Chemistry::ClassExtensions::ToMz(402.18629720155, 2));
        ScansHere->push_back(tempVar3);
        
        Scans = *ScansHere; //.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    TestDataFile::TestDataFile(const std::string &slightlyLargerDataFile) :
        MsDataFile(2,  new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        std::vector<double> mz1 = {50, 60, 70, 80, 90, Chemistry::ClassExtensions::ToMz(630.27216, 2)};
        std::vector<double> intensities1 =  {1, 1, 1, 1, 1, 1};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
        
        MzLibUtil::MzRange tempVar2(0, 10000);
        auto ScansHere = new std::vector<MsDataScan*>();
        auto S = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, "ff",
                                MZAnalyzerType::Unknown, 1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(S);
        
        std::vector<double> mz2 = {50, 60, 70, 76.0393, 133.0608, 147.0764, 190.0822, 247.1037, 257.1244,
                                   258.127, 275.1350, 385.1830, 442.2045, 630.27216};
        std::vector<double> intensities2 =  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
        auto tempVar3 = new MsDataScan (MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000), "f",
                                        MZAnalyzerType::Unknown, 100000, 1, std::vector<std::vector<double>>(), "scan=2",
                                        Chemistry::ClassExtensions::ToMz(630.27216, 2), 2, 1,
                                        Chemistry::ClassExtensions::ToMz(630.27216, 2),
                                        2, DissociationType::HCD, 1,
                                        Chemistry::ClassExtensions::ToMz(630.27216, 2));
        ScansHere->push_back(tempVar3);
        
        Scans = *ScansHere; //.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    TestDataFile::TestDataFile(double precursor, std::vector<double> &products) :
        MsDataFile(2,  new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        std::vector<double> mz1 = {Chemistry::ClassExtensions::ToMz(precursor, 2)};
        std::vector<double> intensities1 =  {1};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
        
        MzLibUtil::MzRange tempVar2(0, 10000);
        auto ScansHere = new std::vector<MsDataScan*>();
        auto S = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, "ff",
                                MZAnalyzerType::Unknown, 1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(S);
        
        auto mz2 = products;
        std::vector<double> intensities2(products.size());
        for (int i = 0; i < (int)intensities2.size(); i++)
        {
            intensities2[i] = 1;
        }
        auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
        auto tempVar3 = new MsDataScan (MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000),
                                        "f", MZAnalyzerType::Unknown, 100000, 1, std::vector<std::vector<double>>(), "scan=2",
                                        Chemistry::ClassExtensions::ToMz(precursor, 2), 2, 1,
                                        Chemistry::ClassExtensions::ToMz(precursor, 2), 2, DissociationType::HCD, 1,
                                        Chemistry::ClassExtensions::ToMz(precursor, 2));
        ScansHere->push_back(tempVar3);
        
        Scans = *ScansHere; //.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    TestDataFile::TestDataFile(std::vector<PeptideWithSetModifications*> &pepWithSetModss, bool additionalMasses) :
        MsDataFile(2,  new MassSpectrometry::SourceFile("(no nativeID format)", "mzML format", "", "SHA-1", "fake.mzML", ""))
    {
        auto ScansHere = new std::vector<MsDataScan*>();
        for (int i = 0; i < (int)pepWithSetModss.size(); i++)
        {
            auto pepWithSetMods = pepWithSetModss[i];
            std::vector<double> mz1 = {Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 3),
                                       Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + Constants::C13MinusC12, 3),
                                       Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + 2 * Constants::C13MinusC12,3),
                                       Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2),
                                       Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + Constants::C13MinusC12, 2),
                                       Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + 2 * Constants::C13MinusC12, 2)};
            std::vector<double> intensities1 =  {1, 0.5, 0.25, 1, 0.5, 0.25};
            auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
            auto tempVar2 = new MsDataScan (MassSpectrum1, 2 * i + 1, 1, true, Polarity::Positive, 2 * i,
                                            new MzLibUtil::MzRange(0, 10000), "gg", MZAnalyzerType::Orbitrap,
                                            1000, 1, std::vector<std::vector<double>>(), "scan=1");
            ScansHere->push_back(tempVar2);
            
            std::vector<double> mz2;
            std::vector<double> intensities2;
            std::vector<double> additionalMassesArray;
            if (additionalMasses)
            {
                additionalMassesArray = {260.08307817722, 397.14199003569, 498.18966850487,
                                         612.23259594625, 683.2697097314, 146.10552769922, 217.14264148437};
            }

            for (auto aok : pepWithSetMods->Fragment(DissociationType::HCD, FragmentationTerminus::Both))
            {
                mz2.push_back(Chemistry::ClassExtensions::ToMz(aok->NeutralMass,  1));
                mz2.push_back(Chemistry::ClassExtensions::ToMz(aok->NeutralMass + Constants::C13MinusC12, 1));
                intensities2.push_back(1);
                intensities2.push_back(1);
            }
#ifdef ORIG
            auto MassSpectrum2 = new MzSpectrum(mz2.OrderBy([&] (std::any b)
            {
                return b;
            })->ToArray(), intensities2.ToArray(), false);
#endif
            std::sort(mz2.begin(), mz2.end() );            
            auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false );

            auto tempVar3 = new MsDataScan (MassSpectrum2, 2 * i + 2, 2, true, Polarity::Positive, 2 * i + 1,
                                            new MzLibUtil::MzRange(0, 10000), "gg", MZAnalyzerType::Orbitrap, 234734, 1,
                                            std::vector<std::vector<double>>(), "scan=2",
                                            Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2), 2, 1,
                                            Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2),
                                            2, DissociationType::HCD, 2 * i + 1,
                                            Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2));
            ScansHere->push_back(tempVar3);
            
            //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
            //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
            //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
            //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
        }
        Scans = *ScansHere; //.ToArray();
    }
    
    TestDataFile::TestDataFile(PeptideWithSetModifications *pepWithSetMods) :
        MsDataFile(2,  new MassSpectrometry::SourceFile("(no nativeID format)", "mzML format", "", "SHA-1", "fake.mzML", ""))
    {
        std::vector<double> mz1 = {Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2),
                                   Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + Constants::C13MinusC12, 2),
                                   Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + 2 * Constants::C13MinusC12, 2)};
        std::vector<double> intensities1 =  {1, 1, 1};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
        
        MzLibUtil::MzRange tempVar2(0, 10000);
        auto ScansHere = new std::vector<MsDataScan*>();
        auto S = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, "ff",
                                MZAnalyzerType::Unknown, 1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(S);
        
        std::vector<double> mz2;
        std::vector<double> intensities2;
        for (auto aok : pepWithSetMods->Fragment(DissociationType::HCD, FragmentationTerminus::Both))
        {
            mz2.push_back(Chemistry::ClassExtensions::ToMz(aok->NeutralMass, 1));
            mz2.push_back(Chemistry::ClassExtensions::ToMz(aok->NeutralMass + Constants::C13MinusC12, 1));
            intensities2.push_back(1);
            intensities2.push_back(1);
        }

#ifdef ORIG
        auto MassSpectrum2 = new MzSpectrum(mz2.OrderBy([&] (std::any b)
        {
            return b;
        })->ToArray(), intensities2.ToArray(), false);
#endif
        std::sort(mz2.begin(), mz2.end() );
        auto MassSpectrum2 = new MzSpectrum ( mz2, intensities2, false );
        
        MzLibUtil::MzRange tempVar3(0, 10000);
        auto scan2 = new MsDataScan(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, &tempVar3, "df",
                                    MZAnalyzerType::Orbitrap, 234734, 1, std::vector<std::vector<double>>(), "scan=2",
                                    Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2), 2, 1,
                                    Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2), 2, DissociationType::HCD, 1,
                                    Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2));
        scan2->ComputeSelectedPeakIntensity(MassSpectrum1);
        scan2->ComputeMonoisotopicPeakIntensity(MassSpectrum1);
        ScansHere->push_back(scan2);

        Scans = *ScansHere; //.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete scan2' statement was not added since scan2 was
        //passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    TestDataFile::TestDataFile(PeptideWithSetModifications *pepWithSetMods, const std::string &v) :
        MsDataFile(2,  new MassSpectrometry::SourceFile( "", "", "", "", "" ))
    {
        if (v == "quadratic")
        {
            // Add three ms1 peaks with charge 2, exact
            std::vector<double> v1 = {Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2),
                                      Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + Constants::C13MinusC12, 2),
                                      Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + 2 * Constants::C13MinusC12, 2)};
            std::vector<double> v2 = {1, 1, 1};
            auto MassSpectrum1 = new MzSpectrum(v1 , v2, false);
            
            std::vector<double> mz2;
            std::vector<double> intensities2;
            for (auto aok : pepWithSetMods->Fragment(DissociationType::HCD, FragmentationTerminus::Both))
            {
                auto t1 = Chemistry::ClassExtensions::ToMz(aok->NeutralMass, 1);
                auto c = 0.0000001;
                mz2.push_back(t1 + c * std::pow(t1, 2));
                auto t2 = Chemistry::ClassExtensions::ToMz(aok->NeutralMass + Constants::C13MinusC12, 1);
                mz2.push_back(t2 + c * std::pow(t2, 2));
                intensities2.push_back(1);
                intensities2.push_back(1);
            }
#ifdef ORIG
            auto MassSpectrum2 = new MzSpectrum(mz2.OrderBy([&] (std::any b)
            {
                return b;
            })->ToArray(), intensities2.ToArray(), false);
#endif
            std::sort(mz2.begin(), mz2.end() );
            auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false );
            
            MzLibUtil::MzRange tempVar2(0, 10000);
            auto scan2 = new MsDataScan(MassSpectrum2, 2, 2, true, Polarity::Positive, 2, &tempVar2, "df",
                                        MZAnalyzerType::Orbitrap, 234734, 1, std::vector<std::vector<double>>(), "scan=2",
                                        Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2), 2, 1,
                                        Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2), 2,
                                        DissociationType::HCD, 1,
                                        Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2));
            scan2->ComputeSelectedPeakIntensity(MassSpectrum1);
            scan2->ComputeMonoisotopicPeakIntensity(MassSpectrum1);
            MzLibUtil::MzRange tempVar3(0, 10000);

            auto ScansHere = new std::vector<MsDataScan*>();
            auto S = new MsDataScan(MassSpectrum1,1, 1, true, Polarity::Positive, 1, &tempVar3, "ff",
                                    MZAnalyzerType::Unknown, 1000,1, std::vector<std::vector<double>>(), "scan=1");
            ScansHere->push_back(S);
            ScansHere->push_back(scan2);

            Scans = *ScansHere; //.ToArray();
            
            //delete scan2;
            //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
            //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
            //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
            //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
        }
    }
    
    TestDataFile::TestDataFile(PeptideWithSetModifications *pepWithSetMods, int charge, double intensity, double rt) :
        MsDataFile(2,  new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        std::vector<double> mz1 = {Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), charge),
                                   Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + 1.003, charge),
                                   Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass() + 2.005, charge)};
        std::vector<double> intensities1 =  {intensity, intensity * 10, intensity / 10};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
        
        MzLibUtil::MzRange tempVar2(0, 10000);
        auto ScansHere = new std::vector<MsDataScan*>();            
        auto S = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, rt, &tempVar2, "ff",
                                MZAnalyzerType::Unknown, 1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(S);
        
        std::vector<double> mz2;
        std::vector<double> intensities2;
        for (auto aok : pepWithSetMods->Fragment(DissociationType::HCD,FragmentationTerminus::Both))
        {
            mz2.push_back(Chemistry::ClassExtensions::ToMz(aok->NeutralMass, 1));
            mz2.push_back(Chemistry::ClassExtensions::ToMz(aok->NeutralMass + 1.003, 1));
            intensities2.push_back(intensity);
            intensities2.push_back(intensity);
        }
#ifdef ORIG
        auto MassSpectrum2 = new MzSpectrum(mz2.OrderBy([&] (std::any b)
        {
            return b;
        })->ToArray(), intensities2.ToArray(), false);
#endif
        std::sort(mz2.begin(), mz2.end() );
        auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false );
        MzLibUtil::MzRange tempVar3(0, 10000);
        auto scan2 = new MsDataScan(MassSpectrum2, 2, 2, true, Polarity::Positive, rt + 0.01, &tempVar3, "df",
                                    MZAnalyzerType::Orbitrap, 234734, 1, std::vector<std::vector<double>>(), "scan=2",
                                    Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2), 2, 1,
                                    Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2), 2,
                                    DissociationType::HCD, 1,
                                    Chemistry::ClassExtensions::ToMz(pepWithSetMods->getMonoisotopicMass(), 2));
        scan2->ComputeSelectedPeakIntensity(MassSpectrum1);
        scan2->ComputeMonoisotopicPeakIntensity(MassSpectrum1);
        ScansHere->push_back(scan2);

        Scans = *ScansHere; //.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete scan2' statement was not added since scan2 was
        //passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    TestDataFile::TestDataFile(int MS3) :
        MsDataFile(2,  new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        std::vector<double> mz1 = {50, 60, 70, 80, 90,
                                   Chemistry::ClassExtensions::ToMz(764.1376, 2)};
        std::vector<double> intensities1 =  {1, 1, 1, 1, 1, 1};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

        MzLibUtil::MzRange tempVar2(0, 10000);
        auto ScansHere = new std::vector<MsDataScan*>();            
        auto S  = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, "ff",
                                 MZAnalyzerType::Unknown, 1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(S);
        
        std::vector<double> mz2 = {52, 62, 72, 147.0764, 257.1244, 258.127, 275.1350, 502};
        std::vector<double> intensities2 =  {1, 1, 1, 1, 1, 1, 1, 1};
        auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
        auto tempVar3 = new MsDataScan (MassSpectrum2, 2, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000),
                                        "f", MZAnalyzerType::Unknown, 100000, 1, std::vector<std::vector<double>>(), "scan=2",
                                        Chemistry::ClassExtensions::ToMz(764.1376, 2), 2, 1,
                                        Chemistry::ClassExtensions::ToMz(764.1376, 2), 2, DissociationType::CID, 1,
                                        Chemistry::ClassExtensions::ToMz(764.1376, 1));
        ScansHere->push_back(tempVar3);
        
        std::vector<double> mz3 = {53, 63, 73, 148.0764, 258.1244, 259.127, 276.1350, 503};
        std::vector<double> intensities3 =  {1, 1, 1, 1, 1, 1, 1, 1};
        auto MassSpectrum3 = new MzSpectrum(mz3, intensities3, false);
        auto tempVar4 = new MsDataScan (MassSpectrum3, 3, 2, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000),
                                        "f", MZAnalyzerType::Unknown, 100000, 1, std::vector<std::vector<double>>(), "scan=3",
                                        Chemistry::ClassExtensions::ToMz(764.1376, 2), 2, 1,
                                        Chemistry::ClassExtensions::ToMz(764.1376, 2), 2, DissociationType::ETD, 1,
                                        Chemistry::ClassExtensions::ToMz(764.1376, 1));
        ScansHere->push_back(tempVar4);
        
        std::vector<double> mz4 = {54, 64, 74, 149.0764, 259.1244, 260.127, 277.1350, 504};
        std::vector<double> intensities4 = {1, 1, 1, 1, 1, 1, 1, 1};
        auto MassSpectrum4 = new MzSpectrum(mz4, intensities4, false);
        auto tempVar5 = new MsDataScan (MassSpectrum4, 4, 3, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000),
                                        "f", MZAnalyzerType::Unknown, 100000, 1, std::vector<std::vector<double>>(), "scan=4",
                                        Chemistry::ClassExtensions::ToMz(275.1350, 1), 1, 1,
                                        Chemistry::ClassExtensions::ToMz(275.1350, 1), 1, DissociationType::HCD, 2,
                                        Chemistry::ClassExtensions::ToMz(275.1350, 1));
        ScansHere->push_back(tempVar5);
        
        std::vector<double> mz5 = {55, 65, 75, 150.0764, 260.1244, 261.127, 278.1350, 505};
        std::vector<double> intensities5 = {1, 1, 1, 1, 1, 1, 1, 1};
        auto MassSpectrum5 = new MzSpectrum(mz5, intensities5, false);
        auto tempVar6 = new MsDataScan(MassSpectrum5, 5, 3, true, Polarity::Positive, 2, new MzLibUtil::MzRange(0, 10000),
                                       "f", MZAnalyzerType::Unknown, 100000, 1, std::vector<std::vector<double>>(), "scan=5",
                                       Chemistry::ClassExtensions::ToMz(257.1244, 1), 1, 1,
                                       Chemistry::ClassExtensions::ToMz(257.1244, 1), 1, DissociationType::HCD, 2,
                                       Chemistry::ClassExtensions::ToMz(257.1244, 1));
        ScansHere->push_back(tempVar6);
        
        Scans = *ScansHere; //.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum5' statement was not added since
        //MassSpectrum5 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum4' statement was not added since
        //MassSpectrum4 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum3' statement was not added since
        //MassSpectrum3 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    TestDataFile::TestDataFile(std::vector<double> &ms2Mz, std::vector<double> &ms2Intensities, double precursorMass,
                               int precursorZ, double rt) :
        MsDataFile(2,  new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        std::vector<double> v1 = {Chemistry::ClassExtensions::ToMz(precursorMass, precursorZ),
                                  Chemistry::ClassExtensions::ToMz(precursorMass + 1.003, precursorZ)};
        std::vector<double> v2 = {1, 1};
        auto ms1 = new MzSpectrum(v1, v2, false);
        auto ms2 = new MzSpectrum(ms2Mz, ms2Intensities, false);
        
        MzLibUtil::MzRange tempVar2(0, 10000);
        auto ScansHere = new std::vector<MsDataScan*>();            
        auto S = new MsDataScan(ms1, 1, 1, true, Polarity::Positive, rt, &tempVar2, "ff", MZAnalyzerType::Unknown, 1000,
                                1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(S);

        MzLibUtil::MzRange tempVar3(0, 10000);
        auto S2 = new MsDataScan(ms2, 1, 2, true, Polarity::Positive, rt + 0.01, &tempVar3, "ff", MZAnalyzerType::Unknown,
                                 1000, 1, std::vector<std::vector<double>>(), "scan=2",
                                 Chemistry::ClassExtensions::ToMz(precursorMass, precursorZ), precursorZ, 1,
                                 Chemistry::ClassExtensions::ToMz(precursorMass, precursorZ), 1.0, DissociationType::HCD, 1,
                                 Chemistry::ClassExtensions::ToMz(precursorMass, precursorZ));
        ScansHere->push_back(S2);
        
        Scans = *ScansHere; //.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete ms2' statement was not added since ms2 was
        //passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete ms1' statement was not added since ms1 was
        //passed to a method or constructor. Handle memory management manually.
    }
    
    std::string TestDataFile::getFilePath() const
    {
        return "TestDataFile";
    }
    
    std::string TestDataFile::getName() const
    {
        return "TestDataFile";
    }
    
    void TestDataFile::ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities)
    {
        MzSpectrum *massSpectrum = new MzSpectrum(mz, intensities, false);
        Scans[0] = new MsDataScan(massSpectrum, Scans[0]->getOneBasedScanNumber(), Scans[0]->getMsnOrder(),
                                  Scans[0]->getIsCentroid(), Scans[0]->getPolarity(),
                                  Scans[0]->getRetentionTime(), Scans[0]->getScanWindowRange(),
                                  Scans[0]->getScanFilter(), Scans[0]->getMzAnalyzer(),
                                  massSpectrum->getSumOfAllY(),
                                  Scans[0]->getInjectionTime(), std::vector<std::vector<double>>(),
                                  Scans[0]->getNativeId());
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete massSpectrum' statement was not added since
        //massSpectrum was passed to a method or constructor. Handle memory management manually.
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
