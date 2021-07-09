#include "CsmSerializationTest.h"
#include "../EngineLayer/CommonParameters.h"
#include "../TaskLayer/XLSearchTask/XLSearchParameters.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/XLSearchTask/XLSearchTask.h"
#include "../TaskLayer/DbForTask.h"
#include "../TaskLayer/EverythingRunnerEngine.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/CrosslinkSearch/CrosslinkedPeptides.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::CrosslinkSearch;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::AminoAcidPolymer;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;


#include "MzLibAssert.h"
#include <filesystem>
#include <iostream>
#include <fstream>

static bool CompareFiles (std::string filename1, std::string filename2)
{
    std::ifstream file1, file2;

    file1.open(filename1);
    if ( !file1.is_open() ) {
        std::cout << "Could not open file " << filename1 << std::endl;
        return false;
    }

    file2.open(filename2);
    if ( !file2.is_open() ) {
        std::cout << "Could not open file " << filename2 << std::endl;
        return false;
    }

    int lineno=0;
    while ( true ) {
        std::string line1, line2;
        if ( std::getline(file1, line1) ) {
            if ( std::getline (file2, line2) ) {
                if ( line1 == line2 ) {
                    continue;
                }
                else {
                    std::cout << "Lines don't match. Line number : " << lineno << std::endl;
                    std::cout << "   File 1: " << filename1 << std::endl << line1 << std::endl; 
                    std::cout << "   File 2: " << filename2 << std::endl << line2 << std::endl;
                    return false;
                }
            }
            else {
                // reach end of file1, but not file1
                return false;
            }
        }
        else {
            if ( !std::getline(file2, line2 )) {
                // reached end of both files
                return true;
            }
            else {
                // reached end of file1, but not file2
                return false;
            }
        }
        lineno++;
    }
    return true;
}


int main ( int argc, char **argv )
{
    int i=0;
    std::cout << i << ". PeriodicTableLoader" << std::endl;
    const std::string elfile="elements.dat";
    const std::string &elr=elfile;
    //Chemistry::PeriodicTable::Load (elr);
    UsefulProteomicsDatabases::PeriodicTableLoader::Load (elr);

    std::cout << ++i << ". CSMSerialization_BSA_DSSO" << std::endl;
    Test::CSMSerializationTest::CSMSerializationTest_BSA_DSSO();

    return 0;
}


namespace Test
{

    void CSMSerializationTest::CSMSerializationTest_BSA_DSSO()
    {
        std::string testdir=std::filesystem::current_path().string();

        //Generate parameters
        PpmTolerance tempVar(10);
        auto commonParameters = new CommonParameters("", DissociationType::EThcD, false, true, 3, 12, true, false, 1, 1,
                                                     200, 0.01, false, true, false, false, nullptr,
                                                     &tempVar, nullptr, -1, new DigestionParams("trypsin", 2, 5) );
            
        auto xlSearchParameters = new XlSearchParameters();
        
        //Create databases contain two protein.
        std::vector<Protein*> proteinList = 
            {
                new Protein("EKVLTSSAR", "Fake01"),
                new Protein("LSQKFPK", "Fake02")
            };
        
        ModificationMotif *motif1;
        ModificationMotif::TryGetMotif("M", &motif1);
        Modification *mod1 = new Modification("Oxidation of M", "", "Common Variable", "", motif1,
                                              "Anywhere.", nullptr, std::make_optional((double)15.99491461957));
        
        ModificationMotif *motif2;
        ModificationMotif::TryGetMotif("C", &motif2);
        Modification *mod2 = new Modification( "Carbamidomethyl of C", "", "Common Fixed", "",  motif2,
                                               "Anywhere.", nullptr, std::make_optional((double)57.02146372068994));
        
        auto variableModifications = std::vector<Modification*> {mod1};
        auto fixedModifications = std::vector<Modification*> {mod2};
        std::vector<Modification*> localizeableModifications;
        
        //Run index engine
        std::vector<std::string>pDbs, nestedIds;
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse,
                                              commonParameters, 30000, false, pDbs, nestedIds);
        
        IndexingResults* indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        
        std::vector<std::vector<int>> tvvec = indexResults->getFragmentIndex();
        std::vector<int> indexedFragments;
        for ( auto p : tvvec ) {
            if ( !p.empty() ) {
                // Do the SelectMany flatmap operation
                for ( auto q: p ) {
                    indexedFragments.push_back(q);
                }
            }
        }
        
        Assert::AreEqual(82, (int)indexedFragments.size());
        Assert::AreEqual(3, (int)indexResults->getPeptideIndex().size());
        
        //Get MS2 scans.
        auto myMsDataFile = new XLTestDataFile();
        CommonParameters tempVar3;
        std::vector<Ms2ScanWithSpecificMass*>listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar3);
        std::sort (listOfSortedms2Scans.begin(), listOfSortedms2Scans.end(), [&]
                   (Ms2ScanWithSpecificMass*l, Ms2ScanWithSpecificMass*r) {
                       return l->getPrecursorMass() < r->getPrecursorMass(); });
        
        //Generate crosslinker, which is DSSO here.
        Crosslinker tempVar4;
        Crosslinker *crosslinker = (&tempVar4)->SelectCrosslinker(CrosslinkerType::DSSO);
        
        std::vector<CrosslinkSpectralMatch*> possiblePsms(listOfSortedms2Scans.size());
        std::vector<std::string> nIds;
        auto pvar = indexResults->getPeptideIndex();
        auto pvar2 = indexResults->getFragmentIndex();

        CrosslinkSearchEngine tempVar5(possiblePsms, listOfSortedms2Scans, pvar,
                                       pvar2, 0, commonParameters,
                                       crosslinker, xlSearchParameters->getRestrictToTopNHits(),
                                       xlSearchParameters->getCrosslinkSearchTopNum(),
                                       xlSearchParameters->getXlQuench_H2O(),
                                       xlSearchParameters->getXlQuench_NH2(),
                                       xlSearchParameters->getXlQuench_Tris(),
                                       nIds);
        (&tempVar5)->Run();
        std::vector<CrosslinkSpectralMatch*> newPsms;
        for ( auto p:  possiblePsms ) {
            if ( p != nullptr ) {
                newPsms.push_back(p);
            }
        }

        for (auto item : newPsms)		{
            item->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
        }
        
        //Test newPsms
        Assert::AreEqual(3, (int)newPsms.size());
        
        //Test Output
        auto task = new XLSearchTask();
        std::vector<std::string> vs, vs2;
        task->WritePepXML_xl(newPsms, proteinList, "", variableModifications, fixedModifications, vs,
                             testdir, "pep.XML", vs2);

        //std::filesystem::rename("singlePsms.tsv", "singlePsms.orig.tsv");
        //std::filesystem::rename("allPsms.tsv", "allPsms.orig.tsv");
        std::filesystem::rename("pep.XML.pep.XM", "pep.orig.xml");
        
        size_t bufsize= 1024*1024;
        char *sbuf  = new char[bufsize];

        int ret = CrosslinkSpectralMatch::Pack(sbuf, bufsize, newPsms);
        std::cout << "After Pack, ret = " << ret << " bufsize = " << bufsize << std::endl;
        std::ofstream output("CSM.out");
        output << sbuf;
        output.close();

        output.open("CsmOrig.out");
        for ( auto psms : newPsms ) {
            output << psms->ToString() << std::endl;
        }
        output.close();

        //std::vector<Protein *> pList = task->getProteinList();
        
        if ( ret > 0 ) {
            std::vector<CrosslinkSpectralMatch*> unpackedPsms;
            int count=-1;
            size_t len=0;
            CrosslinkSpectralMatch::Unpack( sbuf, bufsize, count, len, unpackedPsms, listOfSortedms2Scans, proteinList);
            std::cout << "len = " << len << " veclen = " << unpackedPsms.size() << std::endl;

            output.open("CsmSerialized.out");
            for ( auto psms : unpackedPsms ) {
                output << psms->ToString() << std::endl;
            }
            output.close();

            
            //task->WritePepXML_xl(unpackedPsms, proteinList, "", variableModifications, fixedModifications, vs,
            //                     testdir, "pep.XML", vs2);

            //Assert::IsTrue(CompareFiles("singlePsms.tsv", "singlePsms.orig.tsv"));
            //Assert::IsTrue(CompareFiles("allPsms.tsv", "allPsms.orig.tsv"));
            //Assert::IsTrue(CompareFiles("pep.XML.pep.XM", "pep.orig.xml"));
        }

        delete[] sbuf;
        //std::filesystem::remove("singlePsms.tsv");
        //std::filesystem::remove("singlePsms.orig.tsv");
        //std::filesystem::remove("allPsms.tsv");
        //std::filesystem::remove("allPsms.orig.tsv");
        std::filesystem::remove("pep.XML.pep.XM");
        std::filesystem::remove("pep.orig.xml");

        delete task;
        delete myMsDataFile;
        delete indexEngine;
        delete mod2;
        delete mod1;
        delete xlSearchParameters;
        delete commonParameters;
    }

    XLTestDataFile::XLTestDataFile() : MsDataFile(2, new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        auto mz1 = std::vector<double> {Chemistry::ClassExtensions::ToMz(1994.05,3), Chemistry::ClassExtensions::ToMz(846.4963,1),
                                        Chemistry::ClassExtensions::ToMz(1004.495,1), Chemistry::ClassExtensions::ToMz(1093.544,1),
                                        Chemistry::ClassExtensions::ToMz(1043.561,1)};
        auto intensities1 = std::vector<double> {1, 1, 1, 1, 1};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
        auto ScansHere = new std::vector<MsDataScan*>();
        MzLibUtil::MzRange tempVar2(0, 10000);
        auto s = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, "ff", MZAnalyzerType::Unknown,
                                1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(s);
        
        auto mz2 = std::vector<double> {100, 201.1234, 244.1656, 391.2340, 420.2201, 521.2678, 634.3519, 889.965, 1044.568,
                                        1094.551, 1279.671, 1378.74, 1491.824};
        auto intensities2 = std::vector<double> {100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
        auto tempVar3 = new MsDataScan (MassSpectrum2, 2, 2, true, Polarity::Positive, 1.0, new MzLibUtil::MzRange(0, 10000), "f",
                                        MZAnalyzerType::Unknown, 112, 1.0, std::vector<std::vector<double>>(), "scan=2",
                                        Chemistry::ClassExtensions::ToMz(1994.05,3), 3, 1,
                                        Chemistry::ClassExtensions::ToMz(1994.05,3),
                                        2, DissociationType::HCD, 1, Chemistry::ClassExtensions::ToMz(1994.05,3));
        ScansHere->push_back(tempVar3);
        
        auto mz3 = std::vector<double> {100, 201.1234, 244.1656, 391.2340};
        auto intensities3 = std::vector<double> {100, 1, 1, 1};
        auto MassSpectrum3 = new MzSpectrum(mz3, intensities3, false);
        auto tempVar4 = new MsDataScan (MassSpectrum3, 3, 2, true, Polarity::Positive, 1.0, new MzLibUtil::MzRange(0, 10000),
                                        "f", MZAnalyzerType::Unknown, 103, 1.0, std::vector<std::vector<double>>(), "scan=3",
                                        Chemistry::ClassExtensions::ToMz(846.4963,1), 1, 1,
                                        Chemistry::ClassExtensions::ToMz(846.4963,1), 2, DissociationType::HCD, 1,
                                        Chemistry::ClassExtensions::ToMz(846.4963,1));
        ScansHere->push_back(tempVar4);
        
        auto mz4 = std::vector<double> {100, 201.1234, 244.1656, 391.2340};
        auto intensities4 = std::vector<double> {100, 1, 1, 1};
        auto MassSpectrum4 = new MzSpectrum(mz4, intensities4, false);
        auto tempVar5 = new MsDataScan (MassSpectrum4, 4, 2, true, Polarity::Positive, 1.0, new MzLibUtil::MzRange(0, 10000),
                                        "f", MZAnalyzerType::Unknown, 103, 1.0, std::vector<std::vector<double>>(), "scan=4",
                                        Chemistry::ClassExtensions::ToMz(1004.491,1), 1, 1,
                                        Chemistry::ClassExtensions::ToMz(1004.491,1), 2, DissociationType::HCD, 1,
                                        Chemistry::ClassExtensions::ToMz(1004.491,1));
        ScansHere->push_back(tempVar5);
        
        Scans = *ScansHere;//.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum4' statement was not added since
        //MassSpectrum4 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum3' statement was not added since
        //MassSpectrum3 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    std::string XLTestDataFile::getFilePath() const
    {
        return "XLTestDataFile";
    }
    
    std::string XLTestDataFile::getName() const
    {
        return "XLTestDataFile";
    }
    
    void XLTestDataFile::ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities)
    {
        MzSpectrum *massSpectrum = new MzSpectrum(mz, intensities, false);
        Scans[0] = new MsDataScan(massSpectrum, Scans[0]->getOneBasedScanNumber(), Scans[0]->getMsnOrder(), Scans[0]->getIsCentroid(),
                                  Scans[0]->getPolarity(), Scans[0]->getRetentionTime(), Scans[0]->getScanWindowRange(),
                                  Scans[0]->getScanFilter(), Scans[0]->getMzAnalyzer(), massSpectrum->getSumOfAllY(),
                                  Scans[0]->getInjectionTime(), std::vector<std::vector<double>>(),
                                  Scans[0]->getNativeId() );
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete massSpectrum' statement was not added since
        //massSpectrum was passed to a method or constructor. Handle memory management manually.
    }

}
