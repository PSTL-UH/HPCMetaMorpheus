#include "XLTest.h"
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


#include "Assert.h"
#include <filesystem>
#include <iostream>
#include <fstream>

int main ( int argc, char **argv )
{
    int i=0;
    std::cout << i << ". PeriodicTableLoader" << std::endl;
    const std::string elfile="elements.dat";
    const std::string &elr=elfile;
    //Chemistry::PeriodicTable::Load (elr);
    UsefulProteomicsDatabases::PeriodicTableLoader::Load (elr);

    std::cout << ++i << ". XlTestXlPosCal" << std::endl;
    Test::XLTest::XlTestXlPosCal();

    std::cout << ++i << ". XlTestGenerateIntensityRanks" << std::endl;
    Test::XLTest::XlTestGenerateIntensityRanks();
    
    std::cout << ++i << ". XlTest_BSA_DSSO" << std::endl;
    Test::XLTest::XlTest_BSA_DSSO();

    std::cout << ++i << ". XlTest_BSA_DSS_file" << std::endl;
    Test::XLTest::XlTest_BSA_DSS_file();

    std::cout << ++i << ". XlTest_GenerateUserDefinedCrosslinker" << std::endl;
    Test::XLTest::XlTest_GenerateUserDefinedCrosslinker();

    std::cout << ++i << ". XlTest_DiffCrosslinkSites" << std::endl;
    Test::XLTest::XlTest_DiffCrosslinkSites();

    std::cout << ++i << ". CrosslinkCreateTest" << std::endl;
    Test::XLTest::CrosslinkCreateTest();
    
    std::cout << ++i << ". DeadendPeptideTest" << std::endl;
    Test::XLTest::DeadendPeptideTest();

#ifdef NOT_NOW
    // This test compiles correctly, and runs to large extent, however, we segfault
    // because the peptideindex is not yet written correctly, this is an
    // issue with the cereal interfaces for this class. To be revisited.
    std::cout << ++i << ". XLSearchWithGeneratedIndices" << std::endl;
    Test::XLTest::XLSearchWithGeneratedIndices();
#endif
    std::cout << ++i << ". TestGetPossibleCrosslinkerSites" << std::endl;
    Test::XLTest::TestGetPossibleCrosslinkerSites();
    
    std::cout << ++i << ". TestTheoreticalFragmentsLoop" << std::endl;
    Test::XLTest::TestTheoreticalFragmentsLoop();

    std::cout << ++i << ". TestTheoreticalLoopFragmentsWithMod" << std::endl;
    Test::XLTest::TestTheoreticalLoopFragmentsWithMod();

#ifdef LATER
    std::cout << ++i << ". TestDeadendTris" << std::endl;
    Test::XLTest::TestDeadendTris();

    std::cout << ++i << ". TestTheoreticalFragmentsNonCleavableCrosslink" << std::endl;
    Test::XLTest::TestTheoreticalFragmentsNonCleavableCrosslink();

    std::cout << ++i << ". TestTheoreticalFragmentsCleavableCrosslink" << std::endl;
    Test::XLTest::TestTheoreticalFragmentsCleavableCrosslink();

    std::cout << ++i << ". TestWriteToPercolator" << std::endl;
    Test::XLTest::TestWriteToPercolator();

    std::cout << ++i << ". TestWriteNonSingleCross" << std::endl;
    Test::XLTest::TestWriteNonSingleCross();
#endif        

    return 0;
}

namespace Test
{

    IndexingResults *XLTest::privateindexResults;
    CommonParameters *XLTest::privatecommonParameters;
    TaskLayer::XlSearchParameters *XLTest::privatexlSearchParameters;
    std::vector<Protein*> XLTest::privateproteinList;
    std::vector<Modification*> XLTest::privatevariableModifications;
    std::vector<Modification*> XLTest::privatefixedModifications;
    Crosslinker *XLTest::privatecrosslinker;
    std::vector<PeptideWithSetModifications*> XLTest::privatedigestedList;
    
    IndexingResults *XLTest::getindexResults()
    {
        return privateindexResults;
    }
    
    void XLTest::setindexResults(IndexingResults *value)
    {
        privateindexResults = value;
    }
    
    CommonParameters *XLTest::getcommonParameters()
    {
        return privatecommonParameters;
    }
    
    void XLTest::setcommonParameters(CommonParameters *value)
    {
        privatecommonParameters = value;
    }
    
    XlSearchParameters *XLTest::getxlSearchParameters()
    {
        return privatexlSearchParameters;
    }
    
    void XLTest::setxlSearchParameters(XlSearchParameters *value)
    {
        privatexlSearchParameters = value;
    }
    
    std::vector<Protein*> XLTest::getproteinList()
    {
        return privateproteinList;
    }
    
    void XLTest::setproteinList(const std::vector<Protein*> &value)
    {
        privateproteinList = value;
    }
    
    std::vector<Modification*> XLTest::getvariableModifications()
    {
        return privatevariableModifications;
    }
    
    void XLTest::setvariableModifications(const std::vector<Modification*> &value)
    {
        privatevariableModifications = value;
    }
    
    std::vector<Modification*> XLTest::getfixedModifications()
    {
        return privatefixedModifications;
    }
    
    void XLTest::setfixedModifications(const std::vector<Modification*> &value)
    {
        privatefixedModifications = value;
    }
    
    Crosslinker *XLTest::getcrosslinker()
    {
        return privatecrosslinker;
    }
    
    void XLTest::setcrosslinker(Crosslinker *value)
    {
        privatecrosslinker = value;
    }
    
    std::vector<PeptideWithSetModifications*> XLTest::getdigestedList()
    {
        return privatedigestedList;
    }
    
    void XLTest::setdigestedList(const std::vector<PeptideWithSetModifications*> &value)
    {
        privatedigestedList = value;
    }
    
    void XLTest::XlTestXlPosCal()
    {
        auto prot = new Protein("MNNNKQQQQ", "");
        std::vector<DigestionMotif*> motifs = {new DigestionMotif("K", "", 1, "")};
        Protease *protease = new Protease("New Custom Protease", CleavageSpecificity::Full, "", "", motifs);
        //ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        ProteaseDictionary::insert(protease->getName(), protease);
            
        DigestionParams *digestionParams = new DigestionParams(protease->getName(), 2, 1,  std::numeric_limits<int>::max(),
                                                               1024, InitiatorMethionineBehavior::Retain);
        std::vector<Modification*> variableModifications;

        std::vector<Modification*> tvec1;
        auto ye = prot->Digest(digestionParams, tvec1, variableModifications);
        
        auto pep = ye[0];

        Assert::AreEqual(pep->getBaseSequence(), "MNNNK");
        Crosslinker *crosslinker = new Crosslinker();
        crosslinker->SelectCrosslinker(CrosslinkerType::DSS);
        Assert::AreEqual(crosslinker->getCrosslinkerModSites(), "K");        
        Assert::AreEqual(Residue::GetResidue(crosslinker->getCrosslinkerModSites())->getMonoisotopicMass(), 128.09496301518999, 1e-9);
        auto n = pep->Fragment(DissociationType::HCD, FragmentationTerminus::N);
        auto c = pep->Fragment(DissociationType::HCD, FragmentationTerminus::C);
        Assert::AreEqual((int)n.size(), 4);
        Assert::AreEqual((int)c.size(), 4);
        Assert::AreEqual(c.front()->NeutralMass, 146.10552769899999, 1e-6);

        std::string s = crosslinker->getCrosslinkerModSites();
        std::vector<char> cvec(s.begin(), s.end() );
        auto x = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(cvec, pep);        
        Assert::AreEqual(x[0], 5);

        auto pep2 = ye[2];
        Assert::AreEqual("MNNNKQQQQ", pep2->getBaseSequence());
        auto n2 = pep2->Fragment(DissociationType::HCD, FragmentationTerminus::N);
        auto c2 = pep2->Fragment(DissociationType::HCD, FragmentationTerminus::C);
        Assert::AreEqual((int)n2.size(), 8);
        Assert::AreEqual((int)c2.size(), 8);

        auto x2 = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites(cvec, pep2);//.ToArray();
        Assert::AreEqual(x2[0], 5);
        
        //Test crosslinker with multiple types of mod
        auto protSTC = new Protein("GASTACK", "");
        std::vector<Modification*> tvec2;
        auto peps = protSTC->Digest(digestionParams, tvec2, variableModifications);
        auto pepSTC = peps[0];
        Assert::AreEqual(pepSTC->getBaseSequence(), "GASTACK");

        Crosslinker *crosslinker2 = new Crosslinker("ST", "C", "crosslinkerSTC", false, -18.01056, 0, 0, 0, 0, 0, 0);

        std::string s1 = crosslinker2->getCrosslinkerModSites();
        std::string s2 = crosslinker2->getCrosslinkerModSites2();
        std::vector<char> cvec2( s1.begin(), s1.end() );
        cvec2.insert(cvec2.end(), s2.begin(), s2.end() );
        std::vector<char> cvec3;
        for ( auto p : cvec2 ) {
            bool found = false;
            for ( auto q : cvec3 ) {
                if ( p == q ) {
                    found = true;
                    break;
                }
            }
            if (!found ) {
                cvec3.push_back(p);
            }
        }
        std::string crosslinkerModSitesAll (cvec3.begin(), cvec3.end());
        Assert::AreEqual(crosslinkerModSitesAll, "STC");

        delete crosslinker2;
        delete protSTC;
        delete crosslinker;
        delete digestionParams;
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease
        //was passed to a method or constructor. Handle memory management manually.
        delete prot;
    }


    void XLTest::XlTestGenerateIntensityRanks()
    {
        std::vector<double> intensity = {1.1, 1.1, 0.5, 3.2, 0.5, 6.0};
        std::vector<int> rank = CrosslinkSpectralMatch::GenerateIntensityRanks(intensity);
        std::vector<int> Rank = {4, 3, 6, 2, 5, 1};
        Assert::AreEqual(rank.size(), Rank.size());
        if ( rank.size() == Rank.size() ) {
            for ( int i=0; i< (int)rank.size(); i++ ) {
                Assert::AreEqual(rank[i], Rank[i] );
            }
        }
    }
    
    void XLTest::XlTest_BSA_DSSO()
    {
        std::string testdir=std::filesystem::current_path().string();

        //Generate parameters
        PpmTolerance tempVar(10);
#ifdef ORIG
        DigestionParams tempVar2(minPeptideLength: 5);
        auto commonParameters = new CommonParameters(, DissociationType::EThcD, false, , , , , , , 1,
                                                     , , , , , , ,
                                                     &tempVar, , , &tempVar2);
#endif
        auto commonParameters = new CommonParameters("", DissociationType::EThcD, false, true, 3, 12, true, false, 1, 1,
                                                     200, 0.01, false, true, false, false, nullptr,
                                                     &tempVar, nullptr, -1, new DigestionParams("trypsin", 2, 5) );
            
        auto xlSearchParameters = new XlSearchParameters();
        
        //Create databases contain two protein.
        auto proteinList = std::vector<Protein*>
            {
                new Protein("EKVLTSSAR", "Fake01"),
                new Protein("LSQKFPK", "Fake02")
            };
        
        ModificationMotif *motif1;
        ModificationMotif::TryGetMotif("M", &motif1);
#ifdef ORIG
        Modification *mod1 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1,
                                              _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
#endif
        Modification *mod1 = new Modification("Oxidation of M", "", "Common Variable", "", motif1,
                                              "Anywhere.", nullptr, std::make_optional((double)15.99491461957));
        
        ModificationMotif *motif2;
        ModificationMotif::TryGetMotif("C", &motif2);
#ifdef ORIG
        Modification *mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif2,
                                              _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
#endif
        Modification *mod2 = new Modification( "Carbamidomethyl of C", "", "Common Fixed", "",  motif2,
                                               "Anywhere.", nullptr, std::make_optional((double)57.02146372068994));
        
        auto variableModifications = std::vector<Modification*> {mod1};
        auto fixedModifications = std::vector<Modification*> {mod2};
        auto localizeableModifications = std::vector<Modification*>();
        
        //Run index engine
        std::vector<std::string>pDbs, nestedIds;
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse,
                                              commonParameters, 30000, false, pDbs, nestedIds);
        
        IndexingResults* indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        
#ifdef ORIG
        auto indexedFragments = indexResults->getFragmentIndex().Where([&] (std::any p)    {
                return p != nullptr;
            }).SelectMany([&] (std::any v)   {
                    return v;
                }).ToList();
        
#endif
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
#ifdef ORIG
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar3).OrderBy([&] (std::any b) {
                b::PrecursorMass;
            })->ToArray();
#endif
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
#ifdef ORIG
        auto newPsms = possiblePsms.Where([&] (std::any p)     {
                return p != nullptr;
            }).ToList();
#endif
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
        
        // EDGAR: THese line were already commented out in the C# version
        //Test PsmCross.XlCalculateTotalProductMasses
        //var psmCrossAlpha = new CrosslinkSpectralMatch(digestedList[1], 0, 0, 0, listOfSortedms2Scans[0],
        //                                               commonParameters.DigestionParams,
        //                                               new List<MatchedFragmentIon>());
        //var psmCrossBeta = new CrosslinkSpectralMatch(digestedList[2], 0, 0, 0, listOfSortedms2Scans[0],
        //                                              commonParameters.DigestionParams,
        //                                              new List<MatchedFragmentIon>());
        //var linkPos = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(),
        //                                                                    digestedList[1]);
        //var productMassesAlphaList = CrosslinkedPeptide.XlGetTheoreticalFragments(DissociationType.EThcD, false, crosslinker,
        //                                                                          linkPos, digestedList[2].MonoisotopicMass,
        //                                                                          digestedList[1]);
        //Assert.AreEqual(productMassesAlphaList.First().Value.Count, 50); //TO DO: The number here should be manually verified.
        // EDGAR: END of commented out section

        std::filesystem::remove("singlePsms.tsv");
        std::filesystem::remove("pep.XML.pep.xml");
        std::filesystem::remove("allPsms.tsv");
        
        delete task;
        delete myMsDataFile;
        delete indexEngine;
        delete mod2;
        delete mod1;
        delete xlSearchParameters;
        delete commonParameters;
    }
    

    void XLTest::XlTest_BSA_DSS_file()
    {
        std::string testdir=std::filesystem::current_path().string();
        
        auto task = new XLSearchTask(testdir + "/XlTestData/XLSearchTaskconfig_BSA_DSS_23747.toml");
        
        FileSystem::createDirectory(testdir + "/TESTXlTestData");
        DbForTask *db = new DbForTask(testdir + "/XlTestData/BSA.fasta", false);
        std::string raw = testdir + "/XlTestData/BSA_DSS_23747.mzML";
        std::vector<std::tuple<std::string, TaskLayer::MetaMorpheusTask*>> vec1 = {std::make_tuple("Task", task)};
        std::vector<std::string> svec = {raw};
        std::vector<DbForTask*> dbvec =  {db};
        EverythingRunnerEngine tempVar(vec1, svec, dbvec, testdir + "/TESTXlTestData");
        (&tempVar)->Run();
        std::filesystem::remove(testdir +  "/TESTXlTestData");

        delete task;
        delete db;
    }

    
    void XLTest::XlTest_GenerateUserDefinedCrosslinker()
    {
        XlSearchParameters *xlSearchParameters = new XlSearchParameters();
        xlSearchParameters->setCrosslinkerType(CrosslinkerType::UserDefined);
        xlSearchParameters->setCrosslinkerName("CrossST-C");
        xlSearchParameters->setCrosslinkerResidues("ST");
        xlSearchParameters->setCrosslinkerResidues2("C");
        xlSearchParameters->setCrosslinkerTotalMass(std::make_optional(-18.01056));
        auto crosslinker = XLSearchTask::GenerateUserDefinedCrosslinker(xlSearchParameters);
        
        delete xlSearchParameters;
    }
    
    void XLTest::XlTest_DiffCrosslinkSites()
    {
        std::string testdir=std::experimental::filesystem::current_path().string();
        //Generate parameters
#ifdef ORIG
        DigestionParams *tempVar = new DigestionParams (minPeptideLength: 4);
        auto commonParameters = new CommonParameters(, , false, , , , , , , 1, , , , , , , , , , , &tempVar);
#endif
        DigestionParams *tempVar = new DigestionParams ("trypsin", 2, 4);        
        auto commonParameters = new CommonParameters("", DissociationType::HCD, false, true, 3, 12, true, false, 1, 1,
                                                     200, 0.01, false, true, false, false, nullptr,
                                                     nullptr,  nullptr, -1, tempVar);
        
        auto xlSearchParameters = new XlSearchParameters();
        xlSearchParameters->setCrosslinkerType(CrosslinkerType::UserDefined);
        xlSearchParameters->setCrosslinkerName("CrossST-C");
        xlSearchParameters->setCrosslinkerResidues("ST");
        xlSearchParameters->setCrosslinkerResidues2("C");
        xlSearchParameters->setCrosslinkerTotalMass(std::make_optional(-18.01056));
        
        //Create databases contain two protein.
        auto proteinList = std::vector<Protein*>
            {
                new Protein("VLTAR", "Fake01"),
                new Protein("LCQK", "Fake02")
            };
        
        ModificationMotif *motif1;
        ModificationMotif::TryGetMotif("M", &motif1);
#ifdef ORIG
        Modification *mod1 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1,
                                              _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
#endif
        Modification *mod1 = new Modification("Oxidation of M", "", "Common Variable", "", motif1,
                                              "Anywhere.", nullptr, std::make_optional((double)15.99491461957));
        
        std::vector<Modification*> variableModifications = {mod1};
        std::vector<Modification*> fixedModifications;
        std::vector<Modification*> localizeableModifications;
        
        std::unordered_map<Modification*, unsigned short> modsDictionary;
        
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
        
        //Generate digested peptide lists.
        std::vector<PeptideWithSetModifications*> digestedList;
        for (auto item : proteinList)
        {
            auto digested = item->Digest(commonParameters->getDigestionParams(), fixedModifications, variableModifications);
            digestedList.insert(digestedList.end(), digested.begin(), digested.end());
        }
        
        //Run index engine
        std::vector<std::string>pDbs, nestedIds;
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse,
                                              commonParameters, 30000, false, pDbs, nestedIds );
        
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        
        //Get MS2 scans.
        auto myMsDataFile = new XLTestDataFileDiffSite();
        CommonParameters tempVar2;
#ifdef ORIG
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
#endif
        std::vector<Ms2ScanWithSpecificMass*>listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2);
        std::sort (listOfSortedms2Scans.begin(), listOfSortedms2Scans.end(), [&]
                   (Ms2ScanWithSpecificMass*l, Ms2ScanWithSpecificMass*r) {
                       return l->getPrecursorMass() < r->getPrecursorMass(); });
        
        //Generate crosslinker, which is UserDefined here.
        auto crosslinker = XLSearchTask::GenerateUserDefinedCrosslinker(xlSearchParameters);
        
        //TwoPassCrosslinkSearchEngine.Run().
        std::vector<CrosslinkSpectralMatch*> possiblePsms(listOfSortedms2Scans.size());
        std::vector<std::string> nIds;
        auto pvar = indexResults->getPeptideIndex();
        auto pvar2 = indexResults->getFragmentIndex();
        CrosslinkSearchEngine tempVar3(possiblePsms, listOfSortedms2Scans, pvar,
                                       pvar2, 0, commonParameters, crosslinker,
                                       xlSearchParameters->getRestrictToTopNHits(), xlSearchParameters->getCrosslinkSearchTopNum(),
                                       xlSearchParameters->getXlQuench_H2O(), xlSearchParameters->getXlQuench_NH2(),
                                       xlSearchParameters->getXlQuench_Tris(), nIds);
        (&tempVar3)->Run();
        
#ifdef ORIG
        auto newPsms = possiblePsms.Where([&] (std::any p)    {
                return p != nullptr;
            }).ToList();
#endif
        std::vector<CrosslinkSpectralMatch*> newPsms;
        for ( auto p:  possiblePsms ) {
            if ( p != nullptr ) {
                newPsms.push_back(p);
            }
        }
        
        Assert::AreEqual(1, (int)newPsms.size());
        
        delete myMsDataFile;
        delete indexEngine;
        delete mod1;
        delete xlSearchParameters;
        delete commonParameters;
    }
    
    void XLTest::CrosslinkCreateTest()
    {
        XlSearchParameters tempVar;
        auto var = XLSearchTask::GenerateUserDefinedCrosslinker(&tempVar);
        // slight change from the C# version of the test
        Assert::IsTrue( dynamic_cast<Crosslinker*>(var) != nullptr );
    }
    
    void XLTest::DeadendPeptideTest()
    {
        std::string testdir=std::filesystem::current_path().string();

        std::string myFileXl = testdir + "/XlTestData/BSA_DSSO_ETchD6010.mgf";
        std::string myDatabaseXl = testdir+ "/XlTestData/BSA.fasta";
        std::string outputFolder = testdir+ "/TestDeadendPeptide";
        
        XLSearchTask *xLSearchTask = new XLSearchTask();
        auto tempVar = new CommonParameters( "", DissociationType::HCD, true, true, 3, 12, true, false, 1, 5, 200, 0.01, false,
                                  true, false, false, nullptr, new PpmTolerance(51000));        
        xLSearchTask->setCommonParameters(tempVar);
        
        XLSearchTask *xLSearchTask2 = new XLSearchTask();
        auto tempVar2 = new CommonParameters ( "", DissociationType::HCD, true, true, 3, 12, true, false, 1, 5, 200, 0.01, false,
                                               true, false, false, nullptr, new PpmTolerance(112000));
        xLSearchTask2->setCommonParameters(tempVar2);

        XlSearchParameters* tempVar3 = new XlSearchParameters();
        xLSearchTask2->setXlSearchParameters(tempVar3);
        xLSearchTask2->getXlSearchParameters()->setXlQuench_Tris(false);
        xLSearchTask2->getXlSearchParameters()->setXlQuench_H2O(false);
        xLSearchTask2->getXlSearchParameters()->setXlQuench_NH2(true);
        FileSystem::createDirectory(outputFolder);

        DbForTask *db = new DbForTask(myDatabaseXl, false);
        std::vector<DbForTask *> dbvec = {db};
        std::vector<std::string> filevec = {myFileXl};
        xLSearchTask->RunTask(outputFolder, dbvec, filevec, "test");
        xLSearchTask2->RunTask(outputFolder, dbvec, filevec, "test");

        //std::filesystem::remove_all(testdir +  "/Task Settings");
        //std::filesystem::remove_all(outputFolder);

        delete xLSearchTask2;
        delete xLSearchTask;
    }
    
    void XLTest::XLSearchWithGeneratedIndices()
    {
        std::string testdir=std::filesystem::current_path().string();
        
        XLSearchTask *xlSearchTask = new XLSearchTask();
        CommonParameters tempVar;
        xlSearchTask->setCommonParameters(&tempVar);
        std::string myFile = testdir + "/XlTestData/BSA_DSSO_ETchD6010.mgf";
        std::string myDatabase = testdir + "/XlTestData/BSA.fasta";
        std::string folderPath = testdir+ "/TestXLSearch";
        
        DbForTask *db = new DbForTask(myDatabase, false);
        std::vector<std::tuple<std::string, MetaMorpheusTask *>> taskList = {std::make_tuple("TestXLSearch", xlSearchTask)};

        FileSystem::createDirectory(folderPath);
        
        //creates .params files if they do not exist
        std::vector<DbForTask *> dbVec = {db};
        std::vector<std::string> FileVec = {myFile};
        xlSearchTask->RunTask(folderPath, dbVec, FileVec, "normal");
        //tests .params files
        xlSearchTask->RunTask(folderPath, dbVec, FileVec, "normal");
        
        auto baseDir = FileSystem::getDirectoryName(db->getFilePath());
#ifdef ORIG
        auto directory = new DirectoryInfo(baseDir);
        std::vector<DirectoryInfo*> directories = directory->GetDirectories();
        for (auto possibleFolder : directories)
        {
            if (FileSystem::fileExists(possibleFolder->FullName, "indexEngine.params"))
            {
                File::Delete(possibleFolder->GetFiles().ElementAt(0)->FullName);
            }
        }
#endif
        for ( auto entry: std::filesystem::directory_iterator(baseDir) ) {
            if ( entry.is_directory() ) {
                if ( std::filesystem::exists(entry.path().string() + "/indexEngine.params") ) {
                    std::filesystem::remove(entry.path().string() + "/indexEngine.params");
                }
            }
        }
        
        //tests without .params files
        xlSearchTask->RunTask(folderPath, dbVec, FileVec, "normal");
        
        //auto lines = File::ReadAllLines(folderPath + "/XL_Intralinks.tsv");
        std::ifstream input(folderPath + "/XL_Intralinks.tsv");
        int lines=0;
        if ( input.is_open() ) {
            std::string line;
            while ( getline(input, line ) ) {
                lines++;
            }
        }
        else {
            std::cout << "Could not open file " << folderPath << "/XL_Intralinks.tsv" << std::endl;
        }
        
        Assert::IsTrue(lines == 2);

        std::filesystem::remove_all(folderPath);
        std::filesystem::remove_all(testdir+"/Task Settings");
        
        delete db;
        delete xlSearchTask;
    }
    
    void XLTest::TestGetPossibleCrosslinkerSites()
    {
        std::unordered_map<std::string, Modification*> mMods;
        PeptideWithSetModifications *peptide = new PeptideWithSetModifications("PEPTIDE", mMods);
        std::vector<char> vChar = {'P'};
        std::vector<int> sites = CrosslinkSpectralMatch::GetPossibleCrosslinkerModSites( vChar, peptide);

        std::vector<int> vInt = {1, 3};
        Assert::SequenceEqual(sites, vInt );
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete peptide' statement was not added since peptide was
        //passed to a method or constructor. Handle memory management manually.
    }

    void XLTest::TestTheoreticalFragmentsLoop()
    {
        Protein *p = new Protein("PEPTIDE", "");
        DigestionParams tempVar("trypsin");
        
        std::vector<Modification*> vm1, vm2;
        auto peptide = p->Digest(&tempVar, vm1, vm2).front();
        
        ModificationMotif *motif;
        ModificationMotif::TryGetMotif("X", &motif);
#ifdef ORIG
        auto loopMod = new Modification("Loop", _modificationType: "XLTest", _target: motif, _locationRestriction: "Anywhere.",
                                        _monoisotopicMass: 10000);
#endif
        auto loopMod = new Modification("Loop", "", "XLTest", "", motif, "Anywhere.",nullptr, 
                                        std::make_optional((double)10000));
        std::vector<int> tVec = {3, 5};
        auto loopLocationsWithFragments = CrosslinkedPeptide::XlLoopGetTheoreticalFragments(DissociationType::HCD, loopMod,
                                                                                            tVec, peptide);
        
        Assert::IsTrue((int)loopLocationsWithFragments.size() == 1);
        auto loopLocationWithFragments = *(loopLocationsWithFragments.begin());
        
#ifdef ORIG
        Assert::IsTrue(loopLocationWithFragments->Key->Item1 == 3);
        Assert::IsTrue(loopLocationWithFragments->Key->Item2 == 5);
#endif
        Assert::IsTrue(std::get<0>(std::get<0>(loopLocationWithFragments)) == 3);
        Assert::IsTrue(std::get<1>(std::get<0>(loopLocationWithFragments)) == 5);
        
        auto fragments = std::get<1>(loopLocationWithFragments);
        
#ifdef ORIG
        auto bIons = fragments->Where([&] (std::any v)  {
                return v->ProductType == ProductType::b;
            }).ToList();
        Assert::IsTrue(bIons.Select([&] (std::any v)   {
                    (int)v.NeutralMass;
                }).SequenceEqual(std::vector<int> {97, 226, 10537, 10652}));
#endif
        std::vector<int> bIons;
        for ( auto v : fragments) {
            if ( v->productType == ProductType::b ) {
                bIons.push_back((int)v->NeutralMass);
            }
        }
        std::vector<int> iVec = {97, 226, 10537, 10652};
        Assert::SequenceEqual(bIons, iVec);
        
#ifdef ORIG
        auto yIons = fragments->Where([&] (std::any v)   {
                return v->ProductType == ProductType::y;
            }).ToList();
        Assert::IsTrue(yIons.Select([&] (std::any v)      {
                    (int)v.NeutralMass;
                }).SequenceEqual(std::vector<int> {147, 262, 10573, 10702}));
#endif
        std::vector<int> yIons;
        for ( auto v: fragments  ) {
            if ( v->productType == ProductType::y ) {
                yIons.push_back((int)v->NeutralMass);
            }
        }
        std::vector<int> iVec2 = {147, 262, 10573, 10702};
        Assert::SequenceEqual(yIons, iVec2);
        
        delete loopMod;
        delete p;
    }
    
void XLTest::TestTheoreticalLoopFragmentsWithMod()
    {
        ModificationMotif *tMotif;
        ModificationMotif::TryGetMotif("T", &tMotif);

#ifdef ORIG
        Modification *phospho = new Modification(_originalId: "Phospho", _modificationType: "Mod", _locationRestriction: "Anywhere.",
                                                 _monoisotopicMass: 79.98, _target: tMotif);
#endif
        
        Modification *phospho = new Modification("Phospho", "", "Mod", "", tMotif, "Anywhere.", nullptr,
                                                 std::make_optional(79.98) );
        
        std::vector<int> modPositions = {2, 3, 4, 5, 6};
        
        for (auto modPosition : modPositions)
        {
            std::unordered_map<int, std::vector<Modification*>> oneBasedMods =
                {
                    {
                        modPosition, {phospho}
                    }
                };
            DigestionParams tempVar("trypsin");
#ifdef ORIG
            Protein *p = new Protein("PTTTTTE", "", oneBasedModifications: oneBasedMods);
            auto peptide = p->Digest(&tempVar, std::vector<Modification*>(), std::vector<Modification*>()).Where([&] (std::any v){
                    return v::AllModsOneIsNterminus->Count == 1;
                }).First();
#endif
            std::vector<std::tuple<std::string, std::string>> geneNames;
            Protein *p = new Protein("PTTTTTE", "", "", geneNames, oneBasedMods);

            std::vector<Modification*> vm1, vm2;
            auto peptideVec =  p->Digest(&tempVar, vm1, vm2 );
            PeptideWithSetModifications *peptide;
            for ( auto v : peptideVec ) {
                if ( (int) (v->getAllModsOneIsNterminus().size()) == 1 ){
                    peptide = v;
                    break;
                }
            }
            
            ModificationMotif *motif;
            ModificationMotif::TryGetMotif("X", &motif);
#ifdef ORIG
            auto loopMod = new Modification("Loop", _modificationType: "XLTest", _target: motif, _locationRestriction: "Anywhere.",
                                            _monoisotopicMass: 10000);
#endif
            auto loopMod = new Modification("Loop", "", "XLTest", "", motif, "Anywhere.", nullptr, 
                                            std::make_optional((double)10000));

            std::vector<int> tVec = {3, 5};
            auto loopLocationsWithFragments = CrosslinkedPeptide::XlLoopGetTheoreticalFragments(DissociationType::HCD, loopMod,
                                                                                                tVec, peptide);
            
            Assert::IsTrue((int) loopLocationsWithFragments.size() == 1);
            auto loopLocationWithFragments = *(loopLocationsWithFragments.begin());
            
            Assert::IsTrue(std::get<0>(std::get<0>(loopLocationWithFragments)) == 3);
            Assert::IsTrue(std::get<1>(std::get<0>(loopLocationWithFragments)) == 5);
            
            auto fragments = std::get<1>(loopLocationWithFragments);
            
#ifdef ORIG
            auto bIons = fragments->Where([&] (std::any v)   {
                    return v->ProductType == ProductType::b;
                }).ToList();
            auto yIons = fragments->Where([&] (std::any v)   {
                    return v->ProductType == ProductType::y;
                }).ToList();
#endif
            std::vector<int> bIons;
            for ( auto v : fragments) {
                if ( v->productType == ProductType::b ) {
                    bIons.push_back((int)v->NeutralMass);
                }
            }
            
            std::vector<int> yIons;
            for ( auto v: fragments  ) {
                if ( v->productType == ProductType::y ) {
                    yIons.push_back((int)v->NeutralMass);
                }
            }
            
            if (modPosition == 2)
            {
                //             _
                //            | |
                // PT[Phospho]TTTTE
#ifdef ORIG
                Assert::IsTrue(bIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {97, 278, 10581, 10682}));
                Assert::IsTrue(yIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {147, 248, 10551, 10732}));
#endif
                std::vector<int> vInt1 = {97, 278, 10581, 10682};
                std::vector<int> vInt2 = {147, 248, 10551, 10732};
                Assert::SequenceEqual(bIons, vInt1);
                Assert::SequenceEqual(yIons, vInt2);
                
            }
            else if (modPosition == 3)
            {
                //    __________
                //   |          |
                // PTT[Phospho]TTTE
#ifdef ORIG
                Assert::IsTrue(bIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {97, 198, 10581, 10682}));
                Assert::IsTrue(yIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {147, 248, 10631, 10732}));
#endif
                std::vector<int> vInt1 = {97, 198, 10581, 10682};
                std::vector<int> vInt2 = {147, 248, 10631, 10732};
                Assert::SequenceEqual(bIons, vInt1);
                Assert::SequenceEqual(yIons, vInt2);
            }
            else if (modPosition == 4)
            {
                //    __________
                //   |          |
                // PTTT[Phospho]TTE
#ifdef ORIG
                Assert::IsTrue(bIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {97, 198, 10581, 10682}));
                Assert::IsTrue(yIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {147, 248, 10631, 10732}));
#endif
                std::vector<int> vInt1 = {97, 198, 10581, 10682};
                std::vector<int> vInt2 = {147, 248, 10631, 10732};
                Assert::SequenceEqual(bIons, vInt1);
                Assert::SequenceEqual(yIons, vInt2);
            }
            else if (modPosition == 5)
            {
                //    _
                //   | |
                // PTTTT[Phospho]TE
#ifdef ORIG
                Assert::IsTrue(bIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {97, 198, 10581, 10682}));
                Assert::IsTrue(yIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {147, 248, 10631, 10732}));
#endif
                std::vector<int> vInt1 = {97, 198, 10581, 10682};
                std::vector<int> vInt2 = {147, 248, 10631, 10732};
                Assert::SequenceEqual(bIons, vInt1);
                Assert::SequenceEqual(yIons, vInt2);
            }
            else if (modPosition == 6)
            {
                //    _
                //   | |
                // PTTTTT[Phospho]E
#ifdef ORIG
                Assert::IsTrue(bIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {97, 198, 10501, 10682}));
                Assert::IsTrue(yIons.Select([&] (std::any v) {
                            (int)v.NeutralMass;
                        }).SequenceEqual(std::vector<int> {147, 328, 10631, 10732}));
#endif
                std::vector<int> vInt1 = {97, 198, 10501, 10682};
                std::vector<int> vInt2 = {147, 328, 10631, 10732};
                Assert::SequenceEqual(bIons, vInt1);
                Assert::SequenceEqual(yIons, vInt2);
            }
            
            delete loopMod;
            delete p;
        }
        
        delete phospho;
    }
    
#ifdef LATER
    void XLTest::TestDeadendTris()
    {
        Protein *protein = new Protein("PEPTIDE", "");
        auto csms = std::vector<CrosslinkSpectralMatch*>(1);
        
        // generate the scan with the deadend mod peptide's fragments
        auto scans = std::vector<Ms2ScanWithSpecificMass*>(1);
        ModificationMotif *motif;
        ModificationMotif::TryGetMotif("T", &motif);
        auto crosslinker = new Crosslinker("T", "T", "test", false, 100, 0, 0, 0, 0, 0, 50);
        Modification *deadend = new Modification("TestId", _target: motif, _locationRestriction: "Anywhere.",
                                                 _monoisotopicMass: crosslinker->getDeadendMassTris(), _modificationType: "Test");
        
        DigestionParams tempVar();
        auto deadendPeptide = protein->Digest(&tempVar, std::vector<Modification*> {deadend}, std::vector<Modification*>()).First();
        
        std::vector<double> mz = deadendPeptide->Fragment(DissociationType::HCD, FragmentationTerminus::Both)->Select([&] (std::any p) {
                p::NeutralMass::ToMz(1);
            }).OrderBy([&] (std::any v) {
                    return v;
                })->ToArray();
        std::vector<double> intensities = {1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        
        MzSpectrum *spectrum = new MzSpectrum(mz, intensities, false);
        MsDataScan *sc = new MsDataScan(spectrum, 1, 2, true, Polarity::Positive, 1, spectrum->Range, "", MZAnalyzerType::Orbitrap,
                                        12, 1.0, nullptr, nullptr);
        CommonParameters tempVar2();
        scans[0] = new Ms2ScanWithSpecificMass(sc, Chemistry::ClassExtensions::ToMz(deadendPeptide->MonoisotopicMass,2), 2, "", &tempVar2);
        
        // search the data with the peptide WITHOUT the deadend mod annotated in the search database.
        // the search engine should be able to correctly identify the deadend mod on T
        IndexingEngine tempVar3({protein}, new std::vector<Modification*>(), new std::vector<Modification*>(), 0, DecoyType::None,
                                new CommonParameters(), 1000, false, new std::vector<FileInfo*>(), new std::vector<std::string>());
        auto indexingResults = static_cast<IndexingResults*>((&tempVar3)->Run());
        
        CrosslinkSearchEngine tempVar4(csms, scans, indexingResults->getPeptideIndex(), indexingResults->getFragmentIndex(), 0,
                                       new CommonParameters(), crosslinker, false, 0, false, false, true, new std::vector<std::string>());
        (&tempVar4)->Run();
        
        CrosslinkSpectralMatch *csm = csms.First();
        Assert::IsTrue(csm->getCrossType() == PsmCrossType::DeadEndTris);
        Assert::IsTrue(csm->MatchedFragmentIons->Count == 12);
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete sc' statement was not added since sc was
        //passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete spectrum' statement was not added since
        //spectrum was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete deadend' statement was not added since
        //deadend was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since
        //crosslinker was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protein' statement was not added since
        //protein was passed to a method or constructor. Handle memory management manually.
    }
    
    void XLTest::TestTheoreticalFragmentsNonCleavableCrosslink()
    {
        Crosslinker *c = new Crosslinker("P", "R", "Test", false, 1000, 0, 0, 1000, 5, 5, 5);
        
        Protein *p1 = new Protein("PEPTIDE", "");
        DigestionParams tempVar();
        auto alphaPeptide = p1->Digest(&tempVar, std::vector<Modification*>(), std::vector<Modification*>()).First();
        
        Protein *p2 = new Protein("PRLTEIN", "");
        DigestionParams tempVar2();
        auto betaPeptide = p2->Digest(&tempVar2, std::vector<Modification*>(), std::vector<Modification*>()).Where([&] (std::any v){
                return v->MissedCleavages == 1;
            }).First();
        
        auto theoreticalCrosslinkFragments = CrosslinkedPeptide::XlGetTheoreticalFragments(DissociationType::HCD, c, std::vector<int> {3},
                                                                                           betaPeptide->MonoisotopicMass, alphaPeptide).ToList();
        
        Assert::IsTrue(theoreticalCrosslinkFragments.size() == 1);
        auto loopLocationWithFragments = theoreticalCrosslinkFragments.front();
        
        Assert::IsTrue(loopLocationWithFragments->Item1 == 3);
        
        auto fragments = loopLocationWithFragments->Item2;
        
        auto bIons = fragments->Where([&] (std::any v)  {
                return v->ProductType == ProductType::b;
            }).ToList();
        Assert::IsTrue(bIons.Select([&] (std::any v)  {
                    (int)v.NeutralMass;
                }).SequenceEqual(std::vector<int> {97, 226, 2164, 2265, 2378, 2493}));
        
        auto yIons = fragments->Where([&] (std::any v)  {
                return v->ProductType == ProductType::y;
            }).ToList();
        Assert::IsTrue(yIons.Select([&] (std::any v) {
                    (int)v.NeutralMass;
                }).SequenceEqual(std::vector<int> {147, 262, 375, 476, 2414, 2543}));
        
        delete p2;
        delete p1;
        //C# TO C++ CONVERTER TODO TASK: A 'delete c' statement was not added since c was
        //passed to a method or constructor. Handle memory management manually.
    }
    
    void XLTest::TestTheoreticalFragmentsCleavableCrosslink()
    {
        Crosslinker *c = new Crosslinker("P", "R", "Test", true, 1000, 15, 25, 1000, 5, 5, 5);
        
        Protein *p1 = new Protein("PEPTIDE", "");
        DigestionParams tempVar();
        auto alphaPeptide = p1->Digest(&tempVar, std::vector<Modification*>(), std::vector<Modification*>()).First();
        
        Protein *p2 = new Protein("PRLTEIN", "");
        DigestionParams tempVar2();
        auto betaPeptide = p2->Digest(&tempVar2, std::vector<Modification*>(), std::vector<Modification*>()).Where([&] (std::any v) {
                return v->MissedCleavages == 1;
            }).First();
        
        auto theoreticalCrosslinkFragments = CrosslinkedPeptide::XlGetTheoreticalFragments(DissociationType::HCD, c, std::vector<int> {3},
                                                                                           10000, alphaPeptide).ToList();
        
        Assert::IsTrue(theoreticalCrosslinkFragments.size() == 1);
        
        // cleaved fragments
        auto linkLocationWithFragments = theoreticalCrosslinkFragments[0];
        
        Assert::IsTrue(linkLocationWithFragments.Item1 == 3);
        auto fragmentsWithCleavedXlPieces = linkLocationWithFragments.Item2;
        
        auto bIons = fragmentsWithCleavedXlPieces->Where([&] (std::any v) {
                return v->ProductType == ProductType::b;
            }).ToList();
        Assert::IsTrue(bIons.Select([&] (std::any v) {
                    (int)v.NeutralMass;
                }).SequenceEqual(std::vector<int> {97, 226, 338, 439, 552, 667, 348, 449, 562, 677}));
        
        auto yIons = fragmentsWithCleavedXlPieces->Where([&] (std::any v) {
                return v->ProductType == ProductType::y;
            }).ToList();
        Assert::IsTrue(yIons.Select([&] (std::any v)  {
                    (int)v.NeutralMass;
                }).SequenceEqual(std::vector<int> {147, 262, 375, 476, 588, 717, 598, 727}));
        
        auto signatureIons = fragmentsWithCleavedXlPieces->Where([&] (std::any v)  {
                return v->ProductType == ProductType::M;
            }).ToList();
        Assert::IsTrue(signatureIons.size() == 2);
        Assert::IsTrue(signatureIons.Select([&] (std::any v) {
                    (int)v.NeutralMass;
                }).SequenceEqual(std::vector<int> {814, 824}));
        
        delete p2;
        delete p1;
        //C# TO C++ CONVERTER TODO TASK: A 'delete c' statement was not added since c was
        //passed to a method or constructor. Handle memory management manually.
    }
    
    void XLTest::TestWriteToPercolator()
    {
        std::string testdir=std::filesystem::current_path().string();

        XLSearchTask *xlst = new XLSearchTask();
        XlSearchParameters tempVar();
        xlst->setXlSearchParameters(&tempVar);
        xlst->getXlSearchParameters()->setWriteOutputForPercolator(true);
        
        std::string myFileXl = testdir + "/XlTestData/BSA_DSSO_ETchD6010.mgf";
        std::string myDatabaseXl = testdir + "/XlTestData/BSA.fasta";
        std::string outputFolder = testdir + "/TestXLSearch";
        DbForTask *db = new DbForTask(myDatabaseXl, false);
        
        std::vector<(std::string, MetaMorpheusTask)*> taskList = {("TestPercolator", xlst)};
        
        auto engine = new EverythingRunnerEngine(taskList, std::vector<std::string> {myFileXl}, std::vector<DbForTask*> {db}, outputFolder);
        engine->Run();
        
        auto results = outputFolder, "(TestPercolator\XL_Intralinks_Percolator.txt)");
        auto lines = File::ReadAllLines(results);
        Assert::IsTrue(lines[0] == "SpecId\tLabel\tScannr\tScore\tdScore\tNormRank\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProtein");
        Assert::IsTrue(lines[1] == "T-1-30.6190992666667\t1\t1\t20.0641008915522\t0\t7\t3\t1994.05202313843\t0.664979354397676\t7\t9\t16\t-.EKVLTSSAR2--LSQKFPK4.-\t3336842(211)\t3336842(245)");
        Directory::Delete(outputFolder, true);
        
        delete engine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete db' statement was not added since db was passed
        //to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete xlst' statement was not added since xlst was passed
        //to a method or constructor. Handle memory management manually.
    }
    
    void XLTest::TestWriteNonSingleCross()
    {
        std::string testdir=std::filesystem::current_path().string();

        XLSearchTask *xlst = new XLSearchTask();
        Protein *protein = new Protein("PEPTIDE", "");
        auto csms = std::vector<CrosslinkSpectralMatch*>(1);
        
        // generate the scan with the deadend mod peptide's fragments
        auto scans = std::vector<Ms2ScanWithSpecificMass*>(1);
        ModificationMotif *motif;
        ModificationMotif::TryGetMotif("T", &motif);
        auto crosslinker = new Crosslinker("T", "T", "test", false, 100, 0, 0, 0, 0, 0, 50);
        Modification *deadend = new Modification("TestId", _target: motif, _locationRestriction: "Anywhere.",
                                                 _monoisotopicMass: crosslinker->getDeadendMassTris(), _modificationType: "Test");
        
        DigestionParams tempVar();
        auto deadendPeptide = protein->Digest(&tempVar, std::vector<Modification*> {deadend}, std::vector<Modification*>()).First();
        
        std::vector<double> mz = deadendPeptide->Fragment(DissociationType::HCD, FragmentationTerminus::Both)->Select([&] (std::any p)  {
                p::NeutralMass::ToMz(1);
            }).OrderBy([&] (std::any v)  {
                    return v;
                })->ToArray();
        std::vector<double> intensities = {1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        
        MzSpectrum *spectrum = new MzSpectrum(mz, intensities, false);
        MsDataScan *sc = new MsDataScan(spectrum, 1, 2, true, Polarity::Positive, 1, spectrum->Range, "", MZAnalyzerType::Orbitrap,
                                        12, 1.0, nullptr, nullptr);
        CommonParameters tempVar2();
        scans[0] = new Ms2ScanWithSpecificMass(sc, Chemistry::ClassExtensions::ToMz(deadendPeptide->MonoisotopicMass,2), 2, "", &tempVar2);
        
        IndexingEngine tempVar3({protein}, new std::vector<Modification*>(), new std::vector<Modification*>(), 0, DecoyType::None,
                                new CommonParameters(), 1000, false, new std::vector<FileInfo*>(), new std::vector<std::string>());
        auto indexingResults = static_cast<IndexingResults*>((&tempVar3)->Run());
        
        CrosslinkSearchEngine tempVar4(csms, scans, indexingResults->getPeptideIndex(), indexingResults->getFragmentIndex(), 0,
                                       new CommonParameters(), crosslinker, false, 0, false, false, true, new std::vector<std::string>());
        (&tempVar4)->Run();
        
        csms[0]->SetFdrValues(0, 0, 0.1, 0, 0, 0, 0, 0, 0, false);
        
        xlst->WritePepXML_xl(csms.ToList(), std::vector<Protein*>(), "", {deadend}, {deadend}, std::vector<std::string>(),
                             testdir, "test", std::vector<std::string>());
        File::Delete(testdir + "/test.pep.XML");
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete sc' statement was not added since sc was passed
        //to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete spectrum' statement was not added since
        //spectrum was passed to a method or constructor. Handle memory management manually.
         //C# TO C++ CONVERTER TODO TASK: A 'delete deadend' statement was not added since
        //deadend was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added
        //since crosslinker was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protein' statement was not added since protein
        //was passed to a method or constructor. Handle memory management manually.
        delete xlst;
    }
#endif
    
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
    
    XLTestDataFileDiffSite::XLTestDataFileDiffSite() : MsDataFile(2, new MassSpectrometry::SourceFile("", "", "", "", "" ))
    {
        auto mz1 = std::vector<double> {100, Chemistry::ClassExtensions::ToMz(1030.5956,1)};
        auto intensities1 = std::vector<double> {100, 1};
        auto MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
        auto ScansHere = new std::vector<MsDataScan*>();
        MzLibUtil::MzRange tempVar2(0, 10000);
        auto s = new MsDataScan(MassSpectrum1, 1, 1, true, Polarity::Positive, 1, &tempVar2, "ff", MZAnalyzerType::Unknown,
                                 1000, 1, std::vector<std::vector<double>>(), "scan=1");
        ScansHere->push_back(s);
        
        auto mz2 = std::vector<double> {100, 147.1128, 175.119, 213.1598, 246.1561, 275.1714, 757.4388, 786.4541, 819.4504,
                                        857.4912, 885.4974, 918.5189, 932.5345};
        auto intensities2 = std::vector<double> {100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        auto MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
        auto tempVar3 = new MsDataScan (MassSpectrum2, 2, 2, true, Polarity::Positive, 1.0, new MzLibUtil::MzRange(0, 10000),
                                        "f", MZAnalyzerType::Unknown, 112, 1.0, std::vector<std::vector<double>>(), "scan=2",
                                        Chemistry::ClassExtensions::ToMz(1030.5956,1), 1, 1,
                                        Chemistry::ClassExtensions::ToMz(1030.5956,1), 2, DissociationType::HCD, 1,
                                        Chemistry::ClassExtensions::ToMz(1030.5956,1));
        ScansHere->push_back(tempVar3);
        
        Scans = *ScansHere;//.ToArray();
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum2' statement was not added since
        //MassSpectrum2 was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete MassSpectrum1' statement was not added since
        //MassSpectrum1 was passed to a method or constructor. Handle memory management manually.
    }
    
    std::string XLTestDataFileDiffSite::getFilePath() const
    {
        return "XLTestDataFileDiffSite";
    }
    
    std::string XLTestDataFileDiffSite::getName() const
    {
        return "XLTestDataFileDiffSite";
    }
    
    void XLTestDataFileDiffSite::ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities)
    {
        MzSpectrum *massSpectrum = new MzSpectrum(mz, intensities, false);
        Scans[0] = new MsDataScan(massSpectrum, Scans[0]->getOneBasedScanNumber(), Scans[0]->getMsnOrder(),
                                  Scans[0]->getIsCentroid(), Scans[0]->getPolarity(), Scans[0]->getRetentionTime(),
                                  Scans[0]->getScanWindowRange(), Scans[0]->getScanFilter(),
                                  Scans[0]->getMzAnalyzer(), massSpectrum->getSumOfAllY(), Scans[0]->getInjectionTime(),
                                  std::vector<std::vector<double>>(), Scans[0]->getNativeId() );
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete massSpectrum' statement was not added since massSpectrum was
        //passed to a method or constructor. Handle memory management manually.
    }
}
