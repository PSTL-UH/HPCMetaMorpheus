#include "BinGenerationTest.h"
#include "TestDataFile.h"

#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../TaskLayer/DbForTask.h"
using namespace TaskLayer;

#include "../EngineLayer/EngineLayer.h"
using namespace EngineLayer;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
using namespace UsefulProteomicsDatabases;

#include "../MzML/Mzml.h"
#include "../MzML/MzmlMethods.h"

#include "Assert.h"
#include <experimental/filesystem>


int main ( int argc, char **argv )
{
    int i=0;
    std::cout << i << ". PeriodicTableLoader" << std::endl;
    const std::string elfile="elements.dat";
    const std::string &elr=elfile;
    Chemistry::PeriodicTable::Load (elr);
    //UsefulProteomicsDatabases::PeriodicTableLoader::Load (elr);

    std::cout << ++i << ". TestBinGeneration" << std::endl;
    Test::BinGenerationTest::TestBinGeneration();

    std::cout << ++i << ". TestProteinSplitAcrossFiles" << std::endl;
    Test::BinGenerationTest::TestProteinSplitAcrossFiles();
    
    return 0;
}

namespace Test
{

    void BinGenerationTest::TestBinGeneration()
    {
        std::string testdir=std::experimental::filesystem::current_path().string();
        
        SearchTask *st = new SearchTask();
#ifdef ORIG
        CommonParameters tempVar(, , true, true, 3, 12, true, false, 1, 1, , , , , , , , , , ,
                                 new DigestionParams(minPeptideLength: 5,
                                                     initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain));
#endif
        CommonParameters* tempVar = new CommonParameters( "",  DissociationType::HCD, true, true, 3, 12, true, false, 1, 1,
                                                          200, 0.01, false,  true, false, false, nullptr, nullptr, nullptr, -1,
                                                          new DigestionParams("trypsin", 2, 5,  std::numeric_limits<int>::max(),
                                                                              1024, InitiatorMethionineBehavior::Retain));
        st->setCommonParameters(tempVar);
        SearchParameters *tempVar2 = new SearchParameters();
        st->setSearchParameters(tempVar2);
        st->getSearchParameters()->setDoHistogramAnalysis(true);
        st->getSearchParameters()->setMassDiffAcceptorType(MassDiffAcceptorType::Open);
        st->getSearchParameters()->setDecoyType(DecoyType::None);
        st->getSearchParameters()->setDoParsimony(true);
        st->getSearchParameters()->setDoQuantification(true);
        
        std::string proteinDbFilePath = testdir +  "/BinGenerationTest.xml";
        //std::string proteinDbFilePath = testdir +  "/BinGenerationTest.fasta";
        std::string mzmlFilePath = testdir + "/BinGenerationTest.mzML";
        
        Protein *prot1 = new Protein("MEDEEK", "prot1");
        Protein *prot2 = new Protein("MENEEK", "prot2");
        
        ModificationMotif *motif;
        ModificationMotif::TryGetMotif("D", &motif);
#ifdef ORIG        
        Modification *mod = new Modification(_target: motif, _locationRestriction: "Anywhere.",
                                             _monoisotopicMass: 10);
#endif
        Modification *mod = new Modification ( "", "", "", "", motif, "Anywhere.", nullptr, std::make_optional(10));

        std::vector<Modification*> tvm1, tvm2, tvm3, tvm4;
        auto pep1_0 = prot1->Digest(st->getCommonParameters()->getDigestionParams(),
                                    tvm1, tvm2).front();
        auto pep1_10 = prot1->Digest(st->getCommonParameters()->getDigestionParams(),
                                     tvm3, tvm4).back();
        
        Protein *prot3 = new Protein("MAAADAAAAAAAAAAAAAAA", "prot3");

        std::vector<Modification*> tvm5, tvm6, tvm7, tvm8 = {mod};
        auto pep2_0 = prot3->Digest(st->getCommonParameters()->getDigestionParams(),
                                    tvm5, tvm6).front();
        auto pep2_10 = prot3->Digest(st->getCommonParameters()->getDigestionParams(),
                                     tvm7, tvm8).back();
        
        Protein *prot4 = new Protein("MNNDNNNN", "prot4");
        std::vector<Modification*> tvm9, tvm10 = {mod};
        auto pep3_10 = prot4->Digest(st->getCommonParameters()->getDigestionParams(),
                                     tvm9, tvm10).back();
        
        std::vector<PeptideWithSetModifications*> pepsWithSetMods = {pep1_0, pep1_10, pep2_0, pep2_10, pep3_10};
        MsDataFile *myMsDataFile = new TestDataFile(pepsWithSetMods);
        
        std::vector<Protein*> proteinList = {prot1, prot2, prot3, prot4};
        
        IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlFilePath, false);
        std::unordered_map<std::string, ModDbTuple_set> tempvar;
        ProteinDbWriter::WriteXmlDatabase(tempvar, proteinList, proteinDbFilePath);
        //ProteinDbWriter::WriteFastaDatabase(proteinList, proteinDbFilePath, "|");
        
        std::string output_folder = testdir + "/TestBinGeneration";
        FileSystem::createDirectory(output_folder);
        auto db  = new DbForTask(proteinDbFilePath, false);
        std::vector<TaskLayer::DbForTask*> tdb = {db};
        std::vector<std::string> vs = {mzmlFilePath};
        st->RunTask(output_folder, tdb , vs, "");
        
#ifdef ORIG
        Assert::AreEqual(3, File::ReadLines(testdir + "/MassDifferenceHistogram.tsv").size());
#endif
        int count=0;
        std::ifstream input(output_folder + "/MassDifferenceHistogram.tsv");
        if ( input.is_open() ) {
            std::string line;
            while ( getline(input, line ) ) {
                count++;
            }
        }
        else {
            std::cout << "Could not open file " << output_folder << "/MassDifferenceHistogram.tsv" << std::endl;
        }
        Assert::AreEqual(3, count);
        
        std::experimental::filesystem::remove(proteinDbFilePath);
        std::experimental::filesystem::remove(mzmlFilePath);
        std::experimental::filesystem::remove_all(output_folder);
        std::experimental::filesystem::remove_all(testdir + "/Task Settings");        
        delete myMsDataFile;
        delete prot4;
        delete prot3;
        delete mod;
        delete prot2;
        delete prot1;
        delete st;
    }
    
    void BinGenerationTest::TestProteinSplitAcrossFiles()
    {
        std::string testdir=std::experimental::filesystem::current_path().string();
        
        SearchTask *st = new SearchTask();
#ifdef ORIG
        CommonParameters tempVar(, , true, true, 3, 12, true, false, 1, 1, , , , , , , , , , ,
                                 new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5,
                                                     initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain));
#endif
        CommonParameters *tempVar = new CommonParameters( "",  DissociationType::HCD, true, true, 3, 12, true, false, 1, 1,
                                                          200, 0.01, false,  true, false, false, nullptr, nullptr, nullptr, -1,
                                                          new DigestionParams("trypsin", 0, 5,  std::numeric_limits<int>::max(),
                                                                              1024, InitiatorMethionineBehavior::Retain));
        
        st->setCommonParameters(tempVar);
        SearchParameters *tempVar2  = new SearchParameters();
        st->setSearchParameters(tempVar2);
        st->getSearchParameters()->setDoHistogramAnalysis(true);
        st->getSearchParameters()->setMassDiffAcceptorType(MassDiffAcceptorType::Open);
        st->getSearchParameters()->setMatchBetweenRuns(true);
        st->getSearchParameters()->setDoQuantification(true);
        
        std::string proteinDbFilePath = testdir + "/TestProteinSplitAcrossFiles.xml";
        std::string mzmlFilePath1 = testdir + "/TestProteinSplitAcrossFiles1.mzML";
        std::string mzmlFilePath2 = testdir + "/TestProteinSplitAcrossFiles2.mzML";
        
        ModificationMotif *motif;
        ModificationMotif::TryGetMotif("D", &motif);

#ifdef ORIG
        Modification *mod = new Modification(_originalId: "mod1 on D", _modificationType: "mt", _target: motif,
                                             _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
#endif
        Modification *mod = new Modification("mod1 on D", "", "mt", "", motif, "Anywhere.", nullptr, std::make_optional(10));
        std::unordered_map<int, std::vector<Modification*>> oneBasedModification =
            {
                {
                    3, {mod}
                }
            };
        
        Protein *prot1 = new Protein("MEDEEK", "prot1", "", std::vector<std::tuple<std::string, std::string>>(),
                                     oneBasedModification);
        
        std::vector<Modification*> tvm1, tvm2, tvm3, tvm4;
        auto pep1 = prot1->Digest(st->getCommonParameters()->getDigestionParams(), 
                                  tvm1, tvm2).front();
        auto pep2 = prot1->Digest(st->getCommonParameters()->getDigestionParams(), 
                                  tvm3, tvm4).back();
        
        std::vector<PeptideWithSetModifications*> listForFile1 = {pep1, pep2};
        std::vector<PeptideWithSetModifications*> listForFile2 = {pep2};
        MsDataFile *myMsDataFile1 = new TestDataFile(listForFile1);
        MsDataFile *myMsDataFile2 = new TestDataFile(listForFile2);
        
        std::vector<Protein*> proteinList = {prot1};
        
        IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlFilePath1, false);
        IO::MzML::MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile2, mzmlFilePath2, false);
        std::unordered_map<std::string, ModDbTuple_set> tempvar;
        ProteinDbWriter::WriteXmlDatabase(tempvar, proteinList, proteinDbFilePath);
        
        std::string output_folder = testdir + "/TestProteinSplitAcrossFiles";
        FileSystem::createDirectory(output_folder);

        auto db  = new DbForTask(proteinDbFilePath, false);
        std::vector<TaskLayer::DbForTask*> tdb = {db};
        std::vector<std::string> vs = {mzmlFilePath1, mzmlFilePath2};

        st->RunTask(output_folder, tdb, vs, "");

        //std::experimental::filesystem::remove(proteinDbFilePath);
        //std::experimental::filesystem::remove(mzmlFilePath1);
        //std::experimental::filesystem::remove(mzmlFilePath2);
        //std::experimental::filesystem::remove_all(output_folder);
        //std::experimental::filesystem::remove_all(testdir + "/Task Settings");
        
        delete myMsDataFile2;
        delete myMsDataFile1;
        delete prot1;
        delete mod;
        delete st;
    }
}
