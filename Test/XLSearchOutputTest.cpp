#include "XLSearchOutputTest.h"
#include "../TaskLayer/XLSearchTask/XLSearchTask.h"
#include "../TaskLayer/DbForTask.h"

using namespace TaskLayer;

#include "MzLibAssert.h"
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

    std::cout << ++i << ". WriteTsvTest" << std::endl;
    Test::XLSearchOutputTest::WriteTsvTest();
}

namespace Test
{

    void XLSearchOutputTest::WriteTsvTest()
    {
        std::string testdir=std::filesystem::current_path().string();
        
        std::string outputFolder = testdir + "/XlOutputTest1";
        std::string myFile = testdir + "/XlTestData/BSA_DSS_23747.mzML";
        std::string myDatabase = testdir + "/XlTestData/BSA.fasta";
        
        FileSystem::createDirectory(outputFolder);
        
        XLSearchTask *xLSearch = new XLSearchTask();
        auto db = new DbForTask(myDatabase, false);
        std::vector<DbForTask *> dbvec = {db};
        std::vector<std::string> fvec = {myFile};
        xLSearch->RunTask(outputFolder, dbvec, fvec, "test");
        
        //auto resultsPath = File::ReadAllLines(outputFolder + "/XL_Interlinks.tsv");
        std::vector<std::string> resultsPath;
        std::ifstream fs(outputFolder + "/XL_Interlinks.tsv" ) ;
        if ( fs.is_open() ) {
            std::string line;
            while ( getline(fs, line ) ) {
                resultsPath.push_back(line);
                //std::cout << line << std::endl;
            }
        }
        else {
            std::cout << "Could not open file " << outputFolder + "/XL_Interlinks.tsv" << std::endl;
        }
        if ( resultsPath.size() > 0 ) {
            auto sections = StringHelper::split(resultsPath[1], '\t');
            Assert::AreEqual((int)sections.size(), 45);
        }
        Assert::AreEqual((int)resultsPath.size(), 2);
        std::filesystem::remove_all(outputFolder);
        
        delete xLSearch;
    }
}
