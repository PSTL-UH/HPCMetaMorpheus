#include "CalibrationTests.h"
#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
#include "../TaskLayer/DbForTask.h"
#include <regex>
#include <string>

#include "../TaskLayer/CalibrationTask/CalibrationTask.h"
using namespace TaskLayer;

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
using namespace UsefulProteomicsDatabases;

#include "MzLibAssert.h"
#include <experimental/filesystem>
#include <iostream>
#include <fstream>

int main ( int argc, char **argv )
{
    int i=0;
    std::cout << i << ". PeriodicTableLoader" << std::endl;
    const std::string elfile="elements.dat";
    const std::string &elr=elfile;
    Chemistry::PeriodicTable::Load (elr);
    //UsefulProteomicsDatabases::PeriodicTableLoader::Load (elr);

    std::cout << ++i << ". ExperimentalDesignCalibrationTest" << std::endl;
    Test::CalibrationTests::ExperimentalDesignCalibrationTest();

    return 0;
}

namespace Test
{

    void CalibrationTests::ExperimentalDesignCalibrationTest()
    {
        std::string testdir=std::experimental::filesystem::current_path().string();

        CalibrationTask *calibrationTask = new CalibrationTask();
        std::string outputFolder = testdir + "/TestCalibration";
        std::string myFile =  "TestData/SmallCalibratible_Yeast.mzML";
        std::string myDatabase =  "TestData/smalldb.fasta";
        std::string filePath =  "TestData/ExperimentalDesign.tsv";
        FileSystem::createDirectory(outputFolder);
        
        std::ofstream output(filePath);
        output << "FileName\tCondition\tBiorep\tFraction\tTechrep \n";
        output << "SmallCalibratible_Yeast" << "\t" << "condition" << "\t" << "1" << "\t" << "1" << "\t" << "1";
        output.close();
        
        auto tobj = new DbForTask(myDatabase, false);
        std::vector<DbForTask*> tvec = {tobj};
        std::vector<std::string> svec = {myFile};
        calibrationTask->RunTask(outputFolder, tvec, svec,  "test");
        auto expDesignPath = outputFolder + "/ExperimentalDesign.tsv";
        //auto expDesign = File::ReadAllLines(expDesignPath);
        std::vector<std::string> expDesign;
        std::ifstream input(expDesignPath);
        if ( input.is_open() ) {
            std::string line;
            while ( getline(input, line ) ) {
                expDesign.push_back(line);
            }
        }
        input.close();
        
        Assert::IsTrue(expDesign[1].find("SmallCalibratible_Yeast-calib") != std::string::npos);
        Assert::IsTrue(FileSystem::fileExists(outputFolder + "/SmallCalibratible_Yeast-calib.mzML"));
        Assert::IsTrue(FileSystem::fileExists(outputFolder + "/SmallCalibratible_Yeast-calib.toml"));
        //auto lines = File::ReadAllLines(outputFolder + "/SmallCalibratible_Yeast-calib.toml");
        std::vector<std::string> lines;
        std::ifstream input2(outputFolder + "/SmallCalibratible_Yeast-calib.toml");
        if (input2.is_open() ){
            std::string line;
            while (getline(input2, line) ) {
                lines.push_back(line);
            }
        }
        input2.close();
        
        //auto tolerance = Regex::Match(lines[0], "(\d+\.\d*)")->Value;
        std::regex reg1(R"(\d+\.\d*)");
        std::sregex_iterator sreg1(lines[0].begin(), lines[0].end(), reg1 ) ;
        std::sregex_iterator reg_end;
        if ( sreg1 != reg_end ) {
            std::string tolerance = sreg1->str();
            double tol;
            //Assert::IsTrue(double::TryParse(tolerance, tol) == true);
            try {
                tol = std::stod(tolerance);
            }
            catch(std::invalid_argument &e) {
                Assert::IsTrue(false);
            }
        }
        else {
            Assert::IsTrue(false);
        }
        
        //auto tolerance1 = Regex::Match(lines[1], "(\d+\.\d*)")->Value;
        std::sregex_iterator sreg2(lines[1].begin(), lines[1].end(), reg1 ) ;
        if ( sreg2 != reg_end ) {
            std::string tolerance1 = sreg2->str();
            double tol1;
            //Assert::IsTrue(double::TryParse(tolerance1, tol1) == true);
            try {
                tol1 = std::stod(tolerance1);
            }
            catch(std::invalid_argument &e) {
                Assert::IsTrue(false);
            }
        }
        else {
            Assert::IsTrue(false);
        }


        Assert::IsTrue(lines[0].find("PrecursorMassTolerance") != std::string::npos);
        Assert::IsTrue(lines[1].find("ProductMassTolerance") != std::string::npos);


        std::experimental::filesystem::remove(filePath);
        std::experimental::filesystem::remove_all(outputFolder);
        std::experimental::filesystem::remove_all("Task Settings");
        
        delete calibrationTask;
    }
}
