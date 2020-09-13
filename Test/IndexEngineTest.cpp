#include "IndexEngineTest.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Indexing/IndexingEngine.h"
#include "../EngineLayer/Indexing/IndexingResults.h"
#include "Proteomics/Proteomics.h"
#include "MassSpectrometry/MassSpectrometry.h"
#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
#include "MzLibUtil.h"
#include "Chemistry/Chemistry.h"

using namespace EngineLayer;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

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
    Chemistry::PeriodicTable::Load (elr);
    //UsefulProteomicsDatabases::PeriodicTableLoader::Load (elr);

    std::cout << ++i << ". TestIndexEngine" << std::endl;
    Test::IndexEngineTest::TestIndexEngine();

    std::cout << ++i << ". TestIndexEngineWithWeirdSeq" << std::endl;
    Test::IndexEngineTest::TestIndexEngineWithWeirdSeq();

    return 0;
}

namespace Test
{ 
    void IndexEngineTest::TestIndexEngine()
    {
        auto proteinList = std::vector<Protein*> {new Protein("MNNNKQQQ", "")};
        std::vector<Modification*> variableModifications;
        std::vector<Modification*>fixedModifications;
        std::vector<Modification*> localizeableModifications;
        
        std::unordered_map<Modification*, unsigned short> modsDictionary;
        for (auto mod : fixedModifications)
        {
            modsDictionary.emplace(mod, 0);
        }
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
        
        std::vector<DigestionMotif*> motifs = {new DigestionMotif("K", "", 1, "")};
        Protease *p = new Protease("Custom Protease2",CleavageSpecificity::Full, "", "", motifs);
        //ProteaseDictionary::Dictionary->Add(p->Name, p);
        ProteaseDictionary::insert(p->getName(), p);
        
        auto tempVar = new DigestionParams(p->getName(), 2, 1);
        //CommonParameters *CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: &tempVar);
        auto cp = new CommonParameters("",  DissociationType::HCD, true, true, 3, 12, true, false, 1, 1);
        cp->setDigestionParams(tempVar);

        std::vector<std::string> vs1, vs2;
        auto engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse,
                                         cp, 30000, false, vs1, vs2);
        
        auto results = static_cast<IndexingResults*>(engine->Run());
        
        Assert::AreEqual(5, (int)results->getPeptideIndex().size());
        std::vector<Modification*> vm1;
        auto digestedList = proteinList[0]->Digest(cp->getDigestionParams(), vm1,
                                                   variableModifications);//.ToList();
        
        Assert::AreEqual(5, (int)digestedList.size());
        auto pepIndex = results->getPeptideIndex();
        for (auto fdfd : digestedList)
        {
            //Assert->Contains(fdfd, results->getPeptideIndex());
            bool found = false;
            for ( auto p : pepIndex ) {
                if ( p->Equals(fdfd ) ){
                    found = true;
                }
            }
            Assert::IsTrue( found );
        }
        
        delete engine;
        delete cp;
        delete tempVar;
        delete p;
    }


    void IndexEngineTest::TestIndexEngineWithWeirdSeq()
    {
        auto proteinList = std::vector<Protein*> {new Protein("MQXQ", "")};
        std::vector<Modification*> variableModifications;
        std::vector<Modification*> fixedModifications;
        std::vector<Modification*> localizeableModifications;
        
        std::unordered_map<Modification*, unsigned short> modsDictionary;
        for (auto mod : fixedModifications)
        {
            modsDictionary.emplace(mod, 0);
        }
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
        
        std::vector<DigestionMotif*> motifs = {new DigestionMotif("K", "", 1, "")};
        Protease *protease = new Protease("Custom Protease", CleavageSpecificity::Full, "", "", motifs);
        ProteaseDictionary::insert(protease->getName(), protease);
        auto tempVar = new DigestionParams(protease->getName(), 2, 1, std::numeric_limits<int>::max(), 1024,
                                           InitiatorMethionineBehavior::Retain);
        //CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1);
        auto cp = new CommonParameters("",  DissociationType::HCD, true, true, 3, 12, true, false, 1, 1);
        cp->setDigestionParams(tempVar);
        
        std::vector<std::string> vs1, vs2;
        auto engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1,
                                         DecoyType::Reverse, cp, 30000, false, vs1, vs2 );
                                         
        
        auto results = static_cast<IndexingResults*>(engine->Run());
        
        Assert::AreEqual(1, (int)results->getPeptideIndex().size());
        
        //Assert::IsNaN(results->getPeptideIndex()[0]->MonoisotopicMass);
        Assert::IsTrue(std::isnan(results->getPeptideIndex()[0]->getMonoisotopicMass()));
        Assert::AreEqual(30000000 + 1, (int)results->getFragmentIndex().size());
        
        delete engine;
        delete cp;
        delete protease;
        delete tempVar;
    }

}
