#include "IndexEngineTest.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Indexing/IndexingEngine.h"
#include "Proteomics/Proteomics.h"
#include "MassSpectrometry/MassSpectrometry.h"

using namespace EngineLayer;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
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

    std::cout << ++i << ". TestIndexEngine" << std::endl;
    Test::IndexEngineTest::TestIndexEngine();

#ifdef LATER
    std::cout << ++i << ". TestIndexEngineWithWeirdSeq" << std::endl;
    Test::IndexEngineTest::TestIndexEngineWithWeirdSeq();
#endif

    return 0;
}

namespace Test
{ 
    void IndexEngineTest::TestIndexEngine()
    {
        auto proteinList = std::vector<Protein*> {new Protein("MNNNKQQQ", nullptr)};
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto localizeableModifications = std::vector<Modification*>();
        
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
        
        std::vector<DigestionMotif*> motifs = {new DigestionMotif("K", nullptr, 1, nullptr)};
        Protease *p = new Protease("Custom Protease2",CleavageSpecificity::Full, nullptr, nullptr, motifs);
        ProteaseDictionary::Dictionary->Add(p->Name, p);

        DigestionParams tempVar(protease: p->Name, minPeptideLength: 1);
        CommonParameters *CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: &tempVar);

        std::vector<std::string> vs1, vs2;
        auto engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse,
                                         CommonParameters, 30000, false, vs1, vs2);
        
        auto results = static_cast<IndexingResults*>(engine->Run());
        
        Assert::AreEqual(5, (int)results->getPeptideIndex().size());
        
        auto digestedList = proteinList[0]->Digest(CommonParameters->getDigestionParams(), std::vector<Modification*>(),
                                                   variableModifications);//.ToList();
        
        Assert::AreEqual(5, (int)digestedList.size());
        for (auto fdfd : digestedList)
        {
            Assert->Contains(fdfd, results->getPeptideIndex());
        }
        
        delete engine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement
        //was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete p' statement
        //was not added since p was passed to a method or constructor. Handle memory management manually.
    }

#ifdef LATER    
    void IndexEngineTest::TestIndexEngineWithWeirdSeq()
    {
        auto proteinList = std::vector<Protein*> {new Protein("MQXQ", nullptr)};
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto localizeableModifications = std::vector<Modification*>();
        
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
        
        std::vector<DigestionMotif*> motifs = {new DigestionMotif("K", nullptr, 1, nullptr)};
        Protease *protease = new Protease("Custom Protease", CleavageSpecificity::Full, nullptr, nullptr, motifs);
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        DigestionParams tempVar(protease: protease->Name, minPeptideLength: 1,
                                initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain);
        CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1);
        
        std::vector<std::string> vs1, vs2;
        auto engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1,
                                         DecoyType::Reverse, CommonParameters, 30000, false, vs1, vs2 );
                                         
        
        auto results = static_cast<IndexingResults*>(engine->Run());
        
        Assert::AreEqual(1, (int)results->getPeptideIndex().size());
        
        Assert::IsNaN(results->getPeptideIndex()[0]->MonoisotopicMass);
        Assert::AreEqual(30000000 + 1, (int)results->getFragmentIndex().size());
        
        delete engine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement
        //was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement
        //was not added since protease was passed to a method or constructor. Handle memory management manually.
    }
#endif
}
