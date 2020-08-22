#include "SearchEngineTests.h"
#include "../EngineLayer/CommonParameters.h"
#include "TestDataFile.h"
#include "../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/PrecursorSearchModes/OpenMassDiffAcceptor.h"
#include "../TaskLayer/SearchTask/SearchParameters.h"
#include "../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
#include "../EngineLayer/ProteinScoringAndFdr/FdrCategory.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/MetaMorpheusEngineResults.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::Indexing;
using namespace EngineLayer::ModernSearch;
using namespace EngineLayer::NonSpecificEnzymeSearch;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics;
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

#ifdef LATER
    std::cout << ++i << ". TestClassicSearchEngine " << std::endl;
    Test::SearchEngineTests::TestClassicSearchEngine();

    std::cout << ++i << ". TestClassicSearchEngineWithWeirdPeptide " << std::endl;
    Test::SearchEngineTests::TestClassicSearchEngineWithWeirdPeptide();
#endif
    
    std::cout << ++i << ". TestModernSearchEngine " << std::endl;
    Test::SearchEngineTests::TestModernSearchEngine();
    
#ifdef LATER
    std::cout << ++i << ". TestModernSearchFragmentEdges " << std::endl;
    Test::SearchEngineTests::TestModernSearchFragmentEdges();

    std::cout << ++i << ". TestNonViablePSM " << std::endl;
    Test::SearchEngineTests::TestNonViablePSM();
    
    std::cout << ++i << ". TestModernSearchEngineWithWeirdPeptide " << std::endl;
    Test::SearchEngineTests::TestModernSearchEngineWithWeirdPeptide();
    
    std::cout << ++i << ". TestNonSpecificEnzymeSearchEngineSingleN " << std::endl;
    Test::SearchEngineTests::TestNonSpecificEnzymeSearchEngineSingleN();
    
    std::cout << ++i << ". TestNonSpecificEnzymeSearchEngineSingleCModifications " << std::endl;
    Test::SearchEngineTests::TestNonSpecificEnzymeSearchEngineSingleCModifications();

    std::cout << ++i << ". TestNonSpecificEnzymeSearchEngineSingleNModifications " << std::endl;
    Test::SearchEngineTests::TestNonSpecificEnzymeSearchEngineSingleNModifications();
    
    std::cout << ++i << ". TestNonSpecificEnzymeSearchEngineSingleC " << std::endl;
    Test::SearchEngineTests::TestNonSpecificEnzymeSearchEngineSingleC();
    
    std::cout << ++i << ". TestNonSpecificEnzymeVariableModificationHandlingNTerm " << std::endl;
    Test::SearchEngineTests::TestNonSpecificEnzymeVariableModificationHandlingNTerm();

    std::cout << ++i << ". TestNonSpecificEnzymeVariableModificationHandlingCTerm " << std::endl;
    Test::SearchEngineTests::TestNonSpecificEnzymeVariableModificationHandlingCTerm();

    std::cout << ++i << ". TestSemiSpecificEnzymeEngineSingleN " << std::endl;
    Test::SearchEngineTests::TestSemiSpecificEnzymeEngineSingleN();
    
    std::cout << ++i << ". TestSemiSpecificEnzymeEngineSingleC " << std::endl;
    Test::SearchEngineTests::TestSemiSpecificEnzymeEngineSingleC();

    std::cout << ++i << ". TestClassicSemiProtease " << std::endl;
    Test::SearchEngineTests::TestClassicSemiProtease();
    
    std::cout << ++i << ". TestClassicSemiProteolysis " << std::endl;
    Test::SearchEngineTests::TestClassicSemiProteolysis();
    
    std::cout << ++i << ".  TestClassicSearchOneNterminalModifiedPeptideOneScan " << std::endl;
    Test::SearchEngineTests::TestClassicSearchOneNterminalModifiedPeptideOneScan();
#endif
    return 0;
}

namespace Test
{
#ifdef LATER
    void SearchEngineTests::TestClassicSearchEngine()
    {
        Protease *protease = new Protease("Customized Protease", CleavageSpecificity::Full, nullptr, nullptr,
                                          std::vector<DigestionMotif*> {new DigestionMotif("K", nullptr, 1, "")});
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        DigestionParams tempVar(protease: protease->Name, minPeptideLength: 1);
        CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1);
        
        auto myMsDataFile = new TestDataFile();
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto proteinList = std::vector<Protein*> {new Protein("MNNNKQQQ", nullptr)};
        
        auto searchModes = new SinglePpmAroundZeroSearchMode(5);
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
        
        std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
        ClassicSearchEngine tempVar3(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications,
                                     proteinList, searchModes,
                                     CommonParameters, new std::vector<std::string>());
        (&tempVar3)->Run();
        
        // Single search mode
        Assert::AreEqual(1, allPsmsArray.size());
        
        // One scan
        Assert::AreEqual(1, allPsmsArray.size());
        
        Assert::IsTrue(allPsmsArray[0]->getScore() > 1);
        Assert::AreEqual(2, allPsmsArray[0]->getScanNumber());
        
        Assert::AreEqual("QQQ", allPsmsArray[0]->getBaseSequence());
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
    }
    
    void SearchEngineTests::TestClassicSearchEngineWithWeirdPeptide()
    {
        DigestionParams tempVar(protease: "Customized Protease", maxMissedCleavages: 0, minPeptideLength: 1);
        CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1);
        
        auto myMsDataFile = new TestDataFile();
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto proteinList = std::vector<Protein*> {new Protein("QXQ", nullptr)};
        
        auto searchModes = new OpenSearchMode();
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
        
        std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
        ClassicSearchEngine tempVar3(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes,
                                     CommonParameters, new std::vector<std::string>());
        (&tempVar3)->Run();
        
        // Single search mode
        Assert::AreEqual(1, allPsmsArray.size());
        
        // One Scan
        Assert::AreEqual(1, allPsmsArray.size());
        
        Assert::IsTrue(allPsmsArray[0]->getScore() > 1);
        Assert::AreEqual(2, allPsmsArray[0]->getScanNumber());
        
        Assert::AreEqual("QXQ", allPsmsArray[0]->getBaseSequence());
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
    }
#endif    

    void SearchEngineTests::TestModernSearchEngine()
    {
        SearchParameters *SearchParameters = new SearchParameters();
        SearchParameters->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
        SearchParameters->setSearchTarget(true);
        
        PpmTolerance tempVar(5);
        DigestionParams tempVar2(protease: "Customized Protease", minPeptideLength: 1);
        CommonParameters *CommonParameters = new CommonParameters(precursorMassTolerance: &tempVar, digestionParams: &tempVar2, scoreCutoff: 1);
        
        auto myMsDataFile = new TestDataFile();
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        
        auto proteinList = std::vector<Protein*> {new Protein("MNNNKQQQ", nullptr)};
        
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters,
                                              SearchParameters->getMaxFragmentSize(), false, std::vector<FileInfo*>(), std::vector<std::string>());
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        
        CommonParameters tempVar3();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar3).OrderBy([&] (std::any b) {
                b::PrecursorMass;
            })->ToArray();
        
        MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(),
                                                                             SearchParameters->getMassDiffAcceptorType(),
                                                                             SearchParameters->getCustomMdac());
        
        std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
        ModernSearchEngine tempVar4(allPsmsArray, listOfSortedms2Scans, indexResults->getPeptideIndex(),
                                    indexResults->getFragmentIndex(), 0, CommonParameters,
                                    massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(),
                                    new std::vector<std::string>());
        (&tempVar4)->Run();
        
        // Single search mode
        Assert::AreEqual(1, (int)allPsmsArray.size());
        
        // Single ms2 scan
        Assert::AreEqual(1, (int)allPsmsArray.size());
        Assert::IsTrue(allPsmsArray[0] != nullptr);
        
        Assert::IsTrue(allPsmsArray[0]->getScore() > 1);
        Assert::AreEqual(2, allPsmsArray[0]->getScanNumber());
        
        Assert::AreEqual("QQQ", allPsmsArray[0]->getBaseSequence());
        
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
        delete SearchParameters;
    }

#ifdef LATER    
    void SearchEngineTests::TestModernSearchFragmentEdges()
    {
        // purpose of this test is to check that fragments at the edge of the fragment index in modern search do not fall out of bounds
        // example: fragment mass 400 when max frag mass is 400 may go over the edge of the array because of ppm tolerances
        SearchParameters *SearchParameters = new SearchParameters();
        SearchParameters->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
        SearchParameters->setMaxFragmentSize(1);
        AbsoluteTolerance tempVar(100);
        DigestionParams tempVar2(protease: "trypsin", minPeptideLength: 1);
        CommonParameters *CommonParameters = new CommonParameters(productMassTolerance: &tempVar, digestionParams: &tempVar2, scoreCutoff: 1, addCompIons: true);
        
        auto proteinList = std::vector<Protein*> {new Protein("K", nullptr)};
        
        auto indexEngine = new IndexingEngine(proteinList, std::vector<Modification*>(), std::vector<Modification*>(), 1, DecoyType::Reverse, CommonParameters,
                                              SearchParameters->getMaxFragmentSize(), false, std::vector<FileInfo*>(), std::vector<std::string>());
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        
        TestDataFile tempVar3();
        CommonParameters tempVar4();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(&tempVar3, "", &tempVar4).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
        
        MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(), SearchParameters->getMassDiffAcceptorType(),
                                                                             SearchParameters->getCustomMdac());
        
        std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
        ModernSearchEngine tempVar5(allPsmsArray, listOfSortedms2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(), 0, CommonParameters,
                                    massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(), new std::vector<std::string>());
        (&tempVar5)->Run();
        
        // no assertions... just don't crash...
        
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
        delete SearchParameters;
    }
    
    void SearchEngineTests::TestNonViablePSM()
    {
        //just check if it crashes or not.
        SearchParameters *SearchParameters = new SearchParameters();
        
        DigestionParams tempVar(protease: "Customized Protease");
        CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 255);
        
        MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(), SearchParameters->getMassDiffAcceptorType(),
                                                                             SearchParameters->getCustomMdac());
        
        auto myMsDataFile = new TestDataFile(true); //empty
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        
        auto proteinList = std::vector<Protein*> {new Protein("MNNNKQQQ", nullptr)};
        
        auto searchModes = new SinglePpmAroundZeroSearchMode(5);
        
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters,
                                              SearchParameters->getMaxFragmentSize(), false, std::vector<FileInfo*>(), std::vector<std::string>());
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b) {
                b::PrecursorMass;
            })->ToArray();
        
        std::vector<std::vector<PeptideSpectralMatch*>> allPsmsArrays(1);
        allPsmsArrays[0] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        std::vector<PeptideSpectralMatch*> allPsmsArray = allPsmsArrays[0];
        
        //Classic
        ClassicSearchEngine tempVar3(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes,
                                     CommonParameters, new std::vector<std::string>());
        (&tempVar3)->Run();
        
        //Modern
        ModernSearchEngine tempVar4(allPsmsArray, listOfSortedms2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(), 0, CommonParameters,
                                    massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(), new std::vector<std::string>());
        (&tempVar4)->Run();
        
        //NonSpecific
        NonSpecificEnzymeSearchEngine tempVar5(allPsmsArrays, listOfSortedms2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(),
                                               indexResults->getFragmentIndex(), 0, CommonParameters, massDiffAcceptor,
                                               SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(), new std::vector<std::string>());
        (&tempVar5)->Run();
        
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
        delete SearchParameters;
    }
    
    void SearchEngineTests::TestModernSearchEngineWithWeirdPeptide()
    {
        SearchParameters *SearchParameters = new SearchParameters();
        SearchParameters->setMassDiffAcceptorType(MassDiffAcceptorType::Open);
        SearchParameters->setSearchTarget(true);
        
        DigestionParams tempVar(protease: "Customized Protease", minPeptideLength: 1);
        CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1);
        
        auto myMsDataFile = new TestDataFile();
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto localizeableModifications = std::vector<Modification*>();
        
        auto proteinList = std::vector<Protein*> {new Protein("MNNNKQXQ", nullptr)};
        
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters,
                                              SearchParameters->getMaxFragmentSize(),
                                              false, std::vector<FileInfo*>(), std::vector<std::string>());
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        
        bool DoPrecursorDeconvolution = true;
        bool UseProvidedPrecursorInfo = true;
        double DeconvolutionIntensityRatio = 4;
        int DeconvolutionMaxAssumedChargeState = 10;
        Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
        MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(),
                                                                             SearchParameters->getMassDiffAcceptorType(), SearchParameters->getCustomMdac());
        
        std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
        auto engine = new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults->getPeptideIndex(), indexResults->getFragmentIndex(), 0,
                                             CommonParameters, massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(),
                                             std::vector<std::string>());
        auto searchResults = engine->Run();
        
        // Single search mode
        Assert::AreEqual(1, allPsmsArray.size());
        
        // Single ms2 scan
        Assert::AreEqual(1, allPsmsArray.size());
        Assert::That(allPsmsArray[0] != nullptr);
        
        Assert::IsTrue(allPsmsArray[0]->getScore() > 1);
        Assert::AreEqual(2, allPsmsArray[0]->getScanNumber());
        
        Assert::AreEqual(3, allPsmsArray[0]->getNumDifferentMatchingPeptides());
        
        delete engine;
        delete DeconvolutionMassTolerance;
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
        delete SearchParameters;
    }
    
    void SearchEngineTests::TestNonSpecificEnzymeSearchEngineSingleN()
    {
        SearchParameters *SearchParameters = new SearchParameters();
        SearchParameters->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
        SearchParameters->setLocalFdrCategories(std::vector<FdrCategory> {FdrCategory::FullySpecific, FdrCategory::SemiSpecific, FdrCategory::NonSpecific});
        std::vector<DigestionMotif*> motifs =
            {
                new DigestionMotif("K", nullptr, 0, nullptr),
                new DigestionMotif("K", nullptr, 1, nullptr),
                new DigestionMotif("G", nullptr, 0, nullptr),
                new DigestionMotif("G", nullptr, 1, nullptr)
            };
        Protease *protease = new Protease("single N", CleavageSpecificity::SingleN, nullptr, nullptr, motifs);
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        DigestionParams *dp = new DigestionParams(protease: protease->Name, minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus::N,
                                                  searchModeType: CleavageSpecificity::None);
        PpmTolerance tempVar(5);
        CommonParameters *CommonParameters = new CommonParameters(dissociationType: DissociationType::HCD, precursorMassTolerance: &tempVar,
                                                                  digestionParams: dp, scoreCutoff: 1, addCompIons: true);
        
        auto myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        
        auto proteinList = std::vector<Protein*> {new Protein("GGGGGMNNNKQQQGGGGG", "TestProtein")};
        
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters, 100000,
                                              false, std::vector<FileInfo*>(), std::vector<std::string>());
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        auto peptideIndex = indexResults->getPeptideIndex();
        auto fragmentIndexDict = indexResults->getFragmentIndex();
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
        
        MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(), SearchParameters->getMassDiffAcceptorType(),
                                                                             SearchParameters->getCustomMdac());
        
        std::vector<std::vector<PeptideSpectralMatch*>> allPsmsArrays(3);
        allPsmsArrays[0] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[1] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[2] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        std::vector<PeptideSpectralMatch*> allPsmsArray = allPsmsArrays[2];
        NonSpecificEnzymeSearchEngine tempVar3(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, 0, CommonParameters,
                                               massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(), new std::vector<std::string>());
        (&tempVar3)->Run();
        
        // Single search mode
        Assert::AreEqual(1, allPsmsArray.size());
        
        // Single ms2 scan
        Assert::AreEqual(1, allPsmsArray.size());
        
        Assert::IsTrue(allPsmsArray[0]->getScore() > 4);
        Assert::AreEqual(2, allPsmsArray[0]->getScanNumber());
        DigestionParams tempVar4(protease: protease->Name, minPeptideLength: 1);
        PpmTolerance tempVar5(5);
        CommonParameters = new CommonParameters(digestionParams: &tempVar4, precursorMassTolerance: &tempVar5, scoreCutoff: 1);
        
        std::unordered_map<CompactPeptideBase*, std::unordered_set<PeptideWithSetModifications*>> compactPeptideToProteinPeptideMatching;
        allPsmsArray[0]->ResolveAllAmbiguities();
        Assert::AreEqual("QQQGGGG", allPsmsArray[0]->getBaseSequence());
        
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
        delete dp;
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease
        //was passed to a method or constructor. Handle memory management manually.
        delete SearchParameters;
    }
    
    void SearchEngineTests::TestNonSpecificEnzymeSearchEngineSingleCModifications()
    {
        SearchParameters *SearchParameters = new SearchParameters();
        SearchParameters->setSearchTarget(true);
        SearchParameters->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
        SearchParameters->setLocalFdrCategories(std::vector<FdrCategory> {FdrCategory::NonSpecific});
        DigestionParams *dp = new DigestionParams("singleC", minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus::C,
                                                  searchModeType: CleavageSpecificity::None);
        
        PpmTolerance tempVar(5);
        CommonParameters *CommonParameters = new CommonParameters(dissociationType: DissociationType::HCD, digestionParams: dp, scoreCutoff: 5,
                                                                  precursorMassTolerance: &tempVar, addCompIons: true);
        
        PeptideWithSetModifications *guiltyPwsm = new PeptideWithSetModifications("DQPKLLGIETPLPKKE", nullptr);
        auto fragments = guiltyPwsm->Fragment(CommonParameters->getDissociationType(), FragmentationTerminus::Both);
        
        
        auto myMsDataFile = new TestDataFile(guiltyPwsm->MonoisotopicMass, fragments->Select([&] (std::any x)  {
                    x::NeutralMass::ToMz(1);
                })->ToArray());
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        ModificationMotif motif2;
        ModificationMotif::TryGetMotif("C", motif2);
        Modification *mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif2,
                                              _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
        fixedModifications.push_back(mod2);
        
        auto proteinList = std::vector<Protein*> {new Protein("GGGGGCDQPKLLGIETPLPKKEGGGGG", nullptr)};
        
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::None,
                                              CommonParameters, SearchParameters->getMaxFragmentSize(),true, std::vector<FileInfo*>(),
                                              std::vector<std::string>());
        
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        auto peptideIndex = indexResults->getPeptideIndex();
        auto fragmentIndexDict = indexResults->getFragmentIndex();
        auto precursorIndexDict = indexResults->getPrecursorIndex();
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b)    {
                b::PrecursorMass;
            })->ToArray();
        MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(),
                                                                             SearchParameters->getMassDiffAcceptorType(),
                                                                             SearchParameters->getCustomMdac());
        
        std::vector<std::vector<PeptideSpectralMatch*>> allPsmsArrays(3);
        allPsmsArrays[0] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[1] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[2] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        std::vector<PeptideSpectralMatch*> allPsmsArray = allPsmsArrays[2];
        auto engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters,
                                                        massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(), std::vector<std::string>());
        auto searchResults = engine->Run();
        
        allPsmsArray[0]->ResolveAllAmbiguities();
        //Check that there is no modification hanging out on the n-terminus
        Assert::AreEqual(allPsmsArray[0]->getFullSequence(),guiltyPwsm->FullSequence);
        
        
        proteinList = {new Protein("CDQPKLLGIETPLPKKEGGGGG", nullptr)};
        guiltyPwsm = new PeptideWithSetModifications("C[Common Fixed:Carbamidomethyl on C]DQPKLLGIETPLPKKE",
                                                     std::unordered_map<std::string, Modification*> {{"Carbamidomethyl on C", mod2}});
        fragments = guiltyPwsm->Fragment(CommonParameters->getDissociationType(), FragmentationTerminus::Both);
        myMsDataFile = new TestDataFile(guiltyPwsm->MonoisotopicMass, fragments->Select([&] (std::any x)  {
                    x::NeutralMass::ToMz(1);
                })->ToArray());
        indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::None, CommonParameters,
                                         SearchParameters->getMaxFragmentSize(), true, std::vector<FileInfo*>(), std::vector<std::string>());
        indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        precursorIndexDict = indexResults->getPrecursorIndex();
        peptideIndex = indexResults->getPeptideIndex();
        fragmentIndexDict = indexResults->getFragmentIndex();
        CommonParameters tempVar3();
        listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar3).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
        
        allPsmsArrays[2] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArray = allPsmsArrays[2];
        engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters,
                                                   massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(), std::vector<std::string>());
        searchResults = engine->Run();
        allPsmsArray[0]->ResolveAllAmbiguities();
        //Check that there is a modification hanging out on the protein n-terminus
        Assert::AreEqual(allPsmsArray[0]->getFullSequence(), guiltyPwsm->FullSequence);
        
        proteinList = {new Protein("GGGGGCDQPKLLGIETPLPKKEGG", nullptr)};
        indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::None,
                                         CommonParameters, SearchParameters->getMaxFragmentSize(),
                                         true, std::vector<FileInfo*>(), std::vector<std::string>());
        indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        peptideIndex = indexResults->getPeptideIndex();
        fragmentIndexDict = indexResults->getFragmentIndex();
        precursorIndexDict = indexResults->getPrecursorIndex();
        allPsmsArrays[2] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArray = allPsmsArrays[2];
        engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters,
                                                   massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(), std::vector<std::string>());
        searchResults = engine->Run();
        allPsmsArray[0]->ResolveAllAmbiguities();
        //Check that there is a modification hanging out on the peptide n-terminus
        Assert::AreEqual(allPsmsArray[0]->getFullSequence(), guiltyPwsm->FullSequence);
        
        delete engine;
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete mod2' statement was not added since mod2
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
        delete guiltyPwsm;
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
        delete dp;
        delete SearchParameters;
    }
    
    void SearchEngineTests::TestNonSpecificEnzymeSearchEngineSingleNModifications()
    {
        SearchParameters *SearchParameters = new SearchParameters();
        SearchParameters->setSearchTarget(true);
        SearchParameters->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
        SearchParameters->setLocalFdrCategories(std::vector<FdrCategory> {FdrCategory::NonSpecific});
        DigestionParams *dp = new DigestionParams("singleN", minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus::N,
                                                  searchModeType: CleavageSpecificity::None);
        
        PpmTolerance tempVar(5);
        CommonParameters *CommonParameters = new CommonParameters(dissociationType: DissociationType::HCD, digestionParams: dp,
                                                                  scoreCutoff: 5, precursorMassTolerance: &tempVar,
                                                                  addCompIons: true);
        
        PeptideWithSetModifications *guiltyPwsm = new PeptideWithSetModifications("DQPKLLGIETPLPKKE", nullptr);
        auto fragments = guiltyPwsm->Fragment(CommonParameters->getDissociationType(), FragmentationTerminus::Both);
        
        
        auto myMsDataFile = new TestDataFile(guiltyPwsm->MonoisotopicMass, fragments->Select([&] (std::any x) {
                    x::NeutralMass::ToMz(1);
                })->ToArray());
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        ModificationMotif motif2;
        ModificationMotif::TryGetMotif("C", motif2);
        Modification *mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed",
                                              _target: motif2, _locationRestriction: "Anywhere.",
                                              _monoisotopicMass: 57.02146372068994);
        fixedModifications.push_back(mod2);
        
        auto proteinList = std::vector<Protein*> {new Protein("GGDQPKLLGIETPLPKKECGGGGG", nullptr)};
        
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::None,
                                              CommonParameters, SearchParameters->getMaxFragmentSize(),
                                              true, std::vector<FileInfo*>(), std::vector<std::string>());
        
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        auto peptideIndex = indexResults->getPeptideIndex();
        auto fragmentIndexDict = indexResults->getFragmentIndex();
        auto precursorIndexDict = indexResults->getPrecursorIndex();
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b) {
                b::PrecursorMass;
            })->ToArray();
        MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(),
                                                                             SearchParameters->getMassDiffAcceptorType(),
                                                                             SearchParameters->getCustomMdac());
        
        std::vector<std::vector<PeptideSpectralMatch*>> allPsmsArrays(3);
        allPsmsArrays[0] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[1] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[2] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        std::vector<PeptideSpectralMatch*> allPsmsArray = allPsmsArrays[2];
        auto engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0,
                                                        CommonParameters, massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(),
                                                        std::vector<std::string>());
        auto searchResults = engine->Run();
        
        allPsmsArray[0]->ResolveAllAmbiguities();
        //Check that there is no modification hanging out on the n-terminus
        Assert::AreEqual(allPsmsArray[0]->getFullSequence(), guiltyPwsm->FullSequence);
        
        
        proteinList = {new Protein("GGGGGDQPKLLGIETPLPKKEC", nullptr)};
        guiltyPwsm = new PeptideWithSetModifications("GGDQPKLLGIETPLPKKEC[Common Fixed:Carbamidomethyl on C]", std::unordered_map<std::string, Modification*>   {
                {"Carbamidomethyl on C", mod2}});
        fragments = guiltyPwsm->Fragment(CommonParameters->getDissociationType(), FragmentationTerminus::Both);
        myMsDataFile = new TestDataFile(guiltyPwsm->MonoisotopicMass, fragments->Select([&] (std::any x)
        {
            x::NeutralMass::ToMz(1);
        })->ToArray());
        indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::None, CommonParameters,
                                         SearchParameters->getMaxFragmentSize(), true, std::vector<FileInfo*>(),
                                         std::vector<std::string>());
        indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        precursorIndexDict = indexResults->getPrecursorIndex();
        peptideIndex = indexResults->getPeptideIndex();
        fragmentIndexDict = indexResults->getFragmentIndex();
        CommonParameters tempVar3();
        listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar3).OrderBy([&] (std::any b)
                                                                                                  {
                                                                                                      b::PrecursorMass;
                                                                                                  })->ToArray();
        
        allPsmsArrays[2] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArray = allPsmsArrays[2];
        engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0,
                                                   CommonParameters, massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(),
                                                   std::vector<std::string>());
        searchResults = engine->Run();
        allPsmsArray[0]->ResolveAllAmbiguities();
        //Check that there is a modification hanging out on the protein n-terminus
        Assert::AreEqual(allPsmsArray[0]->getFullSequence(), guiltyPwsm->FullSequence);
        
        proteinList = {new Protein("GGDQPKLLGIETPLPKKECGGGGG", nullptr)};
        indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::None, CommonParameters,
                                         SearchParameters->getMaxFragmentSize(), true, std::vector<FileInfo*>(), std::vector<std::string>());
        indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        peptideIndex = indexResults->getPeptideIndex();
        fragmentIndexDict = indexResults->getFragmentIndex();
        precursorIndexDict = indexResults->getPrecursorIndex();
        allPsmsArrays[2] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArray = allPsmsArrays[2];
        engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0,
                                                   CommonParameters, massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(),
                                                   std::vector<std::string>());
        searchResults = engine->Run();
        allPsmsArray[0]->ResolveAllAmbiguities();
        //Check that there is a modification hanging out on the peptide n-terminus
        Assert::AreEqual(allPsmsArray[0]->getFullSequence(), guiltyPwsm->FullSequence);
        
        delete engine;
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete mod2' statement was not added since mod2
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
        delete guiltyPwsm;
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
        delete dp;
        delete SearchParameters;
    }

    void SearchEngineTests::TestNonSpecificEnzymeSearchEngineSingleC()
    {
        SearchParameters *SearchParameters = new SearchParameters();
        SearchParameters->setSearchTarget(true);
        SearchParameters->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
        SearchParameters->setLocalFdrCategories(std::vector<FdrCategory> {FdrCategory::FullySpecific, FdrCategory::SemiSpecific, FdrCategory::NonSpecific});
        
        std::vector<DigestionMotif*> motifs =
            {
                new DigestionMotif("K", nullptr, 0, nullptr),
                new DigestionMotif("K", nullptr, 1, nullptr),
                new DigestionMotif("G", nullptr, 0, nullptr),
                new DigestionMotif("G", nullptr, 1, nullptr)
            };
        Protease *protease = new Protease("single C", CleavageSpecificity::SingleC, nullptr, nullptr, motifs);
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        DigestionParams *dp = new DigestionParams(protease: protease->Name, minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus::C,
                                                  searchModeType: CleavageSpecificity::None);
        
        PpmTolerance tempVar(5);
        CommonParameters *CommonParameters = new CommonParameters(dissociationType: DissociationType::HCD, digestionParams: dp, scoreCutoff: 4,
                                                                  precursorMassTolerance: &tempVar, addCompIons: true);
        
        auto myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto localizeableModifications = std::vector<Modification*>();
        
        auto proteinList = std::vector<Protein*> {new Protein("GGGGGMNNNKQQQGGGGG", nullptr)};
        
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters,
                                              SearchParameters->getMaxFragmentSize(), false, std::vector<FileInfo*>(), std::vector<std::string>());
        
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        auto peptideIndex = indexResults->getPeptideIndex();
        auto fragmentIndexDict = indexResults->getFragmentIndex();
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
        MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(CommonParameters->getPrecursorMassTolerance(),
                                                                             SearchParameters->getMassDiffAcceptorType(),
                                                                             SearchParameters->getCustomMdac());
        
        std::vector<std::vector<PeptideSpectralMatch*>> allPsmsArrays(3);
        allPsmsArrays[0] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[1] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[2] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        std::vector<PeptideSpectralMatch*> allPsmsArray = allPsmsArrays[2];
        auto engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, 0,
                                                        CommonParameters, massDiffAcceptor, SearchParameters->getMaximumMassThatFragmentIonScoreIsDoubled(),
                                                        std::vector<std::string>());
        auto searchResults = engine->Run();
        
        // Single search mode
        Assert::AreEqual(1, allPsmsArray.size());
        
        //Single ms2 scan
        Assert::AreEqual(1, allPsmsArray.size());
        Assert::That(allPsmsArray[0] != nullptr);
        
        Assert::IsTrue(allPsmsArray[0]->getScore() > 7);
        Assert::AreEqual(2, allPsmsArray[0]->getScanNumber());
        allPsmsArray[0]->ResolveAllAmbiguities();
        Assert::AreEqual("QQQGGGG", allPsmsArray[0]->getBaseSequence());
        
        delete engine;
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
        delete dp;
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease
        //was passed to a method or constructor. Handle memory management manually.
        delete SearchParameters;
    }
    
    void SearchEngineTests::TestNonSpecificEnzymeVariableModificationHandlingNTerm()
    {
        auto protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
        
        std::vector<DigestionMotif*> motifs =
            {
                new DigestionMotif("K", nullptr, 0, nullptr),
                new DigestionMotif("K", nullptr, 1, nullptr),
                new DigestionMotif("G", nullptr, 0, nullptr),
                new DigestionMotif("G", nullptr, 1, nullptr),
                new DigestionMotif("M", nullptr, 0, nullptr),
                new DigestionMotif("M", nullptr, 1, nullptr),
                new DigestionMotif("N", nullptr, 0, nullptr),
                new DigestionMotif("N", nullptr, 1, nullptr),
                new DigestionMotif("Q", nullptr, 0, nullptr),
                new DigestionMotif("Q", nullptr, 1, nullptr)
            };
        auto protease = new Protease("singleN2", CleavageSpecificity::SingleN, nullptr, nullptr, motifs);
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        ModificationMotif motifM;
        ModificationMotif::TryGetMotif("M", motifM);
        auto variableModifications = std::vector<Modification*> {new Modification(_originalId: "16", _target: motifM, _locationRestriction: "Anywhere.",
                                                                                  _monoisotopicMass: 15.994915)};
        DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 5, maxModsForPeptides: 3);
        auto ListOfModifiedPeptides = protein->Digest(digestionParams, std::vector<Modification*>(), variableModifications).ToList();
        Assert::AreEqual(ListOfModifiedPeptides.size(), 192);
        
        auto protein2 = new Protein(std::string((std::string("MGGGGGMNNNKQQQMGGGGMGM")).ToCharArray().Reverse().ToArray()), "TestProtein");
        auto ListOfModifiedPeptides2 = protein2->Digest(digestionParams, std::vector<Modification*>(), variableModifications).ToList();
        Assert::AreEqual(ListOfModifiedPeptides2.size(), 132);
        
        delete protein2;
        //C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease
        //was passed to a method or constructor. Handle memory management manually.
        delete protein;
    }
    
    void SearchEngineTests::TestNonSpecificEnzymeVariableModificationHandlingCTerm()
    {
        auto protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
        std::vector<DigestionMotif*> motifs =
            {
                new DigestionMotif("K", nullptr, 0, nullptr),
                new DigestionMotif("K", nullptr, 1, nullptr),
                new DigestionMotif("G", nullptr, 0, nullptr),
                new DigestionMotif("G", nullptr, 1, nullptr),
                new DigestionMotif("M", nullptr, 0, nullptr),
                new DigestionMotif("M", nullptr, 1, nullptr),
                new DigestionMotif("N", nullptr, 0, nullptr),
                new DigestionMotif("N", nullptr, 1, nullptr),
                new DigestionMotif("Q", nullptr, 0, nullptr),
                new DigestionMotif("Q", nullptr, 1, nullptr)
            };
        auto protease = new Protease("singleC2", CleavageSpecificity::SingleC, nullptr, nullptr, motifs);
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        ModificationMotif motifM;
        ModificationMotif::TryGetMotif("M", motifM);
        auto variableModifications = std::vector<Modification*> {new Modification(_originalId: "16", _target: motifM, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915)};
        DigestionParams *digestionParams = new DigestionParams(protease: protease->Name, minPeptideLength: 5, maxModsForPeptides: 3);
        auto ListOfModifiedPeptides = protein->Digest(digestionParams, std::vector<Modification*>(), variableModifications).ToList();
        Assert::AreEqual(ListOfModifiedPeptides.size(), 132);
        
        auto protein2 = new Protein(std::string((std::string("MGGGGGMNNNKQQQMGGGGMGM")).ToCharArray().Reverse().ToArray()), "TestProtein");
        auto ListOfModifiedPeptides2 = protein2->Digest(digestionParams, std::vector<Modification*>(), variableModifications).ToList();
        Assert::AreEqual(ListOfModifiedPeptides2.size(), 192);
        
        delete protein2;
        //C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease
        //was passed to a method or constructor. Handle memory management manually.
        delete protein;
    }
    
    void SearchEngineTests::TestSemiSpecificEnzymeEngineSingleN()
    {
        auto myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto localizeableModifications = std::vector<Modification*>();
        std::unordered_map<Modification*, unsigned short> modsDictionary;
        for (auto mod : fixedModifications)
        {
            modsDictionary.emplace(mod, 0);
        }
        int ii = 1;
        for (auto mod : variableModifications)
        {
            modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
            ii++;
        }
        for (auto mod : localizeableModifications)
        {
            modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
            ii++;
        }
        
        auto proteinList = std::vector<Protein*> {new Protein("GGGGGMNNNKQQQGGGGGGKKRKG", "TestProtein")};
        
        auto productMassTolerance = new AbsoluteTolerance(0.01);
        auto searchModes = new SinglePpmAroundZeroSearchMode(5);
        std::vector<DigestionMotif*> motifs = {new DigestionMotif("K", nullptr, 1, nullptr)};
        auto protease = new Protease("SingleN", CleavageSpecificity::None, nullptr, nullptr, motifs);
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        DigestionParams tempVar(protease: protease->Name, minPeptideLength: 5, maxModsForPeptides: 2, fragmentationTerminus: FragmentationTerminus::N,
                                searchModeType: CleavageSpecificity::Semi);
        CommonParameters *CommonParameters = new CommonParameters(dissociationType: DissociationType::HCD, productMassTolerance: productMassTolerance,
                                                                  digestionParams: &tempVar, scoreCutoff: 2, addCompIons: true);
        
        std::unordered_set<DigestionParams*> digestParams = {CommonParameters->getDigestionParams()};
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters, 100000,
                                              true, std::vector<FileInfo*>(), std::vector<std::string>());
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        auto peptideIndex = indexResults->getPeptideIndex();
        auto fragmentIndexDict = indexResults->getFragmentIndex();
        auto precursorIndexDict = indexResults->getPrecursorIndex();
        
        bool DoPrecursorDeconvolution = true;
        bool UseProvidedPrecursorInfo = true;
        double DeconvolutionIntensityRatio = 4;
        int DeconvolutionMaxAssumedChargeState = 10;
        Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b) {
                b::PrecursorMass;
            })->ToArray();
        
        std::vector<std::vector<PeptideSpectralMatch*>> allPsmsArrays(2);
        allPsmsArrays[0] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[1] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        std::vector<PeptideSpectralMatch*> allPsmsArray = allPsmsArrays[1];
        auto engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 1,
                                                        CommonParameters, searchModes, 0, std::vector<std::string>());
        auto searchResults = engine->Run();
        
        // Single search mode
        Assert::AreEqual(1, allPsmsArray.size());
        
        // Single ms2 scan
        Assert::AreEqual(1, allPsmsArray.size());
        
        Assert::IsTrue(allPsmsArray[0]->getScore() > 4);
        Assert::AreEqual(2, allPsmsArray[0]->getScanNumber());
        allPsmsArray[0]->ResolveAllAmbiguities();
        Assert::AreEqual("QQQGGGG", allPsmsArray[0]->getBaseSequence());
        
        delete engine;
        delete DeconvolutionMassTolerance;
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete productMassTolerance' statement was not added since
        //productMassTolerance was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
    }
    
    void SearchEngineTests::TestSemiSpecificEnzymeEngineSingleC()
    {
        auto myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto localizeableModifications = std::vector<Modification*>();
        std::unordered_map<Modification*, unsigned short> modsDictionary;
        for (auto mod : fixedModifications)
        {
            modsDictionary.emplace(mod, 0);
        }
        int ii = 1;
        for (auto mod : variableModifications)
        {
            modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
            ii++;
        }
        for (auto mod : localizeableModifications)
        {
            modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
            ii++;
        }
        
        auto proteinList = std::vector<Protein*> {new Protein("GGGGGMKNNNQQQGGGGKGG", nullptr, nullptr, nullptr, nullptr,
                                                              std::vector<ProteolysisProduct*> {new ProteolysisProduct(nullptr, nullptr, "test")})};
        
        auto productMassTolerance = new AbsoluteTolerance(0.01);
        auto searchModes = new SinglePpmAroundZeroSearchMode(5);
        Protease *protease = new Protease("SingleC", CleavageSpecificity::Full, nullptr, nullptr, std::vector<DigestionMotif*>());
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        DigestionParams tempVar(protease: protease->Name, maxMissedCleavages: 5, minPeptideLength: 5, searchModeType: CleavageSpecificity::None,
                                fragmentationTerminus: FragmentationTerminus::C);
        CommonParameters *CommonParameters = new CommonParameters(scoreCutoff: 1, productMassTolerance: productMassTolerance,
                                                                  digestionParams: &tempVar, dissociationType: DissociationType::HCD,
                                                                  addCompIons: true);
        
        std::unordered_set<DigestionParams*> digestParams = {CommonParameters->getDigestionParams()};
        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse,
                                              CommonParameters, 30000, false, std::vector<FileInfo*>(), std::vector<std::string>());
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        auto peptideIndex = indexResults->getPeptideIndex();
        auto fragmentIndexDict = indexResults->getFragmentIndex();
        
        bool DoPrecursorDeconvolution = true;
        bool UseProvidedPrecursorInfo = true;
        double DeconvolutionIntensityRatio = 4;
        int DeconvolutionMaxAssumedChargeState = 10;
        Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
        
        std::vector<std::vector<PeptideSpectralMatch*>> allPsmsArrays(3);
        allPsmsArrays[0] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[1] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        allPsmsArrays[2] = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        std::vector<PeptideSpectralMatch*> allPsmsArray = allPsmsArrays[2];
        auto engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, 1,
                                                        CommonParameters, searchModes, 0, std::vector<std::string>());
        auto searchResults = engine->Run();
        
        // Single search mode
        Assert::AreEqual(1, allPsmsArray.size());
        
        // Single ms2 scan
        Assert::AreEqual(1, allPsmsArray.size());
        
        Assert::IsTrue(allPsmsArray[0]->getScore() > 4);
        Assert::AreEqual(2, allPsmsArray[0]->getScanNumber());
        allPsmsArray[0]->ResolveAllAmbiguities();
        Assert::AreEqual("QQQGGGG", allPsmsArray[0]->getBaseSequence());
        
        delete engine;
        delete DeconvolutionMassTolerance;
        delete indexEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete productMassTolerance' statement was not added since productMassTolerance
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
    }
    
    void SearchEngineTests::TestClassicSemiProtease()
    {
        auto myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto localizeableModifications = std::vector<Modification*>();
        std::unordered_map<Modification*, unsigned short> modsDictionary;
        for (auto mod : fixedModifications)
        {
            modsDictionary.emplace(mod, 0);
        }
        int ii = 1;
        for (auto mod : variableModifications)
        {
            modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
            ii++;
        }
        for (auto mod : localizeableModifications)
        {
            modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
            ii++;
        }
        
        auto proteinList = std::vector<Protein*> {new Protein("MGGGGGMKNNNQQQGGGGKGKKNKKGN", "hello")};
        
        auto productMassTolerance = new AbsoluteTolerance(0.01);
        auto searchModes = new SinglePpmAroundZeroSearchMode(5);
        
        std::vector<DigestionMotif*> motifs1 = {new DigestionMotif("G", nullptr, 1, nullptr)};
        std::vector<DigestionMotif*> motifs2 = {new DigestionMotif("N", nullptr, 1, nullptr)};
        
        auto protease = new Protease("semi-trypsin1", CleavageSpecificity::Semi, nullptr, nullptr, motifs1);
        auto protease2 = new Protease("semi-trypsin2", CleavageSpecificity::Semi, nullptr, nullptr, motifs2);
        
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        ProteaseDictionary::Dictionary->Add(protease2->Name, protease2);
        DigestionParams tempVar(protease: protease->Name, maxMissedCleavages: 5);
        CommonParameters *CommonParameters = new CommonParameters(productMassTolerance: productMassTolerance, digestionParams: &tempVar);
        
        bool DoPrecursorDeconvolution = true;
        bool UseProvidedPrecursorInfo = true;
        double DeconvolutionIntensityRatio = 4;
        int DeconvolutionMaxAssumedChargeState = 10;
        Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);
        
        CommonParameters tempVar2();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar2).OrderBy([&] (std::any b)
                                                                                                       {
                                                                                                           b::PrecursorMass;
                                                                                                       })->ToArray();
        
        std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
        ClassicSearchEngine tempVar3(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes,
                                     CommonParameters, new std::vector<std::string>());
        (&tempVar3)->Run();
        
        //////////////////////////////
        
        DigestionParams tempVar4(protease: protease2->Name, maxMissedCleavages: 5);
        CommonParameters *CommonParameters2 = new CommonParameters(productMassTolerance: productMassTolerance, digestionParams: &tempVar4);
        
        std::unordered_set<DigestionParams*> digestParams2 = {CommonParameters2->getDigestionParams()};
        
        bool DoPrecursorDeconvolution2 = true;
        bool UseProvidedPrecursorInfo2 = true;
        double DeconvolutionIntensityRatio2 = 4;
        int DeconvolutionMaxAssumedChargeState2 = 10;
        Tolerance *DeconvolutionMassTolerance2 = new PpmTolerance(5);
        
        CommonParameters tempVar5();
        auto listOfSortedms2Scans2 = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", &tempVar5).OrderBy([&] (std::any b)   {
                b::PrecursorMass;
            })->ToArray();
        
        std::vector<PeptideSpectralMatch*> allPsmsArray2(listOfSortedms2Scans.size());
        ClassicSearchEngine tempVar6(allPsmsArray2, listOfSortedms2Scans2, variableModifications, fixedModifications, proteinList, searchModes,
                                     CommonParameters2, new std::vector<std::string>());
        (&tempVar6)->Run();
        // Single search mode
        Assert::AreEqual(1, allPsmsArray2.size());
        Assert::AreEqual(allPsmsArray.size(), allPsmsArray2.size());
        
        // Single ms2 scan
        Assert::AreEqual(1, allPsmsArray2.size());
        Assert::AreEqual(allPsmsArray.size(), allPsmsArray2.size());
        
        Assert::IsTrue(allPsmsArray2[0]->getScore() > 4);
        Assert::IsTrue(allPsmsArray[0]->getScore() > 4);
        Assert::AreEqual(2, allPsmsArray2[0]->getScanNumber());
        Assert::AreEqual(allPsmsArray[0]->getScanNumber(), allPsmsArray2[0]->getScanNumber());
        
        Assert::AreEqual("QQQGGGG", allPsmsArray2[0]->getBaseSequence());
        Assert::AreEqual(allPsmsArray[0]->getBaseSequence(), allPsmsArray2[0]->getBaseSequence());
        
        delete DeconvolutionMassTolerance2;
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters2' statement was not added since CommonParameters2
        //was passed to a method or constructor. Handle memory management manually.
        delete DeconvolutionMassTolerance;
        //C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease2' statement was not added since protease2
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete searchModes' statement was not added since searchModes
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete productMassTolerance' statement was not added since
        //productMassTolerance was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete myMsDataFile' statement was not added since myMsDataFile
        //was passed to a method or constructor. Handle memory management manually.
    }

    void SearchEngineTests::TestClassicSemiProteolysis()
    {
        auto variableModifications = std::vector<Modification*>();
        auto fixedModifications = std::vector<Modification*>();
        auto localizeableModifications = std::vector<Modification*>();
        std::unordered_map<Modification*, unsigned short> modsDictionary;
        for (auto mod : fixedModifications)
        {
            modsDictionary.emplace(mod, 0);
        }
        int ii = 1;
        for (auto mod : variableModifications)
        {
            modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
            ii++;
        }
        for (auto mod : localizeableModifications)
        {
            modsDictionary.emplace(mod, static_cast<unsigned short>(ii));
            ii++;
        }
        std::vector<ProteolysisProduct*> protprod = {new ProteolysisProduct(9, 21, "chain")};
        auto proteinList = std::vector<Protein*> {new Protein("MGGGGGMKNNNQQQGGGGKLKGKKNKKGN", "hello", nullptr, nullptr, nullptr, protprod)};
        
        std::vector<DigestionMotif*> motifs1 = {new DigestionMotif("G", nullptr, 1, nullptr)};
        auto protease = new Protease("semi-Trypsin", CleavageSpecificity::Semi, nullptr, nullptr, motifs1);
        ProteaseDictionary::Dictionary->Add(protease->Name, protease);
        DigestionParams *digestParams = new DigestionParams(protease: protease->Name, minPeptideLength: 2);
        
        //expect NNNQQQ, NNNQQ, NNNQ, NNN, NN and LK, KLK
        std::unordered_map<std::string, bool> found =
            {
                {"NNNQQQ", false},
                {"NNNQQ", false},
                {"NNNQ", false},
                {"NNN", false},
                {"NN", false},
                {"LK", false},
                {"KLK", false}
            };
        std::vector<PeptideWithSetModifications*> PWSMs = proteinList[0]->Digest(digestParams, std::vector<Modification*>(), modsDictionary.Keys->ToList());
        for (auto PWSM : PWSMs)
        {
            bool b;
            std::unordered_map<std::string, bool>::const_iterator found_iterator = found.find(PWSM.BaseSequence);
            if (found_iterator != found.end())
            {
                b = found_iterator->second;
                found[PWSM->BaseSequence] = true;
            }
            else
            {
                b = found_iterator->second;
            }
        }
        for (auto kvp : found)
        {
            Assert::IsTrue(kvp.second);
        }
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete digestParams' statement was not added since digestParams
        //was passed to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease
        //was passed to a method or constructor. Handle memory management manually.
    }
    
    void SearchEngineTests::TestClassicSearchOneNterminalModifiedPeptideOneScan()
    {
        //Create the modification and select the motif
        ModificationMotif motif;
        ModificationMotif::TryGetMotif("A", motif);
        Modification *mod = new Modification(_originalId: "acetylation", _modificationType: "testModType", _target: motif,
                                             _chemicalFormula: ChemicalFormula::ParseFormula("C2H2O1"),
                                             _locationRestriction: "Anywhere.");
        
        //Create the protein and add the mod
        std::vector<Modification*> modlist = {mod};
        int modPositionInProtein = 1;
        auto proteinListForSearch = new Protein("AAAGSAAVSGAGTPVAGPTGR", nullptr, oneBasedModifications: std::unordered_map<int, std::vector<Modification*>>    {
                {modPositionInProtein, modlist}
            });
        std::vector<Protein*> proteins = {proteinListForSearch};
        
        //Prepare the scan
        std::vector<double> mz = {101.814659, 101.851852, 101.856758, 101.920807, 101.957397, 101.993172, 102.029053, 102.035851, 102.065628, 102.102516, 102.123955, 102.129654, 102.136391, 102.144806, 102.173798, 102.182594, 102.19574, 102.201569, 102.208534, 102.227333, 102.231323, 102.237877, 102.266922, 110.071281, 124.284775, 125.057686, 129.101776, 136.075348, 141.059647, 141.067139, 147.112045, 152.296143, 157.095963, 159.321686, 159.391266, 159.461609, 159.505859, 159.530792, 159.567841, 159.673218, 159.688232, 159.708939, 159.717117, 159.744415, 159.883804, 169.132446, 175.118637, 175.639145, 185.091217, 185.629135, 185.649612, 198.086334, 202.961533, 225.122726, 230.112839, 232.139282, 247.539169, 256.128204, 257.129761, 269.123596, 283.467743, 287.133331, 304.127319, 306.832855, 307.01059, 307.281647, 313.149384, 314.151093, 329.183044, 358.160583, 360.700592, 365.218567, 370.208008, 372.185425, 373.445496, 382.170685, 383.1745, 400.180878, 401.178406, 415.016876, 425.206726, 430.238647, 431.253052, 438.198273, 439.226074, 442.213654, 453.207336, 454.210938, 467.186005, 470.226654, 471.218048, 472.219635, 487.260254, 488.259033, 498.980499, 524.244202, 525.244385, 528.277832, 541.27124, 542.254883, 543.258118, 558.296448, 559.29895, 602.023743, 637.564087, 638.122498, 638.68103, 639.24469, 641.322876, 642.328674, 651.320923, 665.653137, 665.812622, 678.225586, 678.837524, 684.52533, 699.105225, 702.819702, 710.320374, 718.544617, 731.125793, 731.804626, 732.479858, 732.635559, 745.458618, 754.418152, 755.42041, 770.024719, 835.083557, 855.465454, 883.47406, 884.485962, 885.430664, 894.471252, 912.484192, 913.487305, 983.525391, 984.515137, 1022.541321, 1040.545166, 1041.545166, 1109.568115, 1110.565918, 1122.717896, 1127.574707, 1128.582764, 1143.576172, 1192.367432, 1209.613403, 1226.649048, 1227.651123, 1255.27356, 1270.274414, 1273.796753, 1297.668091, 1298.676636, 1520.323364, 1541.422729, 1546.209839, 1639.540283, 1669.039673, 1670.780518, 1671.45459, 1671.927368, 1672.41272, 1673.799194, 1674.881836, 1687.070557};
        std::vector<double> intensities = {2084.605469, 1503.049316, 1508.572144, 1997.866089, 3907.81958, 3168.024902, 2187.750488, 2006.690186, 1900.748047, 2346.96875, 5075.943359, 1605.881958, 3754.241699, 2652.021484, 2389.109375, 1921.068115, 7624.829102, 2316.127441, 4926.003906, 4061.111328, 4849.558105, 1925.951172, 2153.752441, 1949.42981, 2058.943604, 1687.470215, 1928.76416, 4824.449219, 1376.267334, 1439.787109, 2010.894897, 1483.543213, 1929.985596, 3110.926025, 3372.582275, 2304.003174, 2679.970215, 2560.099854, 1692.449585, 4134.439453, 3104.67627, 3688.150146, 3676.4646, 2814.727051, 1772.865723, 2332.302002, 4913.470703, 1477.579102, 19578.2793, 2073.369141, 2618.062988, 4346.097168, 1499.773193, 1665.242065, 2426.635742, 7056.563477, 1847.738281, 34847.85547, 2724.011963, 6279.092773, 1623.552734, 6072.998535, 3376.6521, 1946.103149, 2381.508301, 1723.251953, 43135.01953, 2463.636719, 2124.577393, 3385.583008, 1918.600098, 3191.296387, 1865.107422, 3907.073975, 2091.009766, 47135.79688, 5113.383301, 19901.54883, 1716.9729, 1691.421753, 2594.325928, 6075.093262, 1948.982544, 1560.634033, 2563.348145, 1857.768311, 96429.07813, 12289.67578, 2152.950928, 1781.416992, 49654.84766, 7577.02832, 19761.38867, 3075.865967, 2130.962158, 6758.932617, 2022.507935, 4839.808105, 2863.790039, 67360.53906, 11148.44824, 23283.18945, 4307.828125, 1818.466187, 4734.982422, 2687.210449, 4251.597168, 3365.632813, 14722.28027, 3734.265869, 1774.239624, 1827.377075, 2492.476074, 9175.493164, 2313.220215, 2536.213867, 2223.19165, 2145.397461, 2189.492676, 3675.263184, 2770.946289, 4485.819336, 2932.245605, 1713.99939, 1857.43042, 45807.97266, 10789.84961, 1848.404297, 1901.754517, 10514.50977, 14093.94141, 5819.866699, 6942.87207, 3736.169189, 34573.74219, 11260.63867, 13557.64746, 3463.163818, 7701.758301, 39518.83203, 15681.95898, 13644.61426, 7358.266602, 5733.412598, 46497.80469, 20523.03516, 1688.138916, 1865.104248, 2219.646729, 9009.290039, 4519.379395, 2556.661621, 1627.402466, 2145.759766, 4486.545898, 3421.69873, 1836.025391, 1980.848999, 2378.1521, 5128.462402, 2796.475586, 2296.582764, 8305.148438, 9054.726563, 6919.573242, 3180.789063, 2163.666504, 1787.661133};

        std::vector<MsDataScan*> msDataScanArray(1);
        MzSpectrum tempVar(mz, intensities, false);
        MzRange tempVar2(400, 1600);
        msDataScanArray[0] = new MsDataScan(massSpectrum: &tempVar, oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true, polarity: Polarity::Positive,
                                            retentionTime: 10.0, scanWindowRange: &tempVar2, scanFilter: "f", mzAnalyzer: MZAnalyzerType::Orbitrap,
                                            totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: nullptr, nativeId: "scan=" + std::to_string(1));
        
        std::vector<Ms2ScanWithSpecificMass*> scans(1);
        double precursorMz = 884.4541;
        int precursorCharge = 2;
        std::string filePath = "path";
        CommonParameters tempVar3();
        scans[0] = new Ms2ScanWithSpecificMass(msDataScanArray[0], precursorMz, precursorCharge, filePath, &tempVar3);
        
        //set digestion parameters and digest the protein
        DigestionParams *testDigestionParams = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1,
                                                                   initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain);
        
        //set the common parameters
        Tolerance *productTolerance = new PpmTolerance(20);
        Tolerance *precursorTolerance = new PpmTolerance(10);
        DissociationType *d = DissociationType::CID;
        CommonParameters *commonParameters = new CommonParameters(, d, , , , , , , , , , , , , , , productTolerance,
                                                                  precursorTolerance, , , testDigestionParams);
        
        //Search the scan against the protein
        std::vector<PeptideSpectralMatch*> allPsmsArray(1);
        ClassicSearchEngine tempVar4(allPsmsArray, scans, new std::vector<Modification*>(), new std::vector<Modification*>(), proteins,
                                     new SinglePpmAroundZeroSearchMode(20), new CommonParameters(, DissociationType::CID),
                                     new std::vector<std::string>());
        MetaMorpheusEngineResults *engineResults = (&tempVar4)->Run();
        
        //Process the results
        std::vector<MatchedFragmentIon*> matchedFragmentIonList = allPsmsArray.SelectMany([&] (std::any p)   {
                p::MatchedFragmentIons;
            }).ToList();
        std::vector<double> neutralFragmentMasses = matchedFragmentIonList.Select([&] (std::any m)  {
                m::NeutralTheoreticalProduct::NeutralMass;
            }).ToList();
        std::vector<double> massToCharges;
        for (auto mass : neutralFragmentMasses)
        {
            massToCharges.push_back(mass.ToMz(1));
        }
        Assert::AreEqual(20, massToCharges.size()());
        
        delete commonParameters;
        delete precursorTolerance;
        delete productTolerance;
        delete testDigestionParams;
        delete proteinListForSearch;
        delete mod;
    }
#endif
    
}
