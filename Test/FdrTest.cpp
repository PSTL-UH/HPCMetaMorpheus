#include "FdrTest.h"
#include "../EngineLayer/PrecursorSearchModes/DotMassDiffAcceptor.h"
#include "../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "TestDataFile.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../TaskLayer/SearchTask/SearchParameters.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../TaskLayer/MetaMorpheusTask.h"
#include "../TaskLayer/SearchTask/SearchTask.h"
using namespace TaskLayer;

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

#include "../EngineLayer/ClassicSearch/ClassicSearchEngine.h"
using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;

#include  "../EngineLayer/FdrAnalysis/FdrAnalysisEngine.h"
using namespace EngineLayer::FdrAnalysis;

#include "../EngineLayer/Indexing/IndexingEngine.h"
#include "../EngineLayer/Indexing/IndexingResults.h"
using namespace EngineLayer::Indexing;

#include "../EngineLayer/ModernSearch/ModernSearchEngine.h"
using namespace EngineLayer::ModernSearch;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
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
    GlobalVariables::GlobalVariables_init();
        
    std::cout << ++i << ". FdrTestMethod" << std::endl;
    Test::FdrTest::FdrTestMethod();

    std::cout << ++i << ". TestDeltaValues" << std::endl;
    Test::FdrTest::TestDeltaValues();

    return 0;
}

namespace Test
{

    void FdrTest::FdrTestMethod()
    {
        auto tempVar = new PpmTolerance(5);
        std::vector<double> vd = {0, 1.0029};
        MassDiffAcceptor *searchModes = new DotMassDiffAcceptor("", vd, tempVar);
        std::vector<std::string> nestedIds;
        
        Protein *p = new Protein("MNKNNKNNNKNNNNK", "");
        DigestionParams *digestionParams = new DigestionParams("trypsin");
        std::vector<Modification*> vm1, vm2;
        auto digested = p->Digest(digestionParams, vm1, vm2);
        
        PeptideWithSetModifications *pep1 = digested[0];
        PeptideWithSetModifications *pep2 = digested[1];
        PeptideWithSetModifications *pep3 = digested[2];
        PeptideWithSetModifications *pep4 = digested[3];

        std::vector<PeptideWithSetModifications*> pepvec = {pep1, pep2, pep3};
        TestDataFile *t = new TestDataFile(pepvec );
        
        MsDataScan *mzLibScan1 = t->GetOneBasedScan(2);
        auto tempVar2 = new CommonParameters();
        std::vector<MassSpectrometry::IsotopicEnvelope*> nEFs;
        Ms2ScanWithSpecificMass *scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, Chemistry::ClassExtensions::ToMz(pep1->getMonoisotopicMass(), 1), 1, "", tempVar2, nEFs);
        std::vector<MatchedFragmentIon*> vMFI;
        PeptideSpectralMatch *psm1 = new PeptideSpectralMatch(pep1, 0, 3, 0, scan1, digestionParams, vMFI);
        
        MsDataScan *mzLibScan2 = t->GetOneBasedScan(4);
        auto tempVar3 = new CommonParameters();
        std::vector<MassSpectrometry::IsotopicEnvelope*> nEFs2;
        Ms2ScanWithSpecificMass *scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, Chemistry::ClassExtensions::ToMz(pep2->getMonoisotopicMass(), 1), 1, "", tempVar3, nEFs2);
        std::vector<MatchedFragmentIon*> vMFI2;
        PeptideSpectralMatch *psm2 = new PeptideSpectralMatch(pep2, 1, 2, 1, scan2, digestionParams, vMFI2);
        
        MsDataScan *mzLibScan3 = t->GetOneBasedScan(6);
        auto tempVar4 = new CommonParameters();
        std::vector<MassSpectrometry::IsotopicEnvelope*> nEFs3;
        Ms2ScanWithSpecificMass *scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, Chemistry::ClassExtensions::ToMz(pep3->getMonoisotopicMass(), 1), 1, "", tempVar4, nEFs3);
        std::vector<MatchedFragmentIon*> vMFI3;
        PeptideSpectralMatch *psm3 = new PeptideSpectralMatch(pep3, 0, 1, 2, scan3, digestionParams, vMFI3);

        std::vector<MatchedFragmentIon*> vMFI4;
        psm3->AddOrReplace(pep4, 1, 1, true, vMFI4);
        
        auto newPsms = std::vector<PeptideSpectralMatch*> {psm1, psm2, psm3};
        for (auto psm : newPsms)
        {
            psm->ResolveAllAmbiguities();
        }
        
        
        CommonParameters *cp = new CommonParameters( "", DissociationType::HCD, true, true, 3, 12, true, false, 1, 5, 200, 0.01, false, true, false, true);
        
        FdrAnalysisEngine *fdr = new FdrAnalysisEngine(newPsms, searchModes->getNumNotches(), cp, nestedIds);
        
        fdr->Run();
        
        Assert::AreEqual(2, searchModes->getNumNotches());
        Assert::AreEqual(0, newPsms[0]->getFdrInfo()->getCumulativeDecoyNotch());
        Assert::AreEqual(1, newPsms[0]->getFdrInfo()->getCumulativeTargetNotch());
        Assert::AreEqual(0, newPsms[1]->getFdrInfo()->getCumulativeDecoyNotch());
        Assert::AreEqual(1, newPsms[1]->getFdrInfo()->getCumulativeTargetNotch());
        Assert::AreEqual(0, newPsms[2]->getFdrInfo()->getCumulativeDecoyNotch());
        Assert::AreEqual(1, newPsms[2]->getFdrInfo()->getCumulativeTargetNotch());
        
        Assert::AreEqual(0, newPsms[0]->getFdrInfo()->getCumulativeDecoy());
        Assert::AreEqual(1, newPsms[0]->getFdrInfo()->getCumulativeTarget());
        Assert::AreEqual(0, newPsms[1]->getFdrInfo()->getCumulativeDecoy());
        Assert::AreEqual(2, newPsms[1]->getFdrInfo()->getCumulativeTarget());
        Assert::AreEqual(0, newPsms[2]->getFdrInfo()->getCumulativeDecoy());
        Assert::AreEqual(3, newPsms[2]->getFdrInfo()->getCumulativeTarget());
        
        delete fdr;
        delete cp;
        delete psm3;
        delete scan3;
        delete psm2;
        delete scan2;
        delete scan1;
        delete psm1;
        delete t;
        delete digestionParams;
        delete p;
        delete searchModes;
    }


    void FdrTest::TestDeltaValues()
    {
        auto tempVar = new DigestionParams("trypsin", 2, 5);
        CommonParameters *cp = new CommonParameters("", DissociationType::HCD, true, true, 3, 12, true, false, 1, 1, 200, 0.01, false, true, true,
                                                    false, nullptr, nullptr, nullptr, -1, tempVar);
        
        SearchParameters *sp = new SearchParameters();
        sp->setMassDiffAcceptorType(MassDiffAcceptorType::Exact);
        auto allmods = GlobalVariables::getAllModsKnown();
        auto listofmodsvar = cp->getListOfModsVariable();
        std::vector<Modification*> variableModifications;
        for ( auto b: allmods ) {
            if ( std::find(listofmodsvar->begin(), listofmodsvar->end(), std::make_tuple(b->getModificationType(), b->getIdWithMotif())) != listofmodsvar->end() ) {
                variableModifications.push_back(b);
            }
        }
        std::vector<Modification*> fixedModifications;
        auto listofmodsfixed =  cp->getListOfModsFixed();
        for ( auto b: allmods ) {
            if ( std::find(listofmodsfixed->begin(), listofmodsfixed->end(), std::make_tuple(b->getModificationType(), b->getIdWithMotif())) != listofmodsfixed->end() ) {
                fixedModifications.push_back(b);
            }
        }
        
        // Generate data for files
        Protein *TargetProtein1 = new Protein("TIDEANTHE", "accession1");
        Protein *TargetProtein2 = new Protein("TIDELVE", "accession2");
        Protein *TargetProtein3 = new Protein("TIDENIE", "accession3");
        Protein *TargetProteinLost = new Protein("PEPTIDEANTHE", "accession4");
        Protein *DecoyProteinFound = new Protein("PETPLEDQGTHE", "accessiond", "",
                                                 std::vector<std::tuple<std::string, std::string>>(),
                                                 std::unordered_map<int, std::vector<Modification*>>(),
                                                 std::vector<ProteolysisProduct*>(),
                                                 "", "", true);

        auto tp1 = TargetProtein1->Digest(cp->getDigestionParams(), fixedModifications, variableModifications)[0];
        auto tp2 = TargetProtein2->Digest(cp->getDigestionParams(), fixedModifications, variableModifications)[0];
        auto tp3 = TargetProtein3->Digest(cp->getDigestionParams(), fixedModifications, variableModifications)[0];
        auto dpf = DecoyProteinFound->Digest(cp->getDigestionParams(), fixedModifications, variableModifications)[0];

        std::vector<PeptideWithSetModifications*> vPWSM = {tp1, tp2, tp3, dpf};
        MsDataFile *myMsDataFile = new TestDataFile( vPWSM );
        
        std::vector<Protein*> proteinList = {TargetProtein1, TargetProtein2, TargetProtein3, TargetProteinLost, DecoyProteinFound};
        
        auto searchModes = new SinglePpmAroundZeroSearchMode(5);
        
        bool DoPrecursorDeconvolution = true;
        bool UseProvidedPrecursorInfo = true;
        double DeconvolutionIntensityRatio = 4;
        int DeconvolutionMaxAssumedChargeState = 10;
        Tolerance *DeconvolutionMassTolerance = new PpmTolerance(5);
        
        auto tempVar2 = new CommonParameters();
        auto listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", tempVar2);
        std::sort( listOfSortedms2Scans.begin(), listOfSortedms2Scans.end(), [&] (auto l, auto r)  {
                return l->getPrecursorMass() < r->getPrecursorMass();
            });
        
        //check better when using delta
        std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
        std::vector<std::string> vs, vs1, vs2, vs3;

        ClassicSearchEngine tempVar3(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications,
                                     proteinList, searchModes,  cp, vs);
        (&tempVar3)->Run();

        auto indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1,
                                              DecoyType::None, cp, 30000, false,     vs1, vs2);
        auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        MassDiffAcceptor *massDiffAcceptor = SearchTask::GetMassDiffAcceptor(cp->getPrecursorMassTolerance(),
                                                                             sp->getMassDiffAcceptorType(),
                                                                             sp->getCustomMdac());
        
        std::vector<PeptideSpectralMatch*> allPsmsArrayModern(listOfSortedms2Scans.size());

        auto v1 = indexResults->getPeptideIndex();
        auto v2 = indexResults->getFragmentIndex();
        ModernSearchEngine tempVar4(allPsmsArrayModern, listOfSortedms2Scans, v1, v2, 0, cp,
                                    massDiffAcceptor, 0, vs3 );
        (&tempVar4)->Run();

        std::vector<std::string> vs4, vs5, vs6, vs7;
        FdrAnalysisEngine tempVar5(allPsmsArray, 1, cp, vs4);
        FdrAnalysisResults *fdrResultsClassicDelta = static_cast<FdrAnalysisResults*>((&tempVar5)->Run());

        FdrAnalysisEngine tempVar6(allPsmsArrayModern, 1, cp, vs5);
        FdrAnalysisResults *fdrResultsModernDelta = static_cast<FdrAnalysisResults*>((&tempVar6)->Run());
        Assert::AreEqual(fdrResultsClassicDelta->getPsmsWithin1PercentFdr(), 3);
        Assert::AreEqual(fdrResultsModernDelta->getPsmsWithin1PercentFdr(), 3);
        
        auto tempVar7 = new DigestionParams ("trypsin", 2,  5);
        cp = new CommonParameters();
        cp->setDigestionParams(tempVar7);
        
        //check worse when using score
        FdrAnalysisEngine tempVar8(allPsmsArray, 1, cp, vs6);
        FdrAnalysisResults *fdrResultsClassic = static_cast<FdrAnalysisResults*>((&tempVar8)->Run());
        FdrAnalysisEngine tempVar9(allPsmsArray, 1, cp, vs7);
        FdrAnalysisResults *fdrResultsModern = static_cast<FdrAnalysisResults*>((&tempVar9)->Run());
        Assert::AreEqual(fdrResultsClassic->getPsmsWithin1PercentFdr(), 0);
        Assert::AreEqual(fdrResultsModern->getPsmsWithin1PercentFdr(), 0);
        
        //check that when delta is bad, we used the score
        // Generate data for files
        Protein *DecoyProtein1 = new Protein("TLEDAGGTHE", "accession1d", "", 
                                             std::vector<std::tuple<std::string, std::string>>(),
                                             std::unordered_map<int, std::vector<Modification*>>(),
                                             std::vector<ProteolysisProduct*>(),
                                             "", "", true);
        Protein *DecoyProtein2 = new Protein("TLEDLVE", "accession2d", "", 
                                             std::vector<std::tuple<std::string, std::string>>(),
                                             std::unordered_map<int, std::vector<Modification*>>(),
                                             std::vector<ProteolysisProduct*>(),
                                             "", "", true);
        Protein *DecoyProtein3 = new Protein("TLEDNIE", "accession3d", "", 
                                             std::vector<std::tuple<std::string, std::string>>(),
                                             std::unordered_map<int, std::vector<Modification*>>(),
                                             std::vector<ProteolysisProduct*>(),
                                             "", "", true);
        Protein *DecoyProteinShiny = new Protein("GGGGGG", "accessionShinyd", "", 
                                                 std::vector<std::tuple<std::string, std::string>>(),
                                                 std::unordered_map<int, std::vector<Modification*>>(),
                                                 std::vector<ProteolysisProduct*>(),
                                                 "", "", true);

        tp1 = TargetProtein1->Digest(cp->getDigestionParams(), fixedModifications, variableModifications)[0];
        tp2 = TargetProtein2->Digest(cp->getDigestionParams(), fixedModifications, variableModifications)[0];
        tp3 = TargetProtein3->Digest(cp->getDigestionParams(), fixedModifications, variableModifications)[0];
        auto dps = DecoyProteinShiny->Digest(cp->getDigestionParams(), fixedModifications, variableModifications)[0];

        std::vector<PeptideWithSetModifications*> vPWSM2 = {tp1, tp2, tp3, dps};
        myMsDataFile = new TestDataFile(vPWSM2);
        
        proteinList = {TargetProtein1, DecoyProtein1, TargetProtein2, DecoyProtein2, TargetProtein3, DecoyProtein3, DecoyProteinShiny};
        
        auto tempVar10 = new CommonParameters();
        listOfSortedms2Scans = MetaMorpheusTask::GetMs2Scans(myMsDataFile, "", tempVar10);
        std::sort( listOfSortedms2Scans.begin(), listOfSortedms2Scans.end(), [&] (auto l, auto r)  {
                return l->getPrecursorMass() < r->getPrecursorMass();
            });
        
        //check no change when using delta
        allPsmsArray = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());
        std::vector<std::string> vs10, vs11, vs12, vs13;
        ClassicSearchEngine tempVar11(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes,
                                      cp, vs10);
        (&tempVar11)->Run();
        
        auto tempVar12 = new DigestionParams( "trypsin", 2,  5);
        cp = new CommonParameters("", DissociationType::HCD, true, true, 3, 12, true, false, 1, 1, 200, 0.01, false, true, true,
                                  false, nullptr, nullptr, nullptr, -1, tempVar12);

        indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::None, cp, 30000, false,
                                         vs11, vs12);
        indexResults = static_cast<IndexingResults*>(indexEngine->Run());
        massDiffAcceptor = SearchTask::GetMassDiffAcceptor(cp->getPrecursorMassTolerance(), sp->getMassDiffAcceptorType(),
                                                           sp->getCustomMdac());
        allPsmsArrayModern = std::vector<PeptideSpectralMatch*>(listOfSortedms2Scans.size());

        auto v3 = indexResults->getPeptideIndex();
        auto v4 = indexResults->getFragmentIndex();
        ModernSearchEngine tempVar13(allPsmsArrayModern, listOfSortedms2Scans, v3, v4, 0,
                                     cp, massDiffAcceptor, 0, vs13);
        (&tempVar13)->Run();
        
        std::vector<std::string> vs14, vs15, vs16, vs17;
        FdrAnalysisEngine tempVar14(allPsmsArray, 1, cp, vs14);
        fdrResultsClassicDelta = static_cast<FdrAnalysisResults*>((&tempVar14)->Run());
        FdrAnalysisEngine tempVar15(allPsmsArrayModern, 1, cp, vs15);
        fdrResultsModernDelta = static_cast<FdrAnalysisResults*>((&tempVar15)->Run());
        Assert::AreEqual(fdrResultsClassicDelta->getPsmsWithin1PercentFdr(),3);
        Assert::AreEqual(fdrResultsModernDelta->getPsmsWithin1PercentFdr(), 3);
        
        auto tempVar16 = new DigestionParams("trypsin", 2, 5);
        cp = new CommonParameters();
        cp->setDigestionParams(tempVar16);
        
        //check no change when using score
        FdrAnalysisEngine tempVar17(allPsmsArray, 1, cp, vs16);
        fdrResultsClassic = static_cast<FdrAnalysisResults*>((&tempVar17)->Run());
        FdrAnalysisEngine tempVar18(allPsmsArrayModern, 1, cp, vs17);
        fdrResultsModern = static_cast<FdrAnalysisResults*>((&tempVar18)->Run());
        Assert::AreEqual(fdrResultsClassic->getPsmsWithin1PercentFdr(), 3);
        Assert::AreEqual(fdrResultsModern->getPsmsWithin1PercentFdr(), 3);
        
        delete tempVar;
        delete tempVar2;
        delete tempVar7;
        delete tempVar16;
        delete sp;
        delete DecoyProteinShiny;
        delete DecoyProtein3;
        delete DecoyProtein2;
        delete DecoyProtein1;
        delete indexEngine;
        delete DeconvolutionMassTolerance;
        delete searchModes;
        delete myMsDataFile;
        delete DecoyProteinFound;
        delete TargetProteinLost;
        delete TargetProtein3;
        delete TargetProtein2;
        delete TargetProtein1;
    }
}
