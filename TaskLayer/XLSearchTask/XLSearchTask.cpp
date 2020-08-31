#include "XLSearchTask.h"
#include "XLSearchParameters.h"
#include "../../EngineLayer/CommonParameters.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"
#include "../MyTaskResults.h"
#include "../MyFileManager.h"
#include "../../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../../EngineLayer/GlobalVariables.h"
#include "../../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../../EngineLayer/CrosslinkSearch/PsmCrossType.h"

#include <sstream>

#include "MassSpectrometry/Enums/DissociationType.h"
#include "UsefulProteomicsDatabases/DecoyType.h"

#include "pepXML/pepXML_v120.h"
#include "stringhelper.h"
#include "time.h"


using namespace EngineLayer;
using namespace EngineLayer::CrosslinkSearch;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace MzLibUtil;
using namespace EngineLayer::FdrAnalysis;
using namespace Proteomics::Fragmentation;

#ifdef TIMING_INFO
#include <sys/time.h>

static double timediff (struct timeval t1, struct timeval t2)
{
    double elapsedtime;
    elapsedtime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedtime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

    return elapsedtime/1000;                            //ms to sec
}
#endif

namespace TaskLayer
{

    XLSearchTask::XLSearchTask() : MetaMorpheusTask(MyTask::XLSearch)
    {
        std::string taskDescr = "";
        DissociationType dissType = DissociationType::HCD;
        int topNpeaks = 200;
        double minRatio = 0.01;
        bool trimMs1Peaks = false;
        bool trimMsMsPeaks = true;
        bool useDeltaScore = false;
        bool calculateEValue = false;
        
        auto tempVar = new EngineLayer::CommonParameters( taskDescr, dissType, true, true, 3, 12, true, false,
                                               1, 3, topNpeaks,
                                               minRatio, trimMs1Peaks, trimMsMsPeaks, useDeltaScore,
                                               calculateEValue, nullptr, new PpmTolerance(10));
        setCommonParameters(tempVar);
        
        auto tempVar2 = new TaskLayer::XlSearchParameters();
        setXlSearchParameters(tempVar2);
    }

    XLSearchTask::XLSearchTask(std::string tomlFile ) : MetaMorpheusTask(MyTask::XLSearch)
    {
        //Check ending that it is a toml file.
        
        Toml trw;
        toml::Value toml_value = trw.tomlReadFile(tomlFile);
        toml::Value* fileParameters = trw.getValue(toml_value, "XlSearchParameters");
        toml::Table tomlTable = fileParameters->as<toml::Table>();              

        auto xlParams = new TaskLayer::XlSearchParameters();

        // parse toml file and set the values
        for (auto const& keyValuePair : tomlTable)
        {
            // we're using the name of the variable here and not a fixed string
            // in case the variable name changes at some point
            if (keyValuePair.first == "DecoyType") {
                auto val = keyValuePair.second.as<std::string>();
                xlParams->setDecoyType(DecoyTypeFromString(val ));
            }
            else if ( keyValuePair.first == "CrosslinkerType") {
                auto val = keyValuePair.second.as<std::string>();
                xlParams->setCrosslinkerType(CrosslinkerTypeFromString(val));
            }
            else if ( keyValuePair.first == "CrosslinkSearchTopNum") {
                xlParams->setCrosslinkSearchTopNum(keyValuePair.second.as<int>() );
            }
            else if  ( keyValuePair.first == "CrosslinkerResidues" ) {
                xlParams->setCrosslinkerResidues(keyValuePair.second.as<std::string>() );
            }
            else if ( keyValuePair.first == "CrosslinkerResidues2" ) {
                xlParams->setCrosslinkerResidues2(keyValuePair.second.as<std::string>() );
            }
            else if ( keyValuePair.first == "IsCleavable" ) {
                xlParams->setIsCleavable ( keyValuePair.second.as<bool>() );
            }
            else if ( keyValuePair.first == "RestrictToTopNHits" ) {
                xlParams->setRestrictToTopNHits(keyValuePair.second.as<bool>() ); 
            }
            else if ( keyValuePair.first == "WriteOutputForPercolator" ) {
                xlParams->setWriteOutputForPercolator(keyValuePair.second.as<bool>() ) ;
            }
            else if ( keyValuePair.first == "WritePepXml") {
                xlParams->setWritePepXml( keyValuePair.second.as<bool>() );
            }
            else if ( keyValuePair.first == "XlQuench_H2O") {
                xlParams->setXlQuench_H2O(keyValuePair.second.as<bool>() );
            }
            else if ( keyValuePair.first == "XlQuench_Tris") {
                xlParams->setXlQuench_Tris(keyValuePair.second.as<bool>() );
            }
            else if ( keyValuePair.first == "XlQuench_NH2") {
                xlParams->setXlQuench_NH2(keyValuePair.second.as<bool>() );
            }
        }
        setXlSearchParameters(xlParams);

        // Do we need to read the common parameters as well? Probably yes.
        std::string taskDescr = "";
        auto tempVar = new EngineLayer::CommonParameters( tomlFile, taskDescr );
        setCommonParameters(tempVar);                
    }

    void XLSearchTask::writeTomlConfig(std::string &filename, std::ofstream &tomlFd )
    {
        if ( !tomlFd.is_open() ) {
            tomlFd.open(filename );
            if ( !tomlFd.is_open() ) {
                std::cout << "XLSearchTask: Could not open file " << filename << std::endl;
                return;
            }
        }

        toml::Value v;
        std::string key = "TaskType", value = "XLSearch";
        v.set ( key, value);
        tomlFd << v;
        
        XlSearchParameters *xlparams = getXlSearchParameters();
        toml::Table search_params;

        auto tvar1 = xlparams->getDecoyType();
        search_params["DecoyType"] = DecoyTypeToString(tvar1);
        auto tvar2 = xlparams->getCrosslinkerType();
        search_params["CrosslinkerType"] = CrosslinkerTypeToString(tvar2);
        search_params["CrosslinkSearchTopNum"] = xlparams->getCrosslinkSearchTopNum();
        //search_params[""] = xlparams->std::string getCrosslinkerName();
        //search_params[""] = xlparams->std::optional<double> getCrosslinkerTotalMass();
        //search_params[""] = xlparams->std::optional<double> getCrosslinkerShortMass();
        //search_params[""] = xlparams->std::optional<double> getCrosslinkerLongMass();
        //search_params[""] = xlparams->std::optional<double> getCrosslinkerLoopMass();
        search_params["CrosslinkerResidues"] = xlparams->getCrosslinkerResidues();
        search_params["CrosslinkerResidues2"] = xlparams->getCrosslinkerResidues2();
        //search_params[""] = xlparams->std::optional<double> getCrosslinkerDeadEndMassH2O();
        //search_params[""] = xlparams->std::optional<double> getCrosslinkerDeadEndMassNH2();
        //search_params[""] = xlparams->std::optional<double> getCrosslinkerDeadEndMassTris();
        search_params["IsCleavable"] = xlparams->getIsCleavable();
        search_params["RestrictToTopNHits"] = xlparams->getRestrictToTopNHits();
        search_params["WriteOutputForPercolator"] = xlparams->getWriteOutputForPercolator();
        search_params["WritePepXml"] = xlparams->getWritePepXml();
        search_params["XlQuench_H2O"] = xlparams->getXlQuench_H2O();
        search_params["XlQuench_Tris"] = xlparams->getXlQuench_Tris();
        search_params["XlQuench_NH2"] = xlparams->getXlQuench_NH2();

        tomlFd << std::endl;
        tomlFd << "[XlSearchParameters]" << std::endl;
        tomlFd << search_params;
        
        // Now write the generic parameters
        MetaMorpheusTask::writeTomlConfig(filename, tomlFd );
        if ( tomlFd.is_open() ) {
            tomlFd.close();
        }
        return;
    }

    
    TaskLayer::XlSearchParameters *XLSearchTask::getXlSearchParameters() const
    {
        return privateXlSearchParameters;
    }
    
    void XLSearchTask::setXlSearchParameters(TaskLayer::XlSearchParameters *value)
    {
        if (privateXlSearchParameters != nullptr ) {
            delete privateXlSearchParameters;
        }
        privateXlSearchParameters = value;
    }
    
    MyTaskResults *XLSearchTask::RunSpecific(const std::string &OutputFolder,
                                             std::vector<DbForTask*> &dbFilenameList,
                                             std::vector<std::string> &currentRawFileList,
                                             const std::string &taskId,
                                             std::vector<FileSpecificParameters*> &fileSettingsList)
    {
        myTaskResults = new MyTaskResults(this);
        std::vector<CrosslinkSpectralMatch*> allPsms;
        
        std::vector<Modification*> variableModifications;
        std::vector<Modification*> fixedModifications;
        std::vector<std::string> localizeableModificationTypes;

#ifdef TIMING_INFO
        struct timeval t1, t1e;
        struct timeval t2, t2e;
        struct timeval t3, t3e;
        struct timeval t4, t4e;
        struct timeval t5, t5e;
        struct timeval t6, t6e;
        struct timeval t7, t7e;
        struct timeval t8, t8e;
        struct timeval t9, t9e;
        double t3total=0.0, t4total=0.0, t5total=0.0, t6total=0.0;

        gettimeofday (&t1, NULL);
#endif
        LoadModifications(taskId, variableModifications, fixedModifications, localizeableModificationTypes);
#ifdef TIMING_INFO
        gettimeofday (&t1e, NULL);
#endif
        
        // load proteins
#ifdef TIMING_INFO
        gettimeofday (&t2, NULL);
#endif
        std::vector<Protein*> proteinList = LoadProteins(taskId, dbFilenameList, true,
                                                         getXlSearchParameters()->getDecoyType(),
                                                         localizeableModificationTypes,
                                                         getCommonParameters());
#ifdef TIMING_INFO
        gettimeofday (&t2e, NULL);
#endif
        
        auto crosslinker = new Crosslinker();
        crosslinker = crosslinker->SelectCrosslinker(getXlSearchParameters()->getCrosslinkerType());
        if (getXlSearchParameters()->getCrosslinkerType() == CrosslinkerType::UserDefined)
        {
            crosslinker = GenerateUserDefinedCrosslinker(getXlSearchParameters());
        }
        
        MyFileManager *myFileManager = new MyFileManager(true);
        
        std::vector<CommonParameters *> fileSpecificCommonParams;
        for ( auto b : fileSettingsList ) {
            fileSpecificCommonParams.push_back(SetAllFileSpecificCommonParams(getCommonParameters(), b));
        }
        
        std::unordered_set<DigestionParams*> ListOfDigestionParams;
        for ( auto p : fileSpecificCommonParams ) {
            ListOfDigestionParams.emplace(p->getDigestionParams());
        }
        
        int completedFiles = 0;
        //std::any indexLock = std::any();
        //std::any psmLock = std::any();
        
        Status("Searching files...", taskId, getVerbose() );
        
        DigestionParams *dPar = getCommonParameters()->getDigestionParams();
        ProseCreatedWhileRunning->append("The following crosslink discovery were used: ");
        ProseCreatedWhileRunning->append("crosslinker name = " + crosslinker->getCrosslinkerName() + "; ");
        ProseCreatedWhileRunning->append("crosslinker type = " + StringHelper::toString(crosslinker->getCleavable()) + "; ");
        ProseCreatedWhileRunning->append("crosslinker mass = " + std::to_string(crosslinker->getTotalMass()) + "; ");
        ProseCreatedWhileRunning->append("crosslinker modification site(s) = " + crosslinker->getCrosslinkerModSites() + "; ");
        
        ProseCreatedWhileRunning->append("protease = " + dPar->getProtease()->ToString() + "; ");
        ProseCreatedWhileRunning->append("maximum missed cleavages = " + std::to_string(dPar->getMaxMissedCleavages()) + "; ");
        ProseCreatedWhileRunning->append("minimum peptide length = " + std::to_string(dPar->getMinPeptideLength()) + "; ");
        ProseCreatedWhileRunning->append(dPar->getMaxPeptideLength() == std::numeric_limits<int>::max() ?
                                         "maximum peptide length = unspecified; " : "maximum peptide length = " +
                                         std::to_string(dPar->getMaxPeptideLength()) + "; ");
        ProseCreatedWhileRunning->append("initiator methionine behavior = " +
                                         std::to_string(static_cast<int>(dPar->getInitiatorMethionineBehavior())) + "; ");
        ProseCreatedWhileRunning->append("max modification isoforms = " +
                                         std::to_string(dPar->getMaxModificationIsoforms()) + "; ");
        std::vector<std::string> vsv;
        for ( auto m : fixedModifications ) {
            vsv.push_back(m->getIdWithMotif());
        }
        std::string del = ", ";
        ProseCreatedWhileRunning->append("fixed modifications = " + StringHelper::join ( vsv, del) + "; ");

        vsv.clear();
        for ( auto m : variableModifications ) {
            vsv.push_back(m->getIdWithMotif());
        }
        ProseCreatedWhileRunning->append("variable modifications = " + StringHelper::join ( vsv, del) + "; ");
        
        ProseCreatedWhileRunning->append("parent mass tolerance(s) = " +
                                         std::to_string(getCommonParameters()->getPrecursorMassTolerance()->getValue()) + "; ");
        ProseCreatedWhileRunning->append("product mass tolerance = " +
                                         std::to_string(getCommonParameters()->getProductMassTolerance()->getValue()) + "; ");
        int c1 = 0;
        for ( auto p: proteinList ) {
            if (p->getIsContaminant()) {
                c1++;
            }
        }
        ProseCreatedWhileRunning->append("The combined search database contained " + std::to_string(proteinList.size()) +
                                         " total entries including " + std::to_string(c1) + " contaminant sequences. ");
        
        for (int spectraFileIndex = 0; spectraFileIndex < (int)currentRawFileList.size(); spectraFileIndex++)
        {
            auto origDataFile = currentRawFileList[spectraFileIndex];
            EngineLayer::CommonParameters *combinedParams = SetAllFileSpecificCommonParams(getCommonParameters(),
                                                                              fileSettingsList[spectraFileIndex]);
            
            auto thisId = std::vector<std::string> {taskId, "Individual Spectra Files", origDataFile};
            //NewCollection(FileSystem::getFileName(origDataFile), thisId);
            
            Status("Loading spectra file...", thisId,  getVerbose() );
#ifdef TIMING_INFO
            gettimeofday (&t3, NULL);
#endif
            MsDataFile *myMsDataFile = myFileManager->LoadFile(origDataFile,
                                                               std::make_optional(combinedParams->getTopNpeaks()),
                                                               std::make_optional(combinedParams->getMinRatio()),
                                                               combinedParams->getTrimMs1Peaks(),
                                                               combinedParams->getTrimMsMsPeaks(), combinedParams);
#ifdef TIMING_INFO
            gettimeofday (&t3e, NULL);
            t3total += timediff(t3, t3e);
#endif
            
            Status("Getting ms2 scans...", thisId, getVerbose());

#ifdef TIMING_INFO
            gettimeofday (&t4, NULL);
#endif
            std::vector<Ms2ScanWithSpecificMass*> arrayOfMs2ScansSortedByMass= GetMs2Scans(myMsDataFile, origDataFile, combinedParams);            
            std::sort(arrayOfMs2ScansSortedByMass.begin(), arrayOfMs2ScansSortedByMass.end(), [&]
                      (Ms2ScanWithSpecificMass* l, Ms2ScanWithSpecificMass* r) {
                          return l->getPrecursorMass() < r->getPrecursorMass();
                      });
#ifdef TIMING_INFO
            gettimeofday (&t4e, NULL);
            t4total += timediff (t4, t4e);
#endif            
            std::vector<CrosslinkSpectralMatch*> newPsms(arrayOfMs2ScansSortedByMass.size());
            for (int currentPartition = 0; currentPartition < getCommonParameters()->getTotalPartitions(); currentPartition++)
            {
                std::vector<PeptideWithSetModifications*> peptideIndex;

                int start = currentPartition * proteinList.size() / combinedParams->getTotalPartitions();
                int count = ((currentPartition + 1) * proteinList.size() / combinedParams->getTotalPartitions()) -
                    (currentPartition * proteinList.size() /combinedParams->getTotalPartitions());
                std::vector<Protein*> proteinListSubset;
                for ( auto p=0; p<count; p++ ) {
                    proteinListSubset.push_back(proteinList[start+p]);
                }
                
                std::vector<std::string> vs = {taskId};
                Status("Getting fragment dictionary...", vs, getVerbose());

                std::vector<std::string> filenameList;
                for ( auto p: dbFilenameList )
                {
                    filenameList.push_back(p->getFilePath());
                }
                std::vector<std::string> vtaskId = {taskId};

#ifdef TIMING_INFO
                gettimeofday (&t5, NULL);
#endif
                auto indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, currentPartition,
                                                      UsefulProteomicsDatabases::DecoyType::Reverse, combinedParams, 30000.0, false,
                                                      filenameList, vtaskId);
                std::vector<std::vector<int>> fragmentIndex;
                std::vector<std::vector<int>> precursorIndex;

                auto allmods = GlobalVariables::getAllModsKnown();
                GenerateIndexes(indexEngine, dbFilenameList, peptideIndex, fragmentIndex, precursorIndex, proteinList,
                                allmods, taskId);

#ifdef TIMING_INFO
                gettimeofday (&t5e, NULL);
                t5total += timediff(t5, t5e );
#endif
                
                Status("Searching files...", taskId, getVerbose());

#ifdef TIMING_INFO
                gettimeofday (&t6, NULL);
#endif
                CrosslinkSearchEngine tempVar(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition,
                                              combinedParams, crosslinker, getXlSearchParameters()->getRestrictToTopNHits(),
                                              getXlSearchParameters()->getCrosslinkSearchTopNum(),
                                              getXlSearchParameters()->getXlQuench_H2O(),
                                              getXlSearchParameters()->getXlQuench_NH2(),
                                              getXlSearchParameters()->getXlQuench_Tris(), thisId);
                (&tempVar)->Run();
#ifdef TIMING_INFO
                gettimeofday (&t6e, NULL);
                t6total += timediff (t6, t6e );
#endif
                std::string s1 = "Done with search " + std::to_string(currentPartition + 1) + "/" +
                    std::to_string(getCommonParameters()->getTotalPartitions()) + "!";
                ProgressEventArgs tempVar2(100, s1, thisId);
                ReportProgress(&tempVar2, getVerbose() );
            }
            
            for ( auto p : newPsms ) {
                if ( p != nullptr ) {
                    allPsms.push_back(p);
                }
            }
            
            completedFiles++;
            std::vector<std::string> vs2 = {taskId, "Individual Spectra Files"};
            ProgressEventArgs tempVar3(completedFiles / currentRawFileList.size(), "Searching...", vs2); 
            ReportProgress(&tempVar3, getVerbose());
        }
        
        std::vector<std::string> vs3 = {taskId, "Individual Spectra Files"};
        ProgressEventArgs tempVar4(100, "Done with all searches!", vs3);
        ReportProgress(&tempVar4, getVerbose());

        std::sort(allPsms.begin(), allPsms.end(), [&] (CrosslinkSpectralMatch *l, CrosslinkSpectralMatch *r) {
                return l->getXLTotalScore() > r->getXLTotalScore();
            });

        std::vector<CrosslinkSpectralMatch *> allPsmsXL;
        for ( auto p : allPsms ) {
            if ( p->getCrossType() == PsmCrossType::Cross ){
                allPsmsXL.push_back(p);
            }
        }
        
        // inter-crosslinks; different proteins are linked
        std::vector<CrosslinkSpectralMatch *> interCsms;
        for ( auto p: allPsmsXL ) {
            if ( p->getProteinAccession() != p->getBetaPeptide()->getProteinAccession() ) {
                interCsms.push_back(p);
            }
        }

        for (auto item : interCsms)
        {
            item->setCrossType(PsmCrossType::Inter);
        }
        
        // intra-crosslinks; crosslinks within a protein
        std::vector<CrosslinkSpectralMatch *> intraCsms;
        for ( auto p: allPsmsXL ) {
            if ( p->getProteinAccession() == p->getBetaPeptide()->getProteinAccession() ) {
                intraCsms.push_back(p);
            }
        }
        
        for (auto item : intraCsms)
        {
            item->setCrossType(PsmCrossType::Intra);
        }
        
        // calculate FDR
#ifdef TIMING_INFO
        gettimeofday (&t7, NULL);
#endif
        DoCrosslinkFdrAnalysis(interCsms);
        DoCrosslinkFdrAnalysis(intraCsms);
        std::vector<std::string> sv1 = {taskId};
        SingleFDRAnalysis(allPsms, sv1 );
#ifdef TIMING_INFO
        gettimeofday (&t7e, NULL);
#endif

        // calculate protein crosslink residue numbers
#ifdef TIMING_INFO
        gettimeofday (&t8, NULL);
#endif
        for (auto csm : allPsmsXL)
        {
            // alpha peptide crosslink residue in the protein
            csm->setXlProteinPos(csm->getOneBasedStartResidueInProtein().value() + csm->getLinkPositions()[0] - 1);
            
            // beta crosslink residue in protein
            csm->getBetaPeptide()->setXlProteinPos(csm->getBetaPeptide()->getOneBasedStartResidueInProtein().value() +
                                                   csm->getBetaPeptide()->getLinkPositions()[0] - 1);
        }
#ifdef TIMING_INFO
        gettimeofday (&t8e, NULL);
#endif        

        // write interlink CSMs
#ifdef TIMING_INFO
        gettimeofday (&t9, NULL);
#endif        
        if (!interCsms.empty())
        {
            std::string file = OutputFolder + "/XL_Interlinks.tsv";
            WritePsmCrossToTsv(interCsms, file, 2);
            std::vector<std::string> vs2 = {taskId};
            FinishedWritingFile(file, vs2, getVerbose());
        }

        int interCsms_size=0;
        for ( auto p: interCsms ) {
            if ( p->getFdrInfo()->getQValue() <= 0.01 && !p->getIsDecoy() && !p->getBetaPeptide()->getIsDecoy() ) {
                interCsms_size++;
            }
        }
        myTaskResults->AddNiceText("Target inter-crosslinks within 1% FDR: " + std::to_string(interCsms_size));
        
        if (getXlSearchParameters()->getWriteOutputForPercolator())
        {
            std::vector<CrosslinkSpectralMatch *> interPsmsXLPercolator ;
            for ( auto p: interCsms ) {
                if (p->getScore() >= 2 && p->getBetaPeptide()->getScore() >= 2 ) {
                    interPsmsXLPercolator.push_back(p);
                }
            }
            std::sort(interPsmsXLPercolator.begin(), interPsmsXLPercolator.end(), [&]
                      (CrosslinkSpectralMatch *l , CrosslinkSpectralMatch *r ) {
                          return l->getScanNumber() < r->getScanNumber();
                      });
            
            std::vector<std::string> vs2a = {taskId};
            WriteCrosslinkToTxtForPercolator(interPsmsXLPercolator, OutputFolder, "XL_Interlinks_Percolator",
                                             crosslinker, vs2a);
        }
        
        // write intralink CSMs
        if (!intraCsms.empty())
        {
            std::string file = OutputFolder + "/XL_Intralinks.tsv";
            WritePsmCrossToTsv(intraCsms, file, 2);
            std::vector<std::string> vs3 = {taskId};
            FinishedWritingFile(file, vs3, getVerbose());
        }

        int intraCsms_size=0;
        for ( auto p: intraCsms ) {
            if ( p->getFdrInfo()->getQValue() <= 0.01 && !p->getIsDecoy() && !p->getBetaPeptide()->getIsDecoy() ) {
                intraCsms_size++;
            }
        }
        myTaskResults->AddNiceText("Target intra-crosslinks within 1% FDR: " +std::to_string(intraCsms_size));
        
        if (getXlSearchParameters()->getWriteOutputForPercolator())
        {
            std::vector<CrosslinkSpectralMatch *> intraPsmsXLPercolator ;
            for ( auto p: intraCsms ) {
                if (p->getScore() >= 2 && p->getBetaPeptide()->getScore() >= 2 ) {
                    intraPsmsXLPercolator.push_back(p);
                }
            }
            std::sort(intraPsmsXLPercolator.begin(), intraPsmsXLPercolator.end(), [&]
                      (CrosslinkSpectralMatch *l , CrosslinkSpectralMatch *r ) {
                          return l->getScanNumber() < r->getScanNumber();
                      });
            
            std::vector<std::string> vs3a = {taskId};
            WriteCrosslinkToTxtForPercolator(intraPsmsXLPercolator, OutputFolder, "XL_Intralinks_Percolator",
                                             crosslinker, vs3a);
        }
        
        // write single peptides
        std::vector<CrosslinkSpectralMatch *> singlePsms;
        for ( auto p: allPsms ) {
            if ( p->getCrossType() == PsmCrossType::Single )  {
                singlePsms.push_back(p);
            }
        }

        if (!singlePsms.empty())
        {
            std::string writtenFileSingle = OutputFolder + "/SinglePeptides" + ".tsv";
            WritePsmCrossToTsv(singlePsms, writtenFileSingle, 1);
            std::vector<std::string> vs4 = {taskId};
            FinishedWritingFile(writtenFileSingle, vs4, getVerbose());
        }

        int singlePsms_size=0;
        for ( auto p: singlePsms ) {
            if (p->getFdrInfo()->getQValue() <= 0.01 && !p->getIsDecoy()) {
                singlePsms_size++;
            }
        }
        myTaskResults->AddNiceText("Target single peptides within 1% FDR: " + std::to_string(singlePsms_size));
        
        // write loops
        std::vector<CrosslinkSpectralMatch *> loopPsms;
        for ( auto p : allPsms ) {
            if ( p->getCrossType() == PsmCrossType::Loop ) {
                loopPsms.push_back(p);
            }
        }

        if (!loopPsms.empty())
        {
            std::string writtenFileLoop = OutputFolder + "/Looplinks" + ".tsv";
            WritePsmCrossToTsv(loopPsms, writtenFileLoop, 1);
            std::vector<std::string> vs4a = {taskId};
            FinishedWritingFile(writtenFileLoop, vs4a, getVerbose());
        }

        int loopPsms_size =0;
        for ( auto p : loopPsms ) {
            if ( p->getFdrInfo()->getQValue() <= 0.01 && !p->getIsDecoy())  {
                loopPsms_size++;
            }
        }
        myTaskResults->AddNiceText("Target loop-linked peptides within 1% FDR: " + std::to_string(loopPsms_size));
        
        // write deadends
        std::vector<CrosslinkSpectralMatch *> deadendPsms;
        for ( auto p : allPsms ) {
            if ( p->getCrossType() == PsmCrossType::DeadEnd ||
                 p->getCrossType() == PsmCrossType::DeadEndH2O ||
                 p->getCrossType() == PsmCrossType::DeadEndNH2 ||
                 p->getCrossType() == PsmCrossType::DeadEndTris) {
                deadendPsms.push_back(p);
            }
        }

        if (!deadendPsms.empty())
        {
            std::string writtenFileDeadend = OutputFolder + "/Deadends" + ".tsv";
            WritePsmCrossToTsv(deadendPsms, writtenFileDeadend, 1);
            std::vector<std::string> vs5a = {taskId};
            FinishedWritingFile(writtenFileDeadend, vs5a, getVerbose());
        }

        int deadendPsms_size =0;
        for ( auto p : deadendPsms ) {
            if ( p->getFdrInfo()->getQValue() <= 0.01 && !p->getIsDecoy())  {
                deadendPsms_size++;
            }
        }
        myTaskResults->AddNiceText("Target deadend peptides within 1% FDR: " + std::to_string(deadendPsms_size));
        
        // write pepXML
        if (getXlSearchParameters()->getWritePepXml())
        {
            std::vector<CrosslinkSpectralMatch*> writeToXml;

            for ( auto p: intraCsms ) {
                if ( !p->getIsDecoy() && !p->getBetaPeptide()->getIsDecoy() && p->getFdrInfo()->getQValue() <= 0.05) {
                    writeToXml.push_back(p);
                }
            }
            for ( auto p: interCsms ) {
                if ( !p->getIsDecoy() && !p->getBetaPeptide()->getIsDecoy() && p->getFdrInfo()->getQValue() <= 0.05 ) {
                    writeToXml.push_back(p);
                }
            }
            for ( auto p: singlePsms ) {
                if (!p->getIsDecoy() && p->getFdrInfo()->getQValue() <= 0.05) {
                    writeToXml.push_back(p);
                }
            }
            for ( auto p: loopPsms ) {
                if ( !p->getIsDecoy() && p->getFdrInfo()->getQValue() <= 0.05 ) {
                    writeToXml.push_back(p);
                }
            }
            for ( auto p: deadendPsms ) {
                if ( !p->getIsDecoy() && p->getFdrInfo()->getQValue() <= 0.05 ) {
                    writeToXml.push_back(p);
                }
            }

            std::sort(writeToXml.begin(), writeToXml.end(), [&] (CrosslinkSpectralMatch *l, CrosslinkSpectralMatch *r ) {
                    return l->getScanNumber() < r->getScanNumber();
                });

            for (auto fullFilePath : currentRawFileList)
            {
                std::string fileNameNoExtension = fullFilePath.substr(0, fullFilePath.find_last_of("."));

                std::vector<std::string> vs6 = {taskId};
                std::vector<CrosslinkSpectralMatch*> tmpwriteToXml;
                for ( auto p: writeToXml ) {
                    if ( p->getFullFilePath() == fullFilePath )  {
                        tmpwriteToXml.push_back(p);
                    }
                }
                WritePepXML_xl(tmpwriteToXml, proteinList, dbFilenameList[0]->getFilePath(), variableModifications,
                               fixedModifications, localizeableModificationTypes, OutputFolder, fileNameNoExtension,
                               vs6);
            }
        }
#ifdef TIMING_INFO
        gettimeofday (&t9e, NULL);

        std::cout << "Load Modifications : " << timediff(t1, t1e ) << " sec \n";
        std::cout << "Load Proteins : " << timediff(t2, t2e ) << " sec \n";
        std::cout << "Load Files : " << t3total << " sec \n";
        std::cout << "GetMs2Scans : " << t4total << " sec \n";
        std::cout << "GenerateIndixes : " << t5total << " sec \n";
        std::cout << "CrosslinkSearch : " << t6total << " sec \n";
        std::cout << "FdrAnalysis : " << timediff(t7, t7e) << " sec \n";
        std::cout << "Calculate Residue numbers : " << timediff(t8, t8e) << " sec \n";
        std::cout << "Write results : " << timediff(t9, t9e) << " sec \n";
#endif        
        delete myFileManager;
        //C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was
        //passed to a method or constructor. Handle memory management manually.
        return myTaskResults;
    }

    void XLSearchTask::SingleFDRAnalysis(std::vector<CrosslinkSpectralMatch*> &items, std::vector<std::string> &taskIds)
    {
        // calculate single PSM FDR

        std::vector<PeptideSpectralMatch*> psms;
        for ( auto p: items ) {
            if ( p->getCrossType() == PsmCrossType::Single ) {
                psms.push_back(dynamic_cast<PeptideSpectralMatch*>(p));
            }
        }
            
        FdrAnalysisEngine tempVar(psms, 0, getCommonParameters(), taskIds);
        (&tempVar)->Run();
        
        // calculate loop PSM FDR
        psms.clear();
        for ( auto p: items ) {
            if ( p->getCrossType() == PsmCrossType::Loop ) {
                psms.push_back(dynamic_cast<PeptideSpectralMatch*>(p));
            }
        }
        
        FdrAnalysisEngine tempVar2(psms, 0, getCommonParameters(), taskIds);
        (&tempVar2)->Run();

        // calculate deadend FDR
        psms.clear();
        for ( auto p: items ) {
            if ( p->getCrossType() == PsmCrossType::DeadEnd || p->getCrossType() == PsmCrossType::DeadEndH2O ||
                 p->getCrossType() == PsmCrossType::DeadEndNH2 || p->getCrossType() == PsmCrossType::DeadEndTris ) {
                psms.push_back(dynamic_cast<PeptideSpectralMatch*>(p));
            }
        }
        
        FdrAnalysisEngine tempVar3(psms, 0, getCommonParameters(), taskIds);
        (&tempVar3)->Run();
    }

    void XLSearchTask::DoCrosslinkFdrAnalysis(std::vector<CrosslinkSpectralMatch*> &csms)
    {
        int cumulativeTarget = 0;
        int cumulativeDecoy = 0;
        
        for (int i = 0; i < (int)csms.size(); i++)
        {
            auto csm = csms[i];
            if (csm->getIsDecoy() || csm->getBetaPeptide()->getIsDecoy())
            {
                cumulativeDecoy++;
            }
            else
            {
                cumulativeTarget++;
            }
            
            double qValue = std::min((double)1.0, static_cast<double>(cumulativeDecoy) / cumulativeTarget);
            csm->SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, 0, 0, 0, 0, 0, 0, false);
        }
        
        double qValueThreshold = 1.0;
        for (int i = csms.size() - 1; i >= 0; i--)
        {
            CrosslinkSpectralMatch *csm = csms[i];
            
            // threshold q-values
            if (csm->getFdrInfo()->getQValue() > qValueThreshold)
            {
                csm->getFdrInfo()->setQValue(qValueThreshold);
            }
            else if (csm->getFdrInfo()->getQValue() < qValueThreshold)
            {
                qValueThreshold = csm->getFdrInfo()->getQValue();
            }
        }
    }
    
    Crosslinker *XLSearchTask::GenerateUserDefinedCrosslinker(TaskLayer::XlSearchParameters *xlSearchParameters)
    {
        std::optional<double> tempVar = xlSearchParameters->getCrosslinkerTotalMass();
        std::optional<double> tempVar2 = xlSearchParameters->getCrosslinkerShortMass();
        std::optional<double> tempVar3 = xlSearchParameters->getCrosslinkerLongMass();
        std::optional<double> tempVar4 = xlSearchParameters->getCrosslinkerLoopMass();
        std::optional<double> tempVar5 = xlSearchParameters->getCrosslinkerDeadEndMassH2O();
        std::optional<double> tempVar6 = xlSearchParameters->getCrosslinkerDeadEndMassNH2();
        std::optional<double> tempVar7 = xlSearchParameters->getCrosslinkerDeadEndMassTris();
        auto crosslinker = new Crosslinker(xlSearchParameters->getCrosslinkerResidues(),
                                           xlSearchParameters->getCrosslinkerResidues2(),
                                           xlSearchParameters->getCrosslinkerName(),
                                           xlSearchParameters->getIsCleavable(),
                                           tempVar.has_value() ? tempVar.value() : NAN,
                                           tempVar2.has_value() ? tempVar2.value() : NAN,
                                           tempVar3.has_value() ? tempVar3.value() : NAN,
                                           tempVar4.has_value() ? tempVar4.value() : NAN,
                                           tempVar5.has_value() ? tempVar5.value() : NAN,
                                           tempVar6.has_value() ? tempVar6.value() : NAN,
                                           tempVar7.has_value() ? tempVar7.value() : NAN);

        return crosslinker;
    }
    
    void XLSearchTask::WritePsmCrossToTsv(std::vector<CrosslinkSpectralMatch*> &items, const std::string &filePath, int writeType)
    {
        if (items.empty())
        {
            return;
        }
        
        //StreamWriter output = StreamWriter(filePath);
        std::ofstream output (filePath);
        std::string header;
        switch (writeType)
        {
            case 1:
                header = CrosslinkSpectralMatch::GetTabSepHeaderSingle();
                break;
            case 2:
                header = CrosslinkSpectralMatch::GetTabSepHeaderCross();
                break;
            default:
                break;
        }
        output << header << std::endl;
        for (auto heh : items)
        {
            output << heh->ToString() << std::endl;
        }
        output.close();
    }
    
    void XLSearchTask::WriteCrosslinkToTxtForPercolator(std::vector<CrosslinkSpectralMatch*> &items, const std::string &outputFolder,
                                                        const std::string &fileName, Crosslinker *crosslinker,
                                                        std::vector<std::string> &nestedIds)
    {
        if (items.empty())
        {
            return;
        }
        std::string writtenFile = outputFolder + "/" + fileName + ".txt";

        {
            std::ofstream output(writtenFile);
            std::string s = "SpecId\tLabel\tScannr\tScore\tdScore\tNormRank\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum";
            std::string s2 = "\tPeptide\tProtein";
            output << s << s2 << std::endl;
            
            for (auto item : items)
            {
                if (item->getBaseSequence() != "" && item->getBetaPeptide()->getBaseSequence() != "" &&
                    item->getProteinAccession() != "" && item->getBetaPeptide()->getProteinAccession() != "")
                {
                    std::string x = "T";
                    int label = 1;
                    if (item->getIsDecoy() || item->getBetaPeptide()->getIsDecoy())
                    {
                        x = "D";
                        label = -1;
                    }
                    
                    output << x + "-" + std::to_string(item->getScanNumber()) + "-" +
                        std::to_string(item->getScanRetentionTime()) + "\t" +
                        std::to_string(label) + "\t" +
                        std::to_string(item->getScanNumber()) + "\t" +
                        std::to_string(item->getXLTotalScore()) + "\t" +
                        std::to_string(item->getDeltaScore()) + "\t" +
                        std::to_string(item->getXlRank()[0] + item->getXlRank()[1])+ "\t" +
                        std::to_string(item->getScanPrecursorCharge()) + "\t" +
                        std::to_string(item->getScanPrecursorMass()) + "\t" +
                        ((item->getPeptideMonisotopicMass().has_value() &&
                          item->getBetaPeptide()->getPeptideMonisotopicMass().has_value() ) ?
                         std::to_string((item->getScanPrecursorMass() - item->getBetaPeptide()->getPeptideMonisotopicMass().value() -
                           item->getPeptideMonisotopicMass().value() - crosslinker->getTotalMass()) /
                                        item->getScanPrecursorMass() * 1E6) : "---") + "\t" +
                        std::to_string(item->getBetaPeptide()->getBaseSequence().length()) + "\t"
                        + std::to_string(item->getBaseSequence().length()) + "\t" +
                        std::to_string(item->getBetaPeptide()->getBaseSequence().length() +
                                       item->getBaseSequence().length()) + "\t" + "-." +
                        item->getBaseSequence() + std::to_string(item->getLinkPositions().front()) + "--" +
                        item->getBetaPeptide()->getBaseSequence() +
                        std::to_string(item->getBetaPeptide()->getLinkPositions().front()) + ".-" + "\t" +
                        std::get<1>(item->getBestMatchingPeptides().front())->getProtein()->getAccession() + "(" +
                        std::to_string(item->getXlProteinPos()) + ")" + "\t" +
                        std::get<1>(item->getBetaPeptide()->getBestMatchingPeptides().front())->getProtein()->getAccession() +
                        "(" + std::to_string(item->getBetaPeptide()->getXlProteinPos()) + ")" << std::endl;
                }
            }
            output.close();
        }
        FinishedWritingFile(writtenFile, nestedIds, getVerbose());
    }

    void XLSearchTask::WritePepXML_xl(std::vector<CrosslinkSpectralMatch*> &items,
                                      std::vector<Protein*> &proteinList,
                                      const std::string &databasePath,
                                      std::vector<Modification*> &variableModifications,
                                      std::vector<Modification*> &fixedModifications,
                                      std::vector<std::string> &localizeableModificationTypes,
                                      const std::string &outputFolder,
                                      const std::string &fileName,
                                      std::vector<std::string> &nestedIds)
    {
        if (items.empty())
        {
            return;
        }
        
        //XmlSerializer *_indexedSerializer = new XmlSerializer(pepXML::Generated::msms_pipeline_analysis::typeid);
        //auto _pepxml = new pepXML::Generated::msms_pipeline_analysis();
        auto _pepxml = new pepXML::msms_pipeline_analysis();
        
        //_pepxml->date = DateTime::Now;
        //_pepxml->summary_xml = items[0]->getFullFilePath() + ".pep.XM";
        time_t timer;
        time(&timer);
        struct tm *tmi = localtime(&timer);

        short zone_hours=0;
        short zone_minutes=0;

        ::xml_schema::date_time *dt = new ::xml_schema::date_time(tmi->tm_year, tmi->tm_mon, tmi->tm_mday, tmi->tm_hour,
                                                                  tmi->tm_min, (double)tmi->tm_sec, zone_hours, zone_minutes);
        _pepxml->date(*dt);
        _pepxml->summary_xml(items[0]->getFullFilePath() + ".pep.XML");
        
        
        std::string proteaseC;
        std::string proteaseNC;
        for ( auto m : getCommonParameters()->getDigestionParams()->getProtease()->getDigestionMotifs() ) {
            proteaseC += m->InducingCleavage;
        }
        for ( auto m : getCommonParameters()->getDigestionParams()->getProtease()->getDigestionMotifs() ) {
            proteaseNC += m->PreventingCleavage;
        }

        auto  tempVar = new Crosslinker();
        Crosslinker *crosslinker = tempVar->SelectCrosslinker(getXlSearchParameters()->getCrosslinkerType());
        if (getXlSearchParameters()->getCrosslinkerType() == CrosslinkerType::UserDefined)
        {
            crosslinker = GenerateUserDefinedCrosslinker(getXlSearchParameters());
        }
        
        //std::string fileNameNoExtension = Path::GetFileNameWithoutExtension(items[0]->getFullFilePath());
        std::string temps = items[0]->getFullFilePath();
        std::string fileNameNoExtension = temps.substr(0, temps.find_last_of("."));
        std::string filePathNoExtension = temps.substr(0, temps.find_last_of("/"));

        std::string s1 = crosslinker->getCrosslinkerModSites();
        std::string s2 = crosslinker->getCrosslinkerModSites2();
        std::vector<char> vs1 (s1.begin(), s1.end() );
        std::vector<char> vs2 (s2.begin(), s2.end() );

        std::vector<char> modChars = vs1;
        for ( auto c : vs2 ) {
            bool found = false;
            for ( auto c2: modChars ) {
                if ( c == c2 ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                modChars.push_back(c);
            }
        }
        std::string modSites(modChars.begin(), modChars.end());
        
        //auto para = std::vector<pepXML::nameValueType*>();
        auto para = new pepXML::search_summary::parameter_sequence();
        {
            pepXML::nameValueType *tempVar2 = new pepXML::nameValueType();
            tempVar2->name("threads");
            tempVar2->value(std::to_string(getCommonParameters()->getMaxThreadsToUsePerFile()));
            para->push_back(*tempVar2);
            delete tempVar2;
            
            pepXML::nameValueType *tempVar3 = new pepXML::nameValueType();
            tempVar3->name("database");
            tempVar3->value(databasePath);
            para->push_back(*tempVar3);
            delete tempVar3;
            
            pepXML::nameValueType *tempVar4 = new pepXML::nameValueType();
            tempVar4->name("MS_data_file");
            tempVar4->value(items[0]->getFullFilePath());
            para->push_back(*tempVar4);
            delete tempVar4;
            
            pepXML::nameValueType *tempVar5 = new pepXML::nameValueType();
            tempVar5->name("Cross-link precursor Mass Tolerance");
            tempVar5->value(std::to_string(getCommonParameters()->getPrecursorMassTolerance()->getValue()));
            para->push_back(*tempVar5);
            delete tempVar5;
            
            pepXML::nameValueType *tempVar6 = new pepXML::nameValueType();
            tempVar6->name("Cross-linker type");
            tempVar6->value(crosslinker->getCrosslinkerName());
            para->push_back(*tempVar6);
            delete tempVar6;

            pepXML::nameValueType *tempVar7 = new pepXML::nameValueType();
            tempVar7->name("Cross-linker mass");
            tempVar7->value(std::to_string(crosslinker->getTotalMass()));
            para->push_back(*tempVar7);
            delete tempVar7;
            
            pepXML::nameValueType *tempVar8 = new pepXML::nameValueType();
            tempVar8->name("Cross-linker cleavable");
            tempVar8->value(StringHelper::toString(crosslinker->getCleavable()));
            para->push_back(*tempVar8);
            delete tempVar8;
            
            pepXML::nameValueType *tempVar9 = new pepXML::nameValueType();
            tempVar9->name("Cross-linker cleavable long mass");
            tempVar9->value(std::to_string(crosslinker->getCleaveMassLong()));
            para->push_back(*tempVar9);
            delete tempVar9;
            
            pepXML::nameValueType *tempVar10 = new pepXML::nameValueType();
            tempVar10->name("Cross-linker cleavable short mass");
            tempVar10->value(std::to_string(crosslinker->getCleaveMassShort()));
            para->push_back(*tempVar10);
            delete tempVar10;
            
            pepXML::nameValueType *tempVar11 = new pepXML::nameValueType();
            tempVar11->name("Cross-linker xl site");
            tempVar11->value(modSites);
            para->push_back(*tempVar11);
            delete tempVar11;
            
            pepXML::nameValueType *tempVar12 = new pepXML::nameValueType();
            tempVar12->name("Generate decoy proteins");
            auto tempVar12a = getXlSearchParameters()->getDecoyType();
            tempVar12->value(DecoyTypeToString(tempVar12a));
            para->push_back(*tempVar12);
            delete tempVar12;
            
            pepXML::nameValueType *tempVar13 = new pepXML::nameValueType();
            tempVar13->name( "MaxMissed Cleavages"); 
            tempVar13->value(std::to_string(getCommonParameters()->getDigestionParams()->getMaxMissedCleavages()));
            para->push_back(*tempVar13);
            delete tempVar13;
            
            pepXML::nameValueType *tempVar14 = new pepXML::nameValueType();
            tempVar14->name( "Protease");
            tempVar14->value(getCommonParameters()->getDigestionParams()->getProtease()->getName());
            para->push_back(*tempVar14);
            delete tempVar14;
            
            pepXML::nameValueType *tempVar15 = new pepXML::nameValueType();
            tempVar15->name( "Initiator Methionine");
            auto tmpVar15a = getCommonParameters()->getDigestionParams()->getInitiatorMethionineBehavior();
            tempVar15->value(InitiatorMethionineBehaviorToString(tmpVar15a));
            para->push_back(*tempVar15);
            delete tempVar15;
            
            pepXML::nameValueType *tempVar16 = new pepXML::nameValueType();
            tempVar16->name( "Max Modification Isoforms");
            tempVar16->value(std::to_string(getCommonParameters()->getDigestionParams()->getMaxModificationIsoforms()));
            para->push_back(*tempVar16);
            delete tempVar16;
            
            pepXML::nameValueType *tempVar17 = new pepXML::nameValueType();
            tempVar17->name( "Min Peptide Len");
            tempVar17->value(std::to_string(getCommonParameters()->getDigestionParams()->getMinPeptideLength()));
            para->push_back(*tempVar17);
            delete tempVar17;
            
            pepXML::nameValueType *tempVar18 = new pepXML::nameValueType();
            tempVar18->name( "Max Peptide Len");
            tempVar18->value(std::to_string(getCommonParameters()->getDigestionParams()->getMaxPeptideLength()));
            para->push_back(*tempVar18);
            delete tempVar18;
            
            pepXML::nameValueType *tempVar19 = new pepXML::nameValueType();
            tempVar19->name( "Product Mass Tolerance");
            tempVar19->value(std::to_string(getCommonParameters()->getProductMassTolerance()->getValue()));
            para->push_back(*tempVar19);
            delete tempVar19;
            
            pepXML::nameValueType *tempVar20 = new pepXML::nameValueType();
            tempVar20->name( "Ions to search");
            std::vector<ProductType> tempVar20a = DissociationTypeCollection::ProductsFromDissociationType[getCommonParameters()->getDissociationType()];
            std::string tempVar20s="";
            for ( auto p: tempVar20a ) {
                tempVar20s += Proteomics::Fragmentation::ProductTypeToString(p) + ", ";
            }
            tempVar20->value(tempVar20s);
            para->push_back(*tempVar20);
            delete tempVar20;
            
            for (auto fixedMod : fixedModifications)
            {
                pepXML::nameValueType *tempVar21 = new pepXML::nameValueType();
                tempVar21->name( "Fixed Modifications: " + fixedMod->getIdWithMotif());
                tempVar21->value(std::to_string(fixedMod->getMonoisotopicMass().value()));
                para->push_back(*tempVar21);
                delete tempVar21;
            }
            for (auto variableMod : variableModifications)
            {
                pepXML::nameValueType *tempVar22 = new pepXML::nameValueType();
                tempVar22->name( "Variable Modifications: " + variableMod->getIdWithMotif());
                tempVar22->value(std::to_string(variableMod->getMonoisotopicMass().value()));
                para->push_back(*tempVar22);
                delete tempVar22;
            }
            
            pepXML::nameValueType *tempVar23 = new pepXML::nameValueType();
            tempVar23->name( "Localize All Modifications");
            tempVar23->value("true");
            para->push_back(*tempVar23);
            delete tempVar23;
        }
        
        pepXML::msms_run_summary *tempVar24 = new pepXML::msms_run_summary();
        tempVar24->base_name( filePathNoExtension);
        tempVar24->raw_data_type("raw");
        tempVar24->raw_data(".mzM");

        tempVar24->sample_enzyme( *(new pepXML::sample_enzyme()));
        tempVar24->sample_enzyme()->name(getCommonParameters()->getDigestionParams()->getProtease()->getName());

        pepXML::specificity *tempVar25 = new pepXML::specificity();
        //Added by Edgar
        tempVar25->sense("C");
        //End Added by Edgar
        tempVar25->cut(proteaseC);
        tempVar25->no_cut(proteaseNC);
        auto t25 = new pepXML::sample_enzyme::specificity_sequence();
        t25->push_back(*tempVar25);
        delete tempVar25;
        tempVar24->sample_enzyme()->specificity(*t25);
        delete t25;
           
        pepXML::search_summary *tempVar26 = new pepXML::search_summary();
        //Added by Edgar:
        tempVar26->search_engine(pepXML::engineType::SEQUEST);
        //END added by Edgar
        tempVar26->base_name( filePathNoExtension);
        tempVar26->search_engine_version(GlobalVariables::getMetaMorpheusVersion());
        tempVar26->precursor_mass_type(pepXML::massType::monoisotopic);
        tempVar26->fragment_mass_type(pepXML::massType::monoisotopic);
        tempVar26->search_id(1);
        tempVar26->search_database(*(new pepXML::search_database()));
        tempVar26->search_database()->local_path(databasePath);
        tempVar26->search_database()->type(pepXML::type::value::AA);
        tempVar26->enzymatic_search_constraint(*(new pepXML::enzymatic_search_constraint()));
        tempVar26->enzymatic_search_constraint()->enzyme(getCommonParameters()->getDigestionParams()->getProtease()->getName());
        tempVar26->enzymatic_search_constraint()->max_num_internal_cleavages(getCommonParameters()->getDigestionParams()->getMaxMissedCleavages());
        //Added by Edgar:
        tempVar26->enzymatic_search_constraint()->min_number_termini(1);
        //END added by Edgar
        tempVar26->parameter(*para);
        delete para;
        
        auto t26 = new pepXML::msms_run_summary::search_summary_sequence();
        t26->push_back(*tempVar26);
        delete tempVar26;
        
        tempVar24->search_summary(*t26);
        delete t26;
        
        auto tt26 = new pepXML::msms_pipeline_analysis::msms_run_summary_sequence();
        tt26->push_back(*tempVar24);
        delete tempVar24;
        
        _pepxml->msms_run_summary(*tt26);
        delete tt26;

        auto ttt26 = new pepXML::msms_run_summary::spectrum_query_sequence(items.size());
        _pepxml->msms_run_summary()[0].spectrum_query(*ttt26);
        delete ttt26;
        
        auto searchHits = std::vector<pepXML::search_hit*>();
        for (int i = 0; i < (int)items.size(); i++)
        {
            //auto mods = std::vector<pepXML::mod_aminoacid_mass*>();
            auto mods = new pepXML::modInfoDataType::mod_aminoacid_mass_sequence();
            PeptideWithSetModifications *alphaPeptide = std::get<1>(items[i]->getBestMatchingPeptides().front());
            
            for (auto modification : alphaPeptide->getAllModsOneIsNterminus() )
            {
                auto mod = new pepXML::mod_aminoacid_mass();
                mod->mass(std::get<1>(modification)->getMonoisotopicMass().value());

                mod->position(std::get<0>(modification) - 1);
                mods->push_back(*mod);
                delete mod;
            }
            
            if (items[i]->getCrossType() == PsmCrossType::Single)
            {
                auto searchHit = new pepXML::search_hit();
                searchHit->hit_rank(1);
                searchHit->peptide(alphaPeptide->getBaseSequence());

                std::stringstream ss;
                ss << alphaPeptide->getPreviousAminoAcid();
                searchHit->peptide_prev_aa(ss.str());
                ss.str("");
                ss << alphaPeptide->getNextAminoAcid();
                searchHit->peptide_next_aa(ss.str());
                
                searchHit->protein(alphaPeptide->getProtein()->getAccession());
                searchHit->num_tot_proteins(1);
                searchHit->calc_neutral_pep_mass(static_cast<float>(items[i]->getScanPrecursorMass()));
                searchHit->massdiff((items[i]->getScanPrecursorMass() - items[i]->getPeptideMonisotopicMass().value()));
                //searchHit->xlink_typeSpecified(true);
                searchHit->xlink_type1().set(pepXML::xlink_type::na);
                searchHit->modification_info(*new pepXML::modInfoDataType());
                searchHit->modification_info()->mod_aminoacid_mass(*mods);

                pepXML::nameValueType *tempVar27 = new pepXML::nameValueType();
                tempVar27->name( "xlTotalScore");
                tempVar27->value(std::to_string(items[i]->getXLTotalScore()));

                pepXML::nameValueType *tempVar28 = new pepXML::nameValueType();
                tempVar28->name( "Qvalue");
                tempVar28->value(std::to_string(items[i]->getFdrInfo()->getQValue()));

                auto t28 = new pepXML::search_hit::search_score_sequence();
                t28->push_back(*tempVar27);
                delete tempVar27;
                t28->push_back(*tempVar28);
                delete tempVar28;
                searchHit->search_score(*t28);
                delete t28;
                searchHits.push_back(searchHit);
                
                //C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since
                //searchHit was passed to a method or constructor. Handle memory management manually.
            }
            else if (items[i]->getCrossType() == PsmCrossType::DeadEnd ||
                     items[i]->getCrossType() == PsmCrossType::DeadEndH2O ||
                     items[i]->getCrossType() == PsmCrossType::DeadEndNH2 ||
                     items[i]->getCrossType() == PsmCrossType::DeadEndTris)
            {
                double crosslinkerDeadEndMass = 0;
                switch (items[i]->getCrossType())
                {
                    case PsmCrossType::DeadEndNH2:
                        crosslinkerDeadEndMass = crosslinker->getDeadendMassNH2();
                        break;
                        
                    case PsmCrossType::DeadEndTris:
                        crosslinkerDeadEndMass = crosslinker->getDeadendMassTris();
                        break;
                        
                    default:
                        crosslinkerDeadEndMass = crosslinker->getDeadendMassH2O();
                        break;
                }
                auto mod = new pepXML::mod_aminoacid_mass();
                mod->mass(crosslinkerDeadEndMass );

                mod->position(items[i]->getLinkPositions().front());
                mods->push_back(*mod);
                delete mod;
                               
                auto searchHit = new pepXML::search_hit();
                searchHit->hit_rank(1);
                searchHit->peptide(alphaPeptide->getBaseSequence());

                std::stringstream ss;
                ss << alphaPeptide->getPreviousAminoAcid();
                searchHit->peptide_prev_aa(ss.str());
                ss.str("");
                ss << alphaPeptide->getNextAminoAcid();
                searchHit->peptide_next_aa(ss.str());

                searchHit->protein(alphaPeptide->getProtein()->getAccession());
                searchHit->num_tot_proteins(1);
                searchHit->calc_neutral_pep_mass(static_cast<float>(items[i]->getScanPrecursorMass()));
                searchHit->massdiff((items[i]->getScanPrecursorMass() - items[i]->getPeptideMonisotopicMass().value() - crosslinkerDeadEndMass));
                //searchHit->xlink_typeSpecified = true;
                searchHit->xlink_type1().set(pepXML::xlink_type::na);
                searchHit->modification_info(*(new pepXML::modInfoDataType()));
                searchHit->modification_info()->mod_aminoacid_mass(*mods);

                pepXML::nameValueType *tempVar29 = new pepXML::nameValueType();
                tempVar29->name( "xlTotalScore");
                tempVar29->value(std::to_string(items[i]->getXLTotalScore()));
                
                pepXML::nameValueType *tempVar30 = new pepXML::nameValueType();
                tempVar30->name( "Qvalue"); 
                tempVar30->value(std::to_string(items[i]->getFdrInfo()->getQValue()));

                auto t30 = new pepXML::search_hit::search_score_sequence();
                t30->push_back(*tempVar29);
                delete tempVar29;
                
                t30->push_back(*tempVar30);
                delete tempVar30;
                
                searchHit->search_score(*t30);
                delete t30;
                
                searchHits.push_back(searchHit);
                //C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since
                //searchHit was passed to a method or constructor. Handle memory management manually.
            }
            else if (items[i]->getCrossType() == PsmCrossType::Inter ||
                     items[i]->getCrossType() == PsmCrossType::Intra ||
                     items[i]->getCrossType() == PsmCrossType::Cross)
            {
                auto betaPeptide = std::get<1>(items[i]->getBetaPeptide()->getBestMatchingPeptides().front());
                auto modsBeta =  new pepXML::modInfoDataType::mod_aminoacid_mass_sequence();
                
                for (auto mod : betaPeptide->getAllModsOneIsNterminus())
                {
                    auto modBeta = new pepXML::mod_aminoacid_mass();
                    modBeta->mass(std::get<1>(mod)->getMonoisotopicMass().value());

                    modBeta->position(std::get<0>(mod) - 1);
                    modsBeta->push_back(*modBeta);
                    delete modBeta;
                }
                
                auto alpha = new pepXML::linked_peptide();
                alpha->peptide(alphaPeptide->getBaseSequence());
                std::stringstream ss;
                ss << alphaPeptide->getPreviousAminoAcid();
                alpha->peptide_prev_aa(ss.str());
                ss.str("");
                ss << alphaPeptide->getNextAminoAcid();
                alpha->peptide_next_aa(ss.str() );

                alpha->protein(alphaPeptide->getProtein()->getAccession());
                alpha->num_tot_proteins(1);
                alpha->calc_neutral_pep_mass(static_cast<float>(items[i]->getPeptideMonisotopicMass().value()));
                alpha->complement_mass(static_cast<float>(items[i]->getScanPrecursorMass() - alphaPeptide->getMonoisotopicMass()));
                alpha->designation("alpha");
                alpha->modification_info(*new pepXML::modInfoDataType());
                alpha->modification_info()->mod_aminoacid_mass(*mods);

                pepXML::nameValueType *tempVar31 = new pepXML::nameValueType();
                tempVar31->name("xlscore");
                tempVar31->value(std::to_string(items[i]->getXLTotalScore()));

                pepXML::nameValueType *tempVar32 = new pepXML::nameValueType();
                tempVar32->name("link");
                tempVar32->value(std::to_string(items[i]->getLinkPositions().front()));
                auto talpha = new pepXML::linked_peptide::xlink_score_sequence();
                talpha->push_back(*tempVar31);
                delete tempVar31;
                
                talpha->push_back(*tempVar32);
                delete tempVar32;
                
                alpha->xlink_score(*talpha);
                delete talpha;

                auto beta = new pepXML::linked_peptide();
                beta->peptide(betaPeptide->getBaseSequence());
                ss.str("");
                ss << betaPeptide->getPreviousAminoAcid();
                beta->peptide_prev_aa(ss.str());
                ss.str("");
                ss << betaPeptide->getNextAminoAcid();
                beta->peptide_next_aa(ss.str());
                beta->protein(betaPeptide->getProtein()->getAccession());
                beta->num_tot_proteins(1);
                beta->calc_neutral_pep_mass(static_cast<float>(betaPeptide->getMonoisotopicMass()));
                beta->complement_mass(static_cast<float>(items[i]->getScanPrecursorMass() - betaPeptide->getMonoisotopicMass()));
                beta->designation("beta");
                beta->modification_info(*(new pepXML::modInfoDataType()));
                beta->modification_info()->mod_aminoacid_mass(*modsBeta);

                pepXML::nameValueType *tempVar33 = new pepXML::nameValueType();
                tempVar33->name( "xlscore");
                tempVar33->value(std::to_string(items[i]->getBetaPeptide()->getScore()));

                pepXML::nameValueType *tempVar34 = new pepXML::nameValueType();
                tempVar34->name( "link");
                tempVar34->value(std::to_string(items[i]->getBetaPeptide()->getLinkPositions().front()));
                auto tbeta = new pepXML::linked_peptide::xlink_score_sequence();
                tbeta->push_back(*tempVar33);
                delete tempVar33;
                
                tbeta->push_back(*tempVar34);
                delete tempVar34;
                
                beta->xlink_score(*tbeta);
                delete tbeta;
                
                //auto cross = {alpha, beta};
                auto cross = new pepXML::xlink::linked_peptide_sequence();
                cross->push_back(*alpha);
                delete alpha;
                cross->push_back(*beta);
                delete beta;
                
                auto searchHit = new pepXML::search_hit();
                searchHit->hit_rank(1);
                searchHit->peptide("-");
                searchHit->peptide_prev_aa("-");
                searchHit->peptide_next_aa("-");
                searchHit->protein("-");
                searchHit->num_tot_proteins(1);
                searchHit->calc_neutral_pep_mass( static_cast<float>(items[i]->getScanPrecursorMass()));
                searchHit->massdiff((items[i]->getScanPrecursorMass() - betaPeptide->getMonoisotopicMass() - alphaPeptide->getMonoisotopicMass() - crosslinker->getTotalMass()));
                //searchHit->xlink_typeSpecified = true;
                searchHit->xlink_type1().set(pepXML::xlink_type::xl);
                searchHit->xlink(*new pepXML::xlink());
                searchHit->xlink()->identifier(crosslinker->getCrosslinkerName());
                searchHit->xlink()->mass(static_cast<float>(crosslinker->getTotalMass()));
                searchHit->xlink()->linked_peptide(*cross);

                pepXML::nameValueType *tempVar35 = new pepXML::nameValueType();
                tempVar35->name( "xlTotalScore");
                tempVar35->value(std::to_string(items[i]->getXLTotalScore()));

                pepXML::nameValueType *tempVar36 = new pepXML::nameValueType();
                tempVar36->name( "Qvalue");
                tempVar36->value(std::to_string(items[i]->getFdrInfo()->getQValue()));
                auto t36 = new pepXML::search_hit::search_score_sequence();
                t36->push_back(*tempVar35);
                delete tempVar35;
                
                t36->push_back(*tempVar36);
                delete tempVar36;
                
                searchHit->search_score(*t36);
                delete t36;
                
                searchHits.push_back(searchHit);
                //C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since
                //searchHit was passed to a method or constructor. Handle memory management manually.
            }
            else if (items[i]->getCrossType() == PsmCrossType::Loop)
            {
                auto thePeptide = new pepXML::linked_peptide();
                //Added by Edgar:
                thePeptide->peptide("-");
                //thePeptide->peptide_prev_aa("-");
                //thePeptide->peptide_next_aa("-");
                thePeptide->protein("-");
                thePeptide->num_tot_proteins(0);
                thePeptide->calc_neutral_pep_mass(static_cast<float>(0));
                thePeptide->complement_mass(static_cast<float>(0));
                //thePeptide->designation("test");
                //End Added by Edgar
                
                pepXML::nameValueType *tempVar37 = new pepXML::nameValueType();
                tempVar37->name( "link");
                tempVar37->value(std::to_string(items[i]->getLinkPositions().front()));

                pepXML::nameValueType *tempVar38 = new pepXML::nameValueType();
                tempVar38->name( "link");
                tempVar38->value(std::to_string(items[i]->getLinkPositions()[1]));

                auto t38 = new pepXML::linked_peptide::xlink_score_sequence();
                t38->push_back(*tempVar37);
                delete tempVar37;
                
                t38->push_back(*tempVar38);
                delete tempVar38;
                
                thePeptide->xlink_score(*t38);
                delete t38;

                //auto cross = {thePeptide};
                auto cross = new pepXML::xlink::linked_peptide_sequence();
                cross->push_back(*thePeptide);
                delete thePeptide;
                
                auto searchHit = new pepXML::search_hit();
                searchHit->hit_rank(1);
                searchHit->peptide(alphaPeptide->getBaseSequence());

                std::stringstream ss;
                ss << alphaPeptide->getPreviousAminoAcid();
                searchHit->peptide_prev_aa(ss.str());
                ss.str("");
                ss << alphaPeptide->getNextAminoAcid();
                searchHit->peptide_next_aa(ss.str() );
                searchHit->protein(alphaPeptide->getProtein()->getAccession());
                searchHit->num_tot_proteins(1);
                searchHit->calc_neutral_pep_mass(static_cast<float>(items[i]->getScanPrecursorMass()));
                searchHit->massdiff((items[i]->getScanPrecursorMass() - alphaPeptide->getMonoisotopicMass() - crosslinker->getLoopMass()));
                //searchHit->xlink_typeSpecified = true;
                searchHit->xlink_type1().set(pepXML::xlink_type::loop);
                searchHit->modification_info(*(new pepXML::modInfoDataType()));
                searchHit->modification_info()->mod_aminoacid_mass(*mods);
                searchHit->xlink() = *(new pepXML::xlink());
                searchHit->xlink()->identifier(crosslinker->getCrosslinkerName());
                searchHit->xlink()->mass(static_cast<float>(crosslinker->getTotalMass()));
                searchHit->xlink()->linked_peptide(*cross);
                delete cross;
                
                pepXML::nameValueType *tempVar39 = new pepXML::nameValueType();
                tempVar39->name( "xlTotalScore");
                tempVar39->value(std::to_string(items[i]->getXLTotalScore()));

                pepXML::nameValueType *tempVar40 = new pepXML::nameValueType();
                tempVar40->name( "Qvalue");
                tempVar40->value(std::to_string(items[i]->getFdrInfo()->getQValue()));
                auto t40 = new pepXML::search_hit::search_score_sequence();
                t40->push_back(*tempVar39);
                delete tempVar39;
                
                t40->push_back(*tempVar40);
                delete tempVar40;
                
                searchHit->search_score(*t40);
                delete t40;
                
                searchHits.push_back(searchHit);
                //C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since
                //searchHit was passed to a method or constructor. Handle memory management manually.
            }
        }
        
        for (int i = 0; i < (int)items.size(); i++)
        {
            auto tempVar41 = new pepXML::spectrum_query();

            tempVar41->spectrum(fileNameNoExtension + "." + std::to_string(items[i]->getScanNumber()));
            tempVar41->start_scan(static_cast<unsigned int>(items[i]->getScanNumber()));
            tempVar41->end_scan(static_cast<unsigned int>(items[i]->getScanNumber()));
            tempVar41->precursor_neutral_mass(static_cast<float>(items[i]->getScanPrecursorMass()));

            tempVar41->assumed_charge(items[i]->getScanPrecursorCharge());
            tempVar41->index(static_cast<unsigned int>(i + 1));
            tempVar41->retention_time_sec(static_cast<float>(items[i]->getScanRetentionTime() * 60));

            auto tempVar42 = new pepXML::search_result();
            auto t42 = new pepXML::search_result::search_hit_sequence();
            t42->push_back(*searchHits[i]);
            delete searchHits[i];
            
            tempVar42->search_hit(*t42);
            delete t42;
            
            auto t41 = new pepXML::spectrum_query::search_result_sequence();
            t41->push_back(*tempVar42);
            delete tempVar42;
            
            tempVar41->search_result(*t41);
            delete t41;
            
            _pepxml->msms_run_summary()[0].spectrum_query()[i] = *tempVar41;
            delete tempVar41;
        }
        
        // Serialize the object model to XML.
        //
        std::string outFileName = outputFolder + "/"+ fileName + ".pep.XM";
            
        xml_schema::namespace_infomap map;
        map[""].name = "";
        map[""].schema = "/home/gabriel/XLMS/mzlib-master/pepXML/pepXML_v120.xsd";

        try{
            std::ofstream ofs (outFileName);
            pepXML::msms_pipeline_analysis_ (ofs, *_pepxml, map);
            ofs.close();
        }

        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
        }
        FinishedWritingFile( outFileName, nestedIds, getVerbose());

        delete _pepxml;
    }
}
