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

using namespace EngineLayer;
using namespace EngineLayer::CrosslinkSearch;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace MzLibUtil;
using namespace EngineLayer::FdrAnalysis;
using namespace Proteomics::Fragmentation;

namespace TaskLayer
{

    XLSearchTask::XLSearchTask() : MetaMorpheusTask(MyTask::XLSearch)
    {
#ifdef ORIG
        EngineLayer::CommonParameters tempVar( , , = true, = true, = 3, = 12, = true, = false, = 1,
                                               3, , , , , , , , new PpmTolerance(10));
#endif
        std::string taskDescr = "";
        DissociationType dissType = DissociationType::HCD;
        int topNpeaks = 200;
        double minRatio = 0.01;
        bool trimMs1Peaks = false;
        bool trimMsMsPeaks = true;
        bool useDeltaScore = false;
        bool calculateEValue = false;
        
        EngineLayer::CommonParameters tempVar( taskDescr, dissType, true, true, 3, 12, true, false,
                                               1, 3, topNpeaks,
                                               minRatio, trimMs1Peaks, trimMsMsPeaks, useDeltaScore,
                                               calculateEValue, nullptr, new PpmTolerance(10));
        setCommonParameters(&tempVar);
        
        TaskLayer::XlSearchParameters tempVar2;
        setXlSearchParameters(&tempVar2);
    }
    
    TaskLayer::XlSearchParameters *XLSearchTask::getXlSearchParameters() const
    {
        return privateXlSearchParameters;
    }
    
    void XLSearchTask::setXlSearchParameters(TaskLayer::XlSearchParameters *value)
    {
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
        LoadModifications(taskId, variableModifications, fixedModifications, localizeableModificationTypes);
        
        // load proteins
        std::vector<Protein*> proteinList = LoadProteins(taskId, dbFilenameList, true,
                                                         getXlSearchParameters()->getDecoyType(),
                                                         localizeableModificationTypes,
                                                         getCommonParameters());
        
        auto crosslinker = new Crosslinker();
        crosslinker = crosslinker->SelectCrosslinker(getXlSearchParameters()->getCrosslinkerType());
        if (getXlSearchParameters()->getCrosslinkerType() == CrosslinkerType::UserDefined)
        {
            crosslinker = GenerateUserDefinedCrosslinker(getXlSearchParameters());
        }
        
        MyFileManager *myFileManager = new MyFileManager(true);
        
#ifdef ORIG
        auto fileSpecificCommonParams = fileSettingsList.Select([&] (std::any b) {
                SetAllFileSpecificCommonParams(getCommonParameters(), b);
            });
#endif
        std::vector<CommonParameters *> fileSpecificCommonParams;
        for ( auto b : fileSettingsList ) {
            fileSpecificCommonParams.push_back(SetAllFileSpecificCommonParams(getCommonParameters(), b));
        }
        
#ifdef ORIG
        std::unordered_set<DigestionParams*> ListOfDigestionParams = std::unordered_set<DigestionParams*>(fileSpecificCommonParams->Select([&] (std::any p) {
                    p::DigestionParams;
                }));
#endif
        std::unordered_set<DigestionParams*> ListOfDigestionParams;
        for ( auto p : fileSpecificCommonParams ) {
            ListOfDigestionParams.emplace(p->getDigestionParams());
        }
        
        int completedFiles = 0;
        //std::any indexLock = std::any();
        //std::any psmLock = std::any();
        
        Status("Searching files...", taskId);
        
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
#ifdef ORIG
        ProseCreatedWhileRunning->append("fixed modifications = " +
                                         std::string::Join(", ", fixedModifications->Select([&] (std::any m) {
                                                     m::IdWithMotif;
                                                 }) + "; "));
#endif
        std::vector<std::string> vsv;
        for ( auto m : fixedModifications ) {
            vsv.push_back(m->getIdWithMotif());
        }
        std::string del = ", ";
        ProseCreatedWhileRunning->append("fixed modifications = " + StringHelper::join ( vsv, del) + "; ");
#ifdef ORIG                
        ProseCreatedWhileRunning->append("variable modifications = " +
                                         std::string::Join(", ", variableModifications->Select([&] (std::any m) {
                                                     m::IdWithMotif;
                                                 })) + "; ");
#endif
        vsv.clear();
        for ( auto m : variableModifications ) {
            vsv.push_back(m->getIdWithMotif());
        }
        ProseCreatedWhileRunning->append("variable modifications = " + StringHelper::join ( vsv, del) + "; ");
        
        ProseCreatedWhileRunning->append("parent mass tolerance(s) = " +
                                         std::to_string(getCommonParameters()->getPrecursorMassTolerance()->getValue()) + "; ");
        ProseCreatedWhileRunning->append("product mass tolerance = " +
                                         std::to_string(getCommonParameters()->getProductMassTolerance()->getValue()) + "; ");
#ifdef ORIG
        ProseCreatedWhileRunning->append("The combined search database contained " + std::to_string(proteinList.size()) +
                                         " total entries including " + proteinList.Where([&] (std::any p) {
                                                 p::IsContaminant;
                                             })->Count() + " contaminant sequences. ");
#endif
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
            NewCollection(FileSystem::getFileName(origDataFile), thisId);
            
            Status("Loading spectra file...", thisId);
            MsDataFile *myMsDataFile = myFileManager->LoadFile(origDataFile,
                                                               std::make_optional(combinedParams->getTopNpeaks()),
                                                               std::make_optional(combinedParams->getMinRatio()),
                                                               combinedParams->getTrimMs1Peaks(),
                                                               combinedParams->getTrimMsMsPeaks(), combinedParams);
            
            Status("Getting ms2 scans...", thisId);
#ifdef ORIG
            std::vector<Ms2ScanWithSpecificMass*> arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy([&] (std::any b)  {
                    b::PrecursorMass;
                })->ToArray();
#endif
            std::vector<Ms2ScanWithSpecificMass*> arrayOfMs2ScansSortedByMass= GetMs2Scans(myMsDataFile, origDataFile, combinedParams);
            std::sort(arrayOfMs2ScansSortedByMass.begin(), arrayOfMs2ScansSortedByMass.end(), [&]
                      (Ms2ScanWithSpecificMass* l, Ms2ScanWithSpecificMass* r) {
                          return l->getPrecursorMass() < r->getPrecursorMass();
                      });
            
            std::vector<CrosslinkSpectralMatch*> newPsms(arrayOfMs2ScansSortedByMass.size());
            for (int currentPartition = 0; currentPartition < getCommonParameters()->getTotalPartitions(); currentPartition++)
            {
                std::vector<PeptideWithSetModifications*> peptideIndex;
                std::vector<Protein*> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.size() /
                                                                               combinedParams->getTotalPartitions(),
                                                                               ((currentPartition + 1) * proteinList.size() /
                                                                                combinedParams->getTotalPartitions()) -
                                                                               (currentPartition * proteinList.size() /
                                                                                combinedParams->getTotalPartitions()));
                
                std::vector<std::string> vs = {taskId};
                Status("Getting fragment dictionary...", vs);

                std::vector<std::string> filenameList;
                for ( auto p: dbFilenameList )
                {
                    filenameList.push_back(p->getFilePath());
                }
                std::vector<std::string> vtaskId = {taskId};
                auto indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, currentPartition,
                                                      UsefulProteomicsDatabases::DecoyType::Reverse, combinedParams, 30000.0, false,
                                                      filenameList, vtaskId);
                std::vector<std::vector<int>> fragmentIndex;
                std::vector<std::vector<int>> precursorIndex;

                auto allmods = GlobalVariables::getAllModsKnown();
                GenerateIndexes(indexEngine, dbFilenameList, peptideIndex, fragmentIndex, precursorIndex, proteinList,
                                allmods, taskId);
                
                Status("Searching files...", taskId);
                CrosslinkSearchEngine tempVar(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition,
                                              combinedParams, crosslinker, getXlSearchParameters()->getRestrictToTopNHits(),
                                              getXlSearchParameters()->getCrosslinkSearchTopNum(),
                                              getXlSearchParameters()->getXlQuench_H2O(),
                                              getXlSearchParameters()->getXlQuench_NH2(),
                                              getXlSearchParameters()->getXlQuench_Tris(), thisId);
                (&tempVar)->Run();

                std::string s1 = "Done with search " + std::to_string(currentPartition + 1) + "/" +
                    std::to_string(getCommonParameters()->getTotalPartitions()) + "!";
                ProgressEventArgs tempVar2(100, s1, thisId);
                ReportProgress(&tempVar2);
                
            }
            
#ifdef ORIG
            allPsms.AddRange(newPsms.Where([&] (std::any p) {
                        return p != nullptr;
                    }));
#endif
            for ( auto p : newPsms ) {
                if ( p != nullptr ) {
                    allPsms.push_back(p);
                }
            }
            
            completedFiles++;
            std::vector<std::string> vs2 = {taskId, "Individual Spectra Files"};
            ProgressEventArgs tempVar3(completedFiles / currentRawFileList.size(), "Searching...", vs2); 
            ReportProgress(&tempVar3);
        }
        
        std::vector<std::string> vs3 = {taskId, "Individual Spectra Files"};
        ProgressEventArgs tempVar4(100, "Done with all searches!", vs3);
        ReportProgress(&tempVar4);

#ifdef ORIG
        allPsms = allPsms.OrderByDescending([&] (std::any p)  {
                p::XLTotalScore;
            }).ToList();
#endif
        std::sort(allPsms.begin(), allPsms.end(), [&] (CrosslinkSpectralMatch *l, CrosslinkSpectralMatch *r) {
                return l->getXLTotalScore() > r->getXLTotalScore();
            });

#ifdef ORIG
        auto allPsmsXL = allPsms.Where([&] (std::any p)    {
                return p->CrossType == PsmCrossType::Cross;
            }).ToList();
#endif
        std::vector<CrosslinkSpectralMatch *> allPsmsXL;
        for ( auto p : allPsms ) {
            if ( p->getCrossType() == PsmCrossType::Cross ){
                allPsmsXL.push_back(p);
            }
        }
                
        
        // inter-crosslinks; different proteins are linked
#ifdef ORIG
        auto interCsms = allPsmsXL.Where([&] (std::any p) {
                !p::ProteinAccession->Equals(p::BetaPeptide::ProteinAccession);
            }).ToList();
#endif
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
#ifdef ORIG
        auto intraCsms = allPsmsXL.Where([&] (std::any p) {
                p::ProteinAccession->Equals(p::BetaPeptide::ProteinAccession);
            }).ToList();
#endif
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
        DoCrosslinkFdrAnalysis(interCsms);
        DoCrosslinkFdrAnalysis(intraCsms);
        std::vector<std::string> sv1 = {taskId};
        SingleFDRAnalysis(allPsms, sv1 );
        
        // calculate protein crosslink residue numbers
        for (auto csm : allPsmsXL)
        {
            // alpha peptide crosslink residue in the protein
            csm->setXlProteinPos(csm->getOneBasedStartResidueInProtein().value() + csm->getLinkPositions()[0] - 1);
            
            // beta crosslink residue in protein
            csm->getBetaPeptide()->setXlProteinPos(csm->getBetaPeptide()->getOneBasedStartResidueInProtein().value() +
                                                   csm->getBetaPeptide()->getLinkPositions()[0] - 1);
        }
        
        // write interlink CSMs
        if (!interCsms.empty())
        {
            std::string file = OutputFolder + "/XL_Interlinks.tsv";
            WritePsmCrossToTsv(interCsms, file, 2);
            std::vector<std::string> vs2 = {taskId};
            FinishedWritingFile(file, vs2);
        }

#ifdef ORIG
        MyTaskResults->AddNiceText("Target inter-crosslinks within 1% FDR: " + interCsms.size()([&] (std::any p)  {
                    return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy && !p::BetaPeptide::IsDecoy;
		}));
#endif
        int interCsms_size=0;
        for ( auto p: interCsms ) {
            if ( p->getFdrInfo()->getQValue() <= 0.01 && !p->getIsDecoy() && !p->getBetaPeptide()->getIsDecoy() ) {
                interCsms_size++;
            }
        }
        myTaskResults->AddNiceText("Target inter-crosslinks within 1% FDR: " + std::to_string(interCsms_size));
        
        if (getXlSearchParameters()->getWriteOutputForPercolator())
        {
#ifdef ORIG
            auto interPsmsXLPercolator = interCsms.Where([&] (std::any p) {
                    return p::Score >= 2 && p::BetaPeptide::Score >= 2;
                }).OrderBy([&] (std::any p)  {
                        p::ScanNumber;
                    }).ToList();
#endif
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
            FinishedWritingFile(file, vs3);
        }

#ifdef ORIG
        myTaskResults->AddNiceText("Target intra-crosslinks within 1% FDR: " + intraCsms.size()([&] (std::any p) {         			return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy && !p::BetaPeptide::IsDecoy;
		}));
#endif
        int intraCsms_size=0;
        for ( auto p: intraCsms ) {
            if ( p->getFdrInfo()->getQValue() <= 0.01 && !p->getIsDecoy() && !p->getBetaPeptide()->getIsDecoy() ) {
                intraCsms_size++;
            }
        }
        myTaskResults->AddNiceText("Target intra-crosslinks within 1% FDR: " +std::to_string(intraCsms_size));
        
        if (getXlSearchParameters()->getWriteOutputForPercolator())
        {
#ifdef ORIG
            auto intraPsmsXLPercolator = intraCsms.Where([&] (std::any p)  {
                    return p::Score >= 2 && p::BetaPeptide::Score >= 2;
                }).OrderBy([&] (std::any p)  {
                        p::ScanNumber;
                    }).ToList();
#endif
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
#ifdef ORIG
        auto singlePsms = allPsms.Where([&] (std::any p)  {
                return p->CrossType == PsmCrossType::Single;
            }).ToList();
#endif
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
            FinishedWritingFile(writtenFileSingle, vs4);
        }
#ifdef ORIG
        myTaskResults->AddNiceText("Target single peptides within 1% FDR: " + singlePsms.size()([&] (std::any p) {         			return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy;
		}));
#endif
        int singlePsms_size=0;
        for ( auto p: singlePsms ) {
            if (p->getFdrInfo()->getQValue() <= 0.01 && !p->getIsDecoy()) {
                singlePsms_size++;
            }
        }
        myTaskResults->AddNiceText("Target single peptides within 1% FDR: " + std::to_string(singlePsms_size));
        
        // write loops
#ifdef ORIG
        auto loopPsms = allPsms.Where([&] (std::any p) {
                return p->CrossType == PsmCrossType::Loop;
            }).ToList();
#endif
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
            FinishedWritingFile(writtenFileLoop, vs4a);
        }
#ifdef ORIG
        myTaskResults->AddNiceText("Target loop-linked peptides within 1% FDR: " + loopPsms.size()([&] (std::any p)  {
                    return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy;
		}));
#endif
        int loopPsms_size =0;
        for ( auto p : loopPsms ) {
            if ( p->getFdrInfo()->getQValue() <= 0.01 && !p->getIsDecoy())  {
                loopPsms_size++;
            }
        }
        myTaskResults->AddNiceText("Target loop-linked peptides within 1% FDR: " + std::to_string(loopPsms_size));
        
        // write deadends
#ifdef ORIG
        auto deadendPsms = allPsms.Where([&] (std::any p)  {
			return p->CrossType == PsmCrossType::DeadEnd || p->CrossType == PsmCrossType::DeadEndH2O ||
                        p->CrossType == PsmCrossType::DeadEndNH2 || p->CrossType == PsmCrossType::DeadEndTris;
            }).ToList();
#endif
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
            FinishedWritingFile(writtenFileDeadend, vs5a);
        }
#ifdef ORIG
        myTaskResults->AddNiceText("Target deadend peptides within 1% FDR: " + deadendPsms.size()([&] (std::any p)   {
                    return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy;
		}));
#endif
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

#ifdef ORIG
            writeToXml.AddRange(intraCsms.Where([&] (std::any p)		{
                        return !p::IsDecoy && !p::BetaPeptide::IsDecoy && p::FdrInfo::QValue <= 0.05;
                    }));
            writeToXml.AddRange(interCsms.Where([&] (std::any p)  {
                        return !p::IsDecoy && !p::BetaPeptide::IsDecoy && p::FdrInfo::QValue <= 0.05;
                    }));
            writeToXml.AddRange(singlePsms.Where([&] (std::any p){
                        return !p::IsDecoy && p::FdrInfo::QValue <= 0.05;
                    }));
            writeToXml.AddRange(loopPsms.Where([&] (std::any p) {
                        return !p::IsDecoy && p::FdrInfo::QValue <= 0.05;
                    }));
            writeToXml.AddRange(deadendPsms.Where([&] (std::any p) {
                        return !p::IsDecoy && p::FdrInfo::QValue <= 0.05;
                    }));
#endif
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

#ifdef ORIG
            writeToXml = writeToXml.OrderBy([&] (std::any p)   {
                    p::ScanNumber;
                }).ToList();
#endif
            std::sort(writeToXml.begin(), writeToXml.end(), [&] (CrosslinkSpectralMatch *l, CrosslinkSpectralMatch *r ) {
                    return l->getScanNumber() < r->getScanNumber();
                });

            for (auto fullFilePath : currentRawFileList)
            {
                std::string fileNameNoExtension = fullFilePath.substr(0, fullFilePath.find_last_of("."));

#ifdef ORIG
                WritePepXML_xl(writeToXml.Where([&] (std::any p) {
                            return p->FullFilePath == fullFilePath;
                        }).ToList(), proteinList, dbFilenameList[0]->getFilePath(), variableModifications,
                    fixedModifications, localizeableModificationTypes, OutputFolder, fileNameNoExtension,
                    std::vector<std::string> {taskId});
#endif
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

        delete myFileManager;
        //C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was
        //passed to a method or constructor. Handle memory management manually.
        return myTaskResults;
    }

    void XLSearchTask::SingleFDRAnalysis(std::vector<CrosslinkSpectralMatch*> &items, std::vector<std::string> &taskIds)
    {
        // calculate single PSM FDR

#ifdef ORIG
        std::vector<PeptideSpectralMatch*> psms = items.Where([&] (std::any p)	{
                return p->CrossType == PsmCrossType::Single;
            })->Select([&] (std::any p)	{
                    dynamic_cast<PeptideSpectralMatch*>(p);
		}).ToList();
#endif
        std::vector<PeptideSpectralMatch*> psms;
        for ( auto p: items ) {
            if ( p->getCrossType() == PsmCrossType::Single ) {
                psms.push_back(dynamic_cast<PeptideSpectralMatch*>(p));
            }
        }
            
        FdrAnalysisEngine tempVar(psms, 0, getCommonParameters(), taskIds);
        (&tempVar)->Run();
        
        // calculate loop PSM FDR
#ifdef ORIG
        psms = items.Where([&] (std::any p)   	{
                return p->CrossType == PsmCrossType::Loop;
            })->Select([&] (std::any p)           {
                    dynamic_cast<PeptideSpectralMatch*>(p);
		}).ToList();
#endif
        psms.clear();
        for ( auto p: items ) {
            if ( p->getCrossType() == PsmCrossType::Loop ) {
                psms.push_back(dynamic_cast<PeptideSpectralMatch*>(p));
            }
        }
        
        FdrAnalysisEngine tempVar2(psms, 0, getCommonParameters(), taskIds);
        (&tempVar2)->Run();

        // calculate deadend FDR
#ifdef ORIG
        psms = items.Where([&] (std::any p)	{
                return p->CrossType == PsmCrossType::DeadEnd || p->CrossType == PsmCrossType::DeadEndH2O ||
                p->CrossType == PsmCrossType::DeadEndNH2 || p->CrossType == PsmCrossType::DeadEndTris;
            })->Select([&] (std::any p)	{
                    dynamic_cast<PeptideSpectralMatch*>(p);
		}).ToList();
#endif
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
        
        StreamWriter output = StreamWriter(filePath);
        std::string header = "";
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
        output.WriteLine(header);
        for (auto heh : items)
        {
            output.WriteLine(heh->ToString());
        }
    }
    
    void XLSearchTask::WriteCrosslinkToTxtForPercolator(std::vector<CrosslinkSpectralMatch*> &items, const std::string &outputFolder, const std::string &fileName, Crosslinker *crosslinker, std::vector<std::string> &nestedIds)
    {
        if (items.empty())
        {
            return;
        }
        auto writtenFile = FileSystem::combine(outputFolder, fileName + ".txt");

        {
            StreamWriter output = StreamWriter(writtenFile);
            output.WriteLine(std::string("SpecId\tLabel\tScannr\tScore\tdScore\tNormRank\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum") + "\tPeptide\tProtein");
            for (auto item : items)
            {
                if (item->getBaseSequence() != "" && item->getBetaPeptide()->getBaseSequence() != "" && item->getProteinAccession() != "" && item->getBetaPeptide()->getProteinAccession() != "")
                {
                    std::string x = "T";
                    int label = 1;
                    if (item->getIsDecoy() || item->getBetaPeptide()->getIsDecoy())
                    {
                        x = "D";
                        label = -1;
                    }
                    
                    output.WriteLine(x + "-" + item->getScanNumber().ToString() + "-" +
                                     item->getScanRetentionTime().ToString() + "\t" +
                                     label.ToString() + "\t" +
                                     item->getScanNumber().ToString() + "\t" +
                                     item->getXLTotalScore().ToString() + "\t" +
                                     item->getDeltaScore().ToString() + "\t" +
                                     (item->getXlRank()[0] + item->getXlRank()[1]).ToString() + "\t" +
                                     item->getScanPrecursorCharge().ToString() + "\t" +
                                     item->getScanPrecursorMass().ToString() + "\t" +
                                     ((item->getPeptideMonisotopicMass().HasValue &&
                                       item->getBetaPeptide()->getPeptideMonisotopicMass().HasValue) ?
                                      ((item->getScanPrecursorMass() - item->getBetaPeptide()->getPeptideMonisotopicMass().Value -
                                        item->getPeptideMonisotopicMass().Value - crosslinker->getTotalMass()) /
                                       item->getScanPrecursorMass() * 1E6).ToString() : "---") + "\t" +
                                     item->getBetaPeptide()->getBaseSequence().length().ToString() + "\t"
                                     + item->getBaseSequence().length().ToString() + "\t" +
                                     (item->getBetaPeptide()->getBaseSequence().length() +
                                      item->getBaseSequence().length()).ToString() + "\t" + "-." +
                                     item->getBaseSequence() + item->getLinkPositions().front().ToString() + "--" +
                                     item->getBetaPeptide()->getBaseSequence() +
                                     item->getBetaPeptide()->getLinkPositions().front().ToString() + ".-" + "\t" +
                                     item->BestMatchingPeptides.First().Peptide.Protein.Accession.ToString() + "(" +
                                     item->getXlProteinPos().ToString() + ")" + "\t" +
                                     item->getBetaPeptide()->BestMatchingPeptides.First().Peptide.Protein.Accession.ToString() +
                                     "(" + item->getBetaPeptide()->getXlProteinPos().ToString() + ")");
                }
            }
        }
        FinishedWritingFile(writtenFile, nestedIds);
    }

    void XLSearchTask::WritePepXML_xl(std::vector<CrosslinkSpectralMatch*> &items, std::vector<Protein*> &proteinList, const std::string &databasePath, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<std::string> &localizeableModificationTypes, const std::string &outputFolder, const std::string &fileName, std::vector<std::string> &nestedIds)
    {
        if (!items.Any())
        {
            return;
        }
        
        XmlSerializer *_indexedSerializer = new XmlSerializer(pepXML::Generated::msms_pipeline_analysis::typeid);
        auto _pepxml = new pepXML::Generated::msms_pipeline_analysis();
        
        _pepxml->date = DateTime::Now;
        _pepxml->summary_xml = items[0]->getFullFilePath() + ".pep.XM";
        
        std::string proteaseC = "";
        std::string proteaseNC = "";
        for (auto x : getCommonParameters()->getDigestionParams()->Protease.DigestionMotifs->Select([&] (std::any m) {
                    m::InducingCleavage;
                }))
        {
            proteaseC += x;
        }
        for (auto x : getCommonParameters()->getDigestionParams()->Protease.DigestionMotifs->Select([&] (std::any m)       {
                    m::PreventingCleavage;
		}))
        {
            proteaseNC += x;
        }

        Crosslinker tempVar();
        Crosslinker *crosslinker = (&tempVar)->SelectCrosslinker(getXlSearchParameters()->getCrosslinkerType());
        if (getXlSearchParameters()->getCrosslinkerType() == CrosslinkerType::UserDefined)
        {
            crosslinker = GenerateUserDefinedCrosslinker(getXlSearchParameters());
        }
        
        std::string fileNameNoExtension = Path::GetFileNameWithoutExtension(items[0]->getFullFilePath());
        std::string filePathNoExtension = Path::ChangeExtension(items[0]->getFullFilePath(), "");

        std::string modSites = crosslinker->getCrosslinkerModSites().ToCharArray().Concat(crosslinker->getCrosslinkerModSites2().ToCharArray())->Distinct().ToString();

        auto para = std::vector<pepXML::Generated::nameValueType*>();
        {
            pepXML::Generated::nameValueType *tempVar2 = new pepXML::Generated::nameValueType();
            tempVar2->name = "threads";
            tempVar2->value = std::to_string(getCommonParameters()->getMaxThreadsToUsePerFile());
            para.push_back(tempVar2);
            pepXML::Generated::nameValueType *tempVar3 = new pepXML::Generated::nameValueType();
            tempVar3->name = "database";
            tempVar3->value = databasePath;
            para.push_back(tempVar3);
            pepXML::Generated::nameValueType *tempVar4 = new pepXML::Generated::nameValueType();
            tempVar4->name = "MS_data_file";
            tempVar4->value = items[0]->getFullFilePath();
            para.push_back(tempVar4);
            pepXML::Generated::nameValueType *tempVar5 = new pepXML::Generated::nameValueType();
            tempVar5->name = "Cross-link precursor Mass Tolerance";

            tempVar5->value = getCommonParameters()->getPrecursorMassTolerance()->ToString();
            para.push_back(tempVar5);
            pepXML::Generated::nameValueType *tempVar6 = new pepXML::Generated::nameValueType();
            tempVar6->name = "Cross-linker type";
            tempVar6->value = crosslinker->getCrosslinkerName();
            para.push_back(tempVar6);
            pepXML::Generated::nameValueType *tempVar7 = new pepXML::Generated::nameValueType();
            tempVar7->name = "Cross-linker mass";
            tempVar7->value = std::to_string(crosslinker->getTotalMass());
            para.push_back(tempVar7);
            pepXML::Generated::nameValueType *tempVar8 = new pepXML::Generated::nameValueType();
            tempVar8->name = "Cross-linker cleavable";
            tempVar8->value = StringHelper::toString(crosslinker->getCleavable());
            para.push_back(tempVar8);
            pepXML::Generated::nameValueType *tempVar9 = new pepXML::Generated::nameValueType();
            tempVar9->name = "Cross-linker cleavable long mass";
            tempVar9->value = std::to_string(crosslinker->getCleaveMassLong());
            para.push_back(tempVar9);
            pepXML::Generated::nameValueType *tempVar10 = new pepXML::Generated::nameValueType();
            tempVar10->name = "Cross-linker cleavable short mass";
            tempVar10->value = std::to_string(crosslinker->getCleaveMassShort());
            para.push_back(tempVar10);
            pepXML::Generated::nameValueType *tempVar11 = new pepXML::Generated::nameValueType();
            tempVar11->name = "Cross-linker xl site";
            tempVar11->value = modSites;
            para.push_back(tempVar11);
            
            pepXML::Generated::nameValueType *tempVar12 = new pepXML::Generated::nameValueType();
            tempVar12->name = "Generate decoy proteins";

            tempVar12->value = getXlSearchParameters()->getDecoyType()->ToString();
            para.push_back(tempVar12);
            pepXML::Generated::nameValueType *tempVar13 = new pepXML::Generated::nameValueType();
            tempVar13->name = "MaxMissed Cleavages";
            
            tempVar13->value = getCommonParameters()->getDigestionParams()->MaxMissedCleavages.ToString();
            para.push_back(tempVar13);
            pepXML::Generated::nameValueType *tempVar14 = new pepXML::Generated::nameValueType();
            tempVar14->name = "Protease";
            tempVar14->value = getCommonParameters()->getDigestionParams()->Protease->Name;
            para.push_back(tempVar14);
            pepXML::Generated::nameValueType *tempVar15 = new pepXML::Generated::nameValueType();
            tempVar15->name = "Initiator Methionine";

            tempVar15->value = getCommonParameters()->getDigestionParams()->InitiatorMethionineBehavior.ToString();
            para.push_back(tempVar15);
            pepXML::Generated::nameValueType *tempVar16 = new pepXML::Generated::nameValueType();
            tempVar16->name = "Max Modification Isoforms";

            tempVar16->value = getCommonParameters()->getDigestionParams()->MaxModificationIsoforms.ToString();
            para.push_back(tempVar16);
            pepXML::Generated::nameValueType *tempVar17 = new pepXML::Generated::nameValueType();
            tempVar17->name = "Min Peptide Len";

            tempVar17->value = getCommonParameters()->getDigestionParams()->MinPeptideLength.ToString();
            para.push_back(tempVar17);
            pepXML::Generated::nameValueType *tempVar18 = new pepXML::Generated::nameValueType();
            tempVar18->name = "Max Peptide Len";

            tempVar18->value = getCommonParameters()->getDigestionParams()->MaxPeptideLength.ToString();
            para.push_back(tempVar18);
            pepXML::Generated::nameValueType *tempVar19 = new pepXML::Generated::nameValueType();
            tempVar19->name = "Product Mass Tolerance";

            tempVar19->value = getCommonParameters()->getProductMassTolerance()->ToString();
            para.push_back(tempVar19);
            pepXML::Generated::nameValueType *tempVar20 = new pepXML::Generated::nameValueType();
            tempVar20->name = "Ions to search";
            tempVar20->value = std::string::Join(", ", DissociationTypeCollection::ProductsFromDissociationType[getCommonParameters()->getDissociationType()]);
            para.push_back(tempVar20);
            
            for (auto fixedMod : fixedModifications)
            {
                pepXML::Generated::nameValueType *tempVar21 = new pepXML::Generated::nameValueType();
                tempVar21->name = "Fixed Modifications: " + fixedMod->IdWithMotif;
                tempVar21->value = fixedMod->MonoisotopicMass.ToString();
                para.push_back(tempVar21);
            }
            for (auto variableMod : variableModifications)
            {
                pepXML::Generated::nameValueType *tempVar22 = new pepXML::Generated::nameValueType();
                tempVar22->name = "Variable Modifications: " + variableMod->IdWithMotif;
                tempVar22->value = variableMod->MonoisotopicMass.ToString();
                para.push_back(tempVar22);
            }
            
            pepXML::Generated::nameValueType *tempVar23 = new pepXML::Generated::nameValueType();
            tempVar23->name = "Localize All Modifications";
            tempVar23->value = "true";
            para.push_back(tempVar23);
        }
        
        pepXML::Generated::msms_pipeline_analysisMsms_run_summary *tempVar24 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summary();
        tempVar24->base_name = filePathNoExtension;
        tempVar24->raw_data_type = "raw";
        tempVar24->raw_data = ".mzM";
        tempVar24->sample_enzyme = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySample_enzyme();
        tempVar24->sample_enzyme->name = getCommonParameters()->getDigestionParams()->Protease->Name;
        pepXML::Generated::msms_pipeline_analysisMsms_run_summarySample_enzymeSpecificity *tempVar25 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySample_enzymeSpecificity();
        tempVar25->cut = proteaseC;
        tempVar25->no_cut = proteaseNC;
        tempVar24->sample_enzyme->specificity = {tempVar25};
        pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summary *tempVar26 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summary();
        tempVar26->base_name = filePathNoExtension;
        tempVar26->search_engine_version = GlobalVariables::getMetaMorpheusVersion();
        tempVar26->precursor_mass_type = pepXML::Generated::massType::monoisotopic;
        tempVar26->fragment_mass_type = pepXML::Generated::massType::monoisotopic;
        tempVar26->search_id = 1;
        tempVar26->search_database = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summarySearch_database();
        tempVar26->search_database->local_path = databasePath;
        tempVar26->search_database->type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summarySearch_databaseType::AA;
        tempVar26->enzymatic_search_constraint = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summaryEnzymatic_search_constraint();
        tempVar26->enzymatic_search_constraint->enzyme = getCommonParameters()->getDigestionParams()->Protease->Name;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
        tempVar26->enzymatic_search_constraint->max_num_internal_cleavages = getCommonParameters()->getDigestionParams()->MaxMissedCleavages.ToString();
        tempVar26->parameter = para.ToArray();
        tempVar24->search_summary = {tempVar26};
        _pepxml->msms_run_summary = {tempVar24};
        
        _pepxml->msms_run_summary[0]->spectrum_query = std::vector<pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_query*>(items.size());
        
        auto searchHits = std::vector<pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit*>();
        for (int i = 0; i < items.size(); i++)
        {
            auto mods = std::vector<pepXML::Generated::modInfoDataTypeMod_aminoacid_mass*>();
            auto alphaPeptide = items[i]->BestMatchingPeptides.First().Peptide;
            
            for (auto modification : alphaPeptide->AllModsOneIsNterminus)
            {
                auto mod = new pepXML::Generated::modInfoDataTypeMod_aminoacid_mass();
                mod->mass = modification->Value->MonoisotopicMass->Value;

                mod->position = (modification->Key - 1).ToString();
                mods.push_back(mod);
                
                //C# TO C++ CONVERTER TODO TASK: A 'delete mod' statement was not added since mod was passed to a
                //method or constructor. Handle memory management manually.
            }
            
            if (items[i]->getCrossType() == PsmCrossType::Single)
            {
                auto searchHit = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit();
                searchHit->hit_rank = 1;
                searchHit->peptide = alphaPeptide->BaseSequence;

                searchHit->peptide_prev_aa = alphaPeptide->PreviousAminoAcid.ToString();

                searchHit->peptide_next_aa = alphaPeptide->NextAminoAcid.ToString();
                searchHit->protein = alphaPeptide->Protein.Accession;
                searchHit->num_tot_proteins = 1;
                searchHit->calc_neutral_pep_mass = static_cast<float>(items[i]->getScanPrecursorMass());
                searchHit->massdiff = std::to_string(items[i]->getScanPrecursorMass() - items[i]->getPeptideMonisotopicMass()->Value);
                searchHit->xlink_typeSpecified = true;
                searchHit->xlink_type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type::na;
                searchHit->modification_info = new pepXML::Generated::modInfoDataType();
                searchHit->modification_info->mod_aminoacid_mass = mods.ToArray();
                pepXML::Generated::nameValueType *tempVar27 = new pepXML::Generated::nameValueType();
                tempVar27->name = "xlTotalScore";

                tempVar27->value = items[i]->getXLTotalScore().ToString();
                pepXML::Generated::nameValueType *tempVar28 = new pepXML::Generated::nameValueType();
                tempVar28->name = "Qvalue";

                tempVar28->value = items[i]->getFdrInfo().getQValue().ToString();
                searchHit->search_score = {tempVar27, tempVar28};
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
                auto mod = new pepXML::Generated::modInfoDataTypeMod_aminoacid_mass();
                mod->mass = crosslinkerDeadEndMass;

                mod->position = items[i]->getLinkPositions().front().ToString();
                mods.push_back(mod);
                auto searchHit = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit();
                searchHit->hit_rank = 1;
                searchHit->peptide = alphaPeptide->BaseSequence;

                searchHit->peptide_prev_aa = alphaPeptide->PreviousAminoAcid.ToString();

                searchHit->peptide_next_aa = alphaPeptide->NextAminoAcid.ToString();
                searchHit->protein = alphaPeptide->Protein.Accession;
                searchHit->num_tot_proteins = 1;
                searchHit->calc_neutral_pep_mass = static_cast<float>(items[i]->getScanPrecursorMass());
                searchHit->massdiff = std::to_string(items[i]->getScanPrecursorMass() - items[i]->getPeptideMonisotopicMass()->Value - crosslinkerDeadEndMass);
                searchHit->xlink_typeSpecified = true;
                searchHit->xlink_type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type::na;
                searchHit->modification_info = new pepXML::Generated::modInfoDataType();
                searchHit->modification_info->mod_aminoacid_mass = mods.ToArray();
                pepXML::Generated::nameValueType *tempVar29 = new pepXML::Generated::nameValueType();
                tempVar29->name = "xlTotalScore";

                tempVar29->value = items[i]->getXLTotalScore().ToString();
                pepXML::Generated::nameValueType *tempVar30 = new pepXML::Generated::nameValueType();
                tempVar30->name = "Qvalue";
                
                tempVar30->value = items[i]->getFdrInfo().getQValue().ToString();
                searchHit->search_score = {tempVar29, tempVar30};
                searchHits.push_back(searchHit);

                //C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since
                //searchHit was passed to a method or constructor. Handle memory management manually.
                //C# TO C++ CONVERTER TODO TASK: A 'delete mod' statement was not added since mod was
                //passed to a method or constructor. Handle memory management manually.
            }
            else if (items[i]->getCrossType() == PsmCrossType::Inter ||
                     items[i]->getCrossType() == PsmCrossType::Intra ||
                     items[i]->getCrossType() == PsmCrossType::Cross)
            {
                auto betaPeptide = items[i]->getBetaPeptide().BestMatchingPeptides.First().Peptide;
                auto modsBeta = std::vector<pepXML::Generated::modInfoDataTypeMod_aminoacid_mass*>();
                
                for (auto mod : betaPeptide->AllModsOneIsNterminus)
                {
                    auto modBeta = new pepXML::Generated::modInfoDataTypeMod_aminoacid_mass();
                    modBeta->mass = mod->Value->MonoisotopicMass->Value;

                    modBeta->position = (mod->Key - 1).ToString();
                    modsBeta.push_back(modBeta);
                    
                    //C# TO C++ CONVERTER TODO TASK: A 'delete modBeta' statement was not added since modBeta
                    //was passed to a method or constructor. Handle memory management manually.
                }
                
                auto alpha = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide();
                alpha->peptide = alphaPeptide->BaseSequence;
                
                alpha->peptide_prev_aa = alphaPeptide->PreviousAminoAcid.ToString();

                alpha->peptide_next_aa = alphaPeptide->NextAminoAcid.ToString();
                alpha->protein = alphaPeptide->Protein.Accession;
                alpha->num_tot_proteins = 1;
                alpha->calc_neutral_pep_mass = static_cast<float>(items[i]->getPeptideMonisotopicMass()->Value);
                alpha->complement_mass = static_cast<float>(items[i]->getScanPrecursorMass() - alphaPeptide->MonoisotopicMass);
                alpha->designation = "alpha";
                alpha->modification_info = new pepXML::Generated::modInfoDataType();
                alpha->modification_info->mod_aminoacid_mass = mods.ToArray();
                pepXML::Generated::nameValueType *tempVar31 = new pepXML::Generated::nameValueType();
                tempVar31->name = "xlscore";

                tempVar31->value = items[i]->getXLTotalScore().ToString();
                pepXML::Generated::nameValueType *tempVar32 = new pepXML::Generated::nameValueType();
                tempVar32->name = "link";

                tempVar32->value = items[i]->getLinkPositions().front().ToString();
                alpha->xlink_score = {tempVar31, tempVar32};
                auto beta = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide();
                beta->peptide = betaPeptide->BaseSequence;

                beta->peptide_prev_aa = betaPeptide->PreviousAminoAcid.ToString();

                beta->peptide_next_aa = betaPeptide->NextAminoAcid.ToString();
                beta->protein = betaPeptide->Protein.Accession;
                beta->num_tot_proteins = 1;
                beta->calc_neutral_pep_mass = static_cast<float>(betaPeptide->MonoisotopicMass);
                beta->complement_mass = static_cast<float>(items[i]->getScanPrecursorMass() - betaPeptide->MonoisotopicMass);
                beta->designation = "beta";
                beta->modification_info = new pepXML::Generated::modInfoDataType();
                beta->modification_info->mod_aminoacid_mass = modsBeta.ToArray();
                pepXML::Generated::nameValueType *tempVar33 = new pepXML::Generated::nameValueType();
                tempVar33->name = "xlscore";

                tempVar33->value = items[i]->getBetaPeptide().getScore().ToString();
                pepXML::Generated::nameValueType *tempVar34 = new pepXML::Generated::nameValueType();
                tempVar34->name = "link";

                tempVar34->value = items[i]->getBetaPeptide().getLinkPositions().front().ToString();
                beta->xlink_score = {tempVar33, tempVar34};
                auto cross = {alpha, beta};
                auto searchHit = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit();
                searchHit->hit_rank = 1;
                searchHit->peptide = "-";
                searchHit->peptide_prev_aa = "-";
                searchHit->peptide_next_aa = "-";
                searchHit->protein = "-";
                searchHit->num_tot_proteins = 1;
                searchHit->calc_neutral_pep_mass = static_cast<float>(items[i]->getScanPrecursorMass());
                searchHit->massdiff = std::to_string(items[i]->getScanPrecursorMass() - betaPeptide->MonoisotopicMass - alphaPeptide->MonoisotopicMass - crosslinker->getTotalMass());
                searchHit->xlink_typeSpecified = true;
                searchHit->xlink_type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type::xl;
                searchHit->xlink = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink();
                searchHit->xlink->identifier = crosslinker->getCrosslinkerName();
                searchHit->xlink->mass = static_cast<float>(crosslinker->getTotalMass());
                searchHit->xlink->linked_peptide = cross;
                pepXML::Generated::nameValueType *tempVar35 = new pepXML::Generated::nameValueType();
                tempVar35->name = "xlTotalScore";

                tempVar35->value = items[i]->getXLTotalScore().ToString();
                pepXML::Generated::nameValueType *tempVar36 = new pepXML::Generated::nameValueType();
                tempVar36->name = "Qvalue";

                tempVar36->value = items[i]->getFdrInfo().getQValue().ToString();
                searchHit->search_score = {tempVar35, tempVar36};
                searchHits.push_back(searchHit);
                
                //C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since
                //searchHit was passed to a method or constructor. Handle memory management manually.
                delete beta;
                delete alpha;
            }
            else if (items[i]->getCrossType() == PsmCrossType::Loop)
            {
                auto thePeptide = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide();
                pepXML::Generated::nameValueType *tempVar37 = new pepXML::Generated::nameValueType();
                tempVar37->name = "link";

                tempVar37->value = items[i]->getLinkPositions().front().ToString();
                pepXML::Generated::nameValueType *tempVar38 = new pepXML::Generated::nameValueType();
                tempVar38->name = "link";

                tempVar38->value = items[i]->getLinkPositions()[1].ToString();
                thePeptide->xlink_score = {tempVar37, tempVar38};
                auto cross = {thePeptide};
                auto searchHit = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit();
                searchHit->hit_rank = 1;
                searchHit->peptide = alphaPeptide->BaseSequence;

                searchHit->peptide_prev_aa = alphaPeptide->PreviousAminoAcid.ToString();

                searchHit->peptide_next_aa = alphaPeptide->NextAminoAcid.ToString();
                searchHit->protein = alphaPeptide->Protein.Accession;
                searchHit->num_tot_proteins = 1;
                searchHit->calc_neutral_pep_mass = static_cast<float>(items[i]->getScanPrecursorMass());
                searchHit->massdiff = std::to_string(items[i]->getScanPrecursorMass() - alphaPeptide->MonoisotopicMass - crosslinker->getLoopMass());
                searchHit->xlink_typeSpecified = true;
                searchHit->xlink_type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type::loop;
                searchHit->modification_info = new pepXML::Generated::modInfoDataType();
                searchHit->modification_info->mod_aminoacid_mass = mods.ToArray();
                searchHit->xlink = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink();
                searchHit->xlink->identifier = crosslinker->getCrosslinkerName();
                searchHit->xlink->mass = static_cast<float>(crosslinker->getTotalMass());
                searchHit->xlink->linked_peptide = cross;
                pepXML::Generated::nameValueType *tempVar39 = new pepXML::Generated::nameValueType();
                tempVar39->name = "xlTotalScore";

                tempVar39->value = items[i]->getXLTotalScore().ToString();
                pepXML::Generated::nameValueType *tempVar40 = new pepXML::Generated::nameValueType();
                tempVar40->name = "Qvalue";

                tempVar40->value = items[i]->getFdrInfo().getQValue().ToString();
                searchHit->search_score = {tempVar39, tempVar40};
                searchHits.push_back(searchHit);
                
                //C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since
                //searchHit was passed to a method or constructor. Handle memory management manually.
                delete thePeptide;
            }
        }
        
        for (int i = 0; i < items.size(); i++)
        {
            pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_query *tempVar41 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_query();

            tempVar41->spectrum = fileNameNoExtension + "." + items[i]->getScanNumber().ToString();
            tempVar41->start_scan = static_cast<unsigned int>(items[i]->getScanNumber());
            tempVar41->end_scan = static_cast<unsigned int>(items[i]->getScanNumber());
            tempVar41->precursor_neutral_mass = static_cast<float>(items[i]->getScanPrecursorMass());

            tempVar41->assumed_charge = items[i]->getScanPrecursorCharge().ToString();
            tempVar41->index = static_cast<unsigned int>(i + 1);
            tempVar41->retention_time_sec = static_cast<float>(items[i]->getScanRetentionTime() * 60);
            pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_result *tempVar42 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_result();
            tempVar42->search_hit = {searchHits[i]};
            tempVar41->search_result = {tempVar42};
            _pepxml->msms_run_summary[0].spectrum_query[i] = tempVar41;
        }
        
        TextWriter *writer = new StreamWriter(FileSystem::combine(outputFolder, fileName + ".pep.XM"));
        _indexedSerializer->Serialize(writer, _pepxml);
        writer->Close();
        FinishedWritingFile(FileSystem::combine(outputFolder, fileName + ".pep.XM"), nestedIds);
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete writer' statement was not added since writer was passed
        //to a method or constructor. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete _pepxml' statement was not added since _pepxml was passed
        //to a method or constructor. Handle memory management manually.
        delete _indexedSerializer;
    }
}
