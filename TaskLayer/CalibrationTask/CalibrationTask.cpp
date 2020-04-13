#include "CalibrationTask.h"
#include "CalibrationParameters.h"
#include "../../EngineLayer/CommonParameters.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"
#include "../MyTaskResults.h"
#include "../MyFileManager.h"
#include "../../EngineLayer/GlobalVariables.h"
#include "../../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/SingleAbsoluteAroundZeroMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../../EngineLayer/Calibration/CalibrationEngine.h"

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"

#include "MzML/Mzml.h"
#include "MzML/MzmlMethods.h"

using namespace UsefulProteomicsDatabases;
using namespace IO::MzML;

#include "stringhelper.h"

namespace TaskLayer
{

    CalibrationTask::CalibrationTask() : MetaMorpheusTask(MyTask::Calibrate)
    {
        
#ifdef ORIG
        EngineLayer::CommonParameters tempVar(, , false, , , , , , , , , , , false, , , new PpmTolerance(25), new PpmTolerance(15));
#endif
        std::string taskDescr = "";
        DissociationType dissType = DissociationType::HCD;
        bool useProvidedPrecursorInfo = true;
        double deconvIntensityRatio = 3;
        int deconvolutionMaxAssumedChargeState = 12;
        bool reportAllAmbiguity = true;
        bool addCompIons = false;
        int totalPartitions = 1;
        double scoreCutoff = 5;
        int topNpeaks = 200;
        double minRatio = 0.01;
        bool trimMs1Peaks = false;
        bool useDeltaScore = false;
        bool calculateEValue = false;
        EngineLayer::CommonParameters tempVar(taskDescr, dissType, false, useProvidedPrecursorInfo, deconvIntensityRatio,
                                              deconvolutionMaxAssumedChargeState, reportAllAmbiguity, addCompIons, totalPartitions,
                                              scoreCutoff, topNpeaks, minRatio, trimMs1Peaks, false, useDeltaScore, calculateEValue,
                                              new PpmTolerance(25), new PpmTolerance(15));
        setCommonParameters(&tempVar);
        
        TaskLayer::CalibrationParameters *tempVar2 = new TaskLayer::CalibrationParameters();
        setCalibrationParameters(tempVar2);
    }
    
    TaskLayer::CalibrationParameters *CalibrationTask::getCalibrationParameters() const
    {
        return privateCalibrationParameters;
    }
    
    void CalibrationTask::setCalibrationParameters(TaskLayer::CalibrationParameters *value)
    {
        privateCalibrationParameters = value;
    }
    
    MyTaskResults *CalibrationTask::RunSpecific(const std::string &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::string> &currentRawFileList, const std::string &taskId, std::vector<FileSpecificParameters*> &fileSettingsList)
    {
        std::vector<Modification*> variableModifications;
        std::vector<Modification*> fixedModifications;
        std::vector<std::string> localizeableModificationTypes;
        LoadModifications(taskId, variableModifications, fixedModifications, localizeableModificationTypes);
        
        // load proteins
        CommonParameters *cPar = getCommonParameters();
        std::vector<Protein*> proteinList = LoadProteins(taskId, dbFilenameList, true, DecoyType::Reverse, localizeableModificationTypes,
                                                         cPar);
        
        // write prose settings
        DigestionParams *dPar = getCommonParameters()->getDigestionParams();
        ProseCreatedWhileRunning->append("The following calibration settings were used: ");
        ProseCreatedWhileRunning->append("protease = " + dPar->getProtease()->ToString() + "; ");
        ProseCreatedWhileRunning->append("maximum missed cleavages = " + std::to_string(dPar->getMaxMissedCleavages())   + "; ");
        ProseCreatedWhileRunning->append("minimum peptide length = " + std::to_string(dPar->getMinPeptideLength()) + "; ");
        ProseCreatedWhileRunning->append(dPar->getMaxPeptideLength() == std::numeric_limits<int>::max() ?
                                         "maximum peptide length = unspecified; " : "maximum peptide length = " +
                                         std::to_string(dPar->getMaxPeptideLength()) + "; ");
        ProseCreatedWhileRunning->append("initiator methionine behavior = " +
                                         std::to_string(static_cast<int>(dPar->getInitiatorMethionineBehavior())) + "; ");
#ifdef ORIG
        ProseCreatedWhileRunning->append("fixed modifications = " + std::string::Join(", ", fixedModifications->Select([&] (std::any m) {
                        m::IdWithMotif;
                    })) + "; ");
#endif

        std::vector<std::string> vsv;
        for ( auto m : fixedModifications ) {
            vsv.push_back(m->getIdWithMotif());
        }
        std::string del = ", ";
        ProseCreatedWhileRunning->append("fixed modifications = " + StringHelper::join ( vsv, del) + "; ");
                                                                                         
        
#ifdef ORIG
        ProseCreatedWhileRunning->append("variable modifications = " + std::string::Join(", ", variableModifications->Select([&] (std::any m)  {
                        m::IdWithMotif;
                    })) + "; ");
#endif
        vsv.clear();
        for ( auto m : variableModifications ) {
            vsv.push_back(m->getIdWithMotif());
        }
        ProseCreatedWhileRunning->append("variable modifications = " + StringHelper::join ( vsv, del) + "; ");
        
        ProseCreatedWhileRunning->append("max mods per peptide = " + std::to_string(dPar->getMaxModsForPeptide()) + "; ");
        ProseCreatedWhileRunning->append("max modification isoforms = " +
                                         std::to_string(dPar->getMaxModificationIsoforms()) + "; ");
        ProseCreatedWhileRunning->append("precursor mass tolerance = " +
                                         std::to_string(getCommonParameters()->getPrecursorMassTolerance()->getValue()) + "; ");
        ProseCreatedWhileRunning->append("product mass tolerance = " +
                                         std::to_string(getCommonParameters()->getProductMassTolerance()->getValue()) + ". ");
#ifdef ORIG
        ProseCreatedWhileRunning->append("The combined search database contained " + proteinList.size()([&] (std::any p) {
                    !p::IsDecoy;
		}) + " non-decoy protein entries including " + proteinList.size()([&] (std::any p) {
			p::IsContaminant;
                    }) + " contaminant sequences. ");
#endif
        int c1=0, c2=0;
        for ( auto p: proteinList ) {
            if ( !p->getIsDecoy() ) {
                c1++;
            }
            if ( p->getIsContaminant() ) {
                c2++;
            }                  
        }
        ProseCreatedWhileRunning->append("The combined search database contained " + std::to_string(c1) +
                                         " non-decoy protein entries including " + std::to_string(c2) +
                                         " contaminant sequences. ");
                
        // start the calibration task
        std::vector<std::string> vs = {taskId};
        Status("Calibrating...", vs);
        myTaskResults = new MyTaskResults(this);
        myTaskResults->NewSpectra = std::vector<std::string>();
        myTaskResults->NewFileSpecificTomls = std::vector<std::string>();
        
        auto myFileManager = new MyFileManager(true);
        std::vector<std::string> spectraFilesAfterCalibration;
        
        for (int spectraFileIndex = 0; spectraFileIndex < (int)currentRawFileList.size(); spectraFileIndex++)
        {
            if (GlobalVariables::getStopLoops())
            {
                break;
            }
            
            bool couldNotFindEnoughDatapoints = false;
            
            // get filename stuff
            auto originalUncalibratedFilePath = currentRawFileList[spectraFileIndex];
            std::string originalUncalibratedFilenameWithoutExtension = Path::GetFileNameWithoutExtension(originalUncalibratedFilePath);
            std::string calibratedFilePath = FileSystem::combine(OutputFolder, originalUncalibratedFilenameWithoutExtension +
                                                                 CalibSuffix + ".mzM");
            
            // mark the file as in-progress
            std::vector<std::string> vs2 = {taskId, "Individual Spectra Files", originalUncalibratedFilePath};
            StartingDataFile(originalUncalibratedFilePath, vs2);
            
            EngineLayer::CommonParameters *combinedParams = SetAllFileSpecificCommonParams(getCommonParameters(), fileSettingsList[spectraFileIndex]);
            
            // load the file
            std::vector<std::string> vs3 = {taskId, "Individual Spectra Files"};
            Status("Loading spectra file...", vs3);
            
            auto myMsDataFile = myFileManager->LoadFile(originalUncalibratedFilePath,
                                                        std::make_optional(getCommonParameters()->getTopNpeaks()),
                                                        std::make_optional(getCommonParameters()->getMinRatio()),
                                                        getCommonParameters()->getTrimMs1Peaks(),
                                                        getCommonParameters()->getTrimMsMsPeaks(),
                                                        getCommonParameters());
            
            // get datapoints to fit calibration function to
            std::vector<std::string> vs4 = {taskId, "Individual Spectra Files"};
            Status("Acquiring calibration data points...", vs4);
            DataPointAquisitionResults *acquisitionResults = nullptr;
            
            for (int i = 1; i <= 5; i++)
            {
                acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications,
                                                               fixedModifications, proteinList, taskId, combinedParams,
                                                               combinedParams->getPrecursorMassTolerance(),
                                                               combinedParams->getProductMassTolerance());
                
                // enough data points to calibrate?
                if ((int)acquisitionResults->Psms.size() >= NumRequiredPsms                 &&
                    (int)acquisitionResults->getMs1List().size() > NumRequiredMs1Datapoints &&
                    (int)acquisitionResults->getMs2List().size() > NumRequiredMs2Datapoints)
                {
                    break;
                }
                
                if (i == 1) // failed round 1
                {
                    PpmTolerance tempVar(20);
                    getCommonParameters()->setPrecursorMassTolerance(&tempVar);
                    PpmTolerance tempVar2(50);
                    getCommonParameters()->setProductMassTolerance(&tempVar2);
                }
                else if (i == 2) // failed round 2
                {
                    PpmTolerance tempVar3(30);
                    getCommonParameters()->setPrecursorMassTolerance(&tempVar3);
                    PpmTolerance tempVar4(100);
                    getCommonParameters()->setProductMassTolerance(&tempVar4);
                }
                else if (i == 3) // failed round 3
                {
                    PpmTolerance tempVar5(40);
                    getCommonParameters()->setPrecursorMassTolerance(&tempVar5);
                    PpmTolerance tempVar6(150);
                    getCommonParameters()->setProductMassTolerance(&tempVar6);
                }
                else // failed round 4
                {
                    if ((int)acquisitionResults->Psms.size() < NumRequiredPsms)
                    {
                        Warn("Calibration failure! Could not find enough high-quality PSMs. Required " + std::to_string(NumRequiredPsms) +
                             ", saw " + std::to_string(acquisitionResults->Psms.size()));
                    }
                    if ((int)acquisitionResults->getMs1List().size() < NumRequiredMs1Datapoints)
                    {
                        Warn("Calibration failure! Could not find enough MS1 datapoints. Required " + std::to_string(NumRequiredMs1Datapoints) +
                             ", saw " + std::to_string(acquisitionResults->getMs1List().size()));
                    }
                    if ((int)acquisitionResults->getMs2List().size() < NumRequiredMs2Datapoints)
                    {
                        Warn("Calibration failure! Could not find enough MS2 datapoints. Required " + std::to_string(NumRequiredMs2Datapoints) +
                             ", saw " + std::to_string(acquisitionResults->getMs2List().size()));
                    }
                    
                    couldNotFindEnoughDatapoints = true;
                    std::vector<std::string> vsx = {taskId, "Individual Spectra Files", originalUncalibratedFilePath};
                    FinishedDataFile(originalUncalibratedFilePath, vsx);
                    break;
                }
                
                Warn("Could not find enough PSMs to calibrate with; opening up tolerances to " +
                     std::to_string(std::round(getCommonParameters()->getPrecursorMassTolerance()->getValue() * std::pow(10, 2)) / std::pow(10, 2)) +
                     " ppm precursor and " + std::to_string(std::round(getCommonParameters()->getProductMassTolerance()->getValue() * std::pow(10, 2)) / std::pow(10, 2)) +
                     " ppm product");
            }
            
            if (couldNotFindEnoughDatapoints)
            {
                spectraFilesAfterCalibration.push_back(Path::GetFileNameWithoutExtension(currentRawFileList[spectraFileIndex]));
                std::vector<std::string> vs = {taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension};
                ProgressEventArgs tempVar7(100, "Failed to calibrate!", vs);
                ReportProgress(&tempVar7);
                continue;
            }
            
            // stats before calibration
            int prevPsmCount = acquisitionResults->Psms.size();
            double preCalibrationPrecursorErrorIqr = acquisitionResults->PsmPrecursorIqrPpmError;
            double preCalibrationProductErrorIqr = acquisitionResults->PsmProductIqrPpmError;
            
            // generate calibration function and shift data points
            std::vector<std::string> vs5 = {taskId, "Individual Spectra Files"};
            Status("Calibrating...", vs5);
            std::vector<std::string> vs6 = {taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension};
            CalibrationEngine *engine = new CalibrationEngine(myMsDataFile, acquisitionResults, getCommonParameters(), vs6);
            engine->Run();
            
            //update file
            myMsDataFile = engine->getCalibratedDataFile();
            
            // do another search to evaluate calibration results
            std::vector<std::string> vs7 = {taskId, "Individual Spectra Files"};
            Status("Post-calibration search...", vs7 );
            acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications,
                                                           fixedModifications, proteinList, taskId, combinedParams,
                                                           combinedParams->getPrecursorMassTolerance(),
                                                           combinedParams->getProductMassTolerance());
            
            //generate calibration function and shift data points AGAIN because it's fast and contributes new data
            std::vector<std::string> vs8 = {taskId, "Individual Spectra Files"};
            Status("Calibrating...", vs8);
            std::vector<std::string> vs9 = {taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension};
            engine = new CalibrationEngine(myMsDataFile, acquisitionResults, getCommonParameters(), vs9);
            engine->Run();
            
            //update file
            myMsDataFile = engine->getCalibratedDataFile();
            
            // write the calibrated mzML file
            //arguments:  MsDataFile *myMsDataFile, const std::string &outputFile, bool writeIndexed
            MzmlMethods::CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, calibratedFilePath, false);
            myFileManager->DoneWithFile(originalUncalibratedFilePath);
            
            // stats after calibration
            int postCalibrationPsmCount = acquisitionResults->Psms.size();
            double postCalibrationPrecursorErrorIqr = acquisitionResults->PsmPrecursorIqrPpmError;
            double postCalibrationProductErrorIqr = acquisitionResults->PsmProductIqrPpmError;
            
            // did the data improve? (not used for anything yet...)
            bool improvement = ImprovGlobal(preCalibrationPrecursorErrorIqr, preCalibrationProductErrorIqr, prevPsmCount,
                                            postCalibrationPsmCount, postCalibrationPrecursorErrorIqr,
                                            postCalibrationProductErrorIqr);
            
            // write toml settings for the calibrated file
            auto newTomlFileName = FileSystem::combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".toml");
            
            auto fileSpecificParams = new FileSpecificParameters();
            
            // carry over file-specific parameters from the uncalibrated file to the calibrated one
            if (fileSettingsList[spectraFileIndex] != nullptr)
            {
                fileSpecificParams = fileSettingsList[spectraFileIndex]->Clone();
            }
            
            //suggest 4 * interquartile range as the ppm tolerance
            PpmTolerance tempVar8((4.0 * postCalibrationPrecursorErrorIqr) + std::abs(acquisitionResults->PsmPrecursorMedianPpmError));
            fileSpecificParams->setPrecursorMassTolerance(&tempVar8);
            PpmTolerance tempVar9((4.0 * postCalibrationProductErrorIqr) + std::abs(acquisitionResults->PsmProductMedianPpmError));
            fileSpecificParams->setProductMassTolerance(&tempVar9);
            
#ifdef ORIG
            Toml::WriteFile(fileSpecificParams, newTomlFileName, tomlConfig);
#endif

            toml::Table fileSpecificParamsTable = fileSpecificParams->getFileSpecificParametersTomlTable();

            Toml trw;
            trw.tomlWriteNewFile(newTomlFileName, fileSpecificParamsTable);

            
            std::vector<std::string> vs10 = {taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension};
            FinishedWritingFile(newTomlFileName, vs10);
            
            // finished calibrating this file
            spectraFilesAfterCalibration.push_back(Path::GetFileNameWithoutExtension(calibratedFilePath));
            std::vector<std::string> vs11 = {taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension};
            FinishedWritingFile(calibratedFilePath, vs11);
            myTaskResults->NewSpectra.push_back(calibratedFilePath);
            myTaskResults->NewFileSpecificTomls.push_back(newTomlFileName);
            std::vector<std::string> vs12 = {taskId, "Individual Spectra Files", originalUncalibratedFilePath};
            FinishedDataFile(originalUncalibratedFilePath, vs12);
            std::vector<std::string> vs13 = {taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension};
            ProgressEventArgs tempVar10(100, "Done!", vs13);
            ReportProgress(&tempVar10);
            
            //C# TO C++ CONVERTER TODO TASK: A 'delete fileSpecificParams' statement was not added since
            //fileSpecificParams was passed to a method or constructor. Handle memory management manually.
            delete engine;
        }
        
        // re-write experimental design (if it has been defined) with new calibrated file names
        std::string assumedPathToExperDesign = Directory::GetParent(currentRawFileList.front())->FullName;
        assumedPathToExperDesign = FileSystem::combine(assumedPathToExperDesign, GlobalVariables::getExperimentalDesignFileName());
        
        if (FileSystem::fileExists(assumedPathToExperDesign))
        {
            WriteNewExperimentalDesignFile(assumedPathToExperDesign, OutputFolder, spectraFilesAfterCalibration);
        }
        
        // finished calibrating all files for the task
        std::vector<std::string> vs14 = {taskId, "Individual Spectra Files"};
        ProgressEventArgs tempVar11(100, "Done!", vs14);
        ReportProgress(&tempVar11);
        
        delete myFileManager;
        return myTaskResults;
    }
    
    const std::string CalibrationTask::CalibSuffix = "-calib";
    
    bool CalibrationTask::ImprovGlobal(double prevPrecTol, double prevProdTol, int prevPsmCount, int thisRoundPsmCount, double thisRoundPrecTol, double thisRoundProdTol)
    {
        if (thisRoundPsmCount > prevPsmCount)
        {
            return true;
        }
        
        auto precRatio = thisRoundPrecTol / prevPrecTol;
        auto prodRatio = thisRoundProdTol / prevProdTol;
        
        if (thisRoundPsmCount == prevPsmCount)
        {
            return precRatio + prodRatio < 2; // Take any improvement in ratios
        }
        
        auto countRatio = static_cast<double>(thisRoundPsmCount) / prevPsmCount;
        return countRatio > 0.9 && precRatio + prodRatio < 1.8;
    }
    
    DataPointAquisitionResults *CalibrationTask::GetDataAcquisitionResults(MsDataFile *myMsDataFile, const std::string &currentDataFile, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<Protein*> &proteinList, const std::string &taskId, EngineLayer::CommonParameters *combinedParameters, Tolerance *initPrecTol, Tolerance *initProdTol)
    {
        std::string fileNameWithoutExtension = Path::GetFileNameWithoutExtension(currentDataFile);
        SinglePpmAroundZeroSearchMode tempVar(initPrecTol->getValue());
        MassDiffAcceptor *searchMode = dynamic_cast<PpmTolerance*>(initPrecTol) != nullptr ? static_cast<MassDiffAcceptor*>(&tempVar):
            new SingleAbsoluteAroundZeroSearchMode(initPrecTol->getValue());
        
#ifdef ORIG
        auto listOfSortedms2Scans = GetMs2Scans(myMsDataFile, currentDataFile, combinedParameters).OrderBy([&] (std::any b)  {
                b::PrecursorMass;
            })->ToArray();
#endif
        std::vector<Ms2ScanWithSpecificMass*> listOfSortedms2Scans =  GetMs2Scans(myMsDataFile, currentDataFile, combinedParameters);
        std::sort(listOfSortedms2Scans.begin(), listOfSortedms2Scans.end(),[&] (Ms2ScanWithSpecificMass*l, Ms2ScanWithSpecificMass* r) {
                return l->getPrecursorMass() > r->getPrecursorMass(); });
        std::vector<PeptideSpectralMatch*> allPsmsArray(listOfSortedms2Scans.size());
        
        std::vector<std::string> vs1 = {taskId, "Individual Spectra Files", fileNameWithoutExtension};
        Log("Searching with searchMode: " + searchMode->ToProseString(), vs1);
        std::vector<std::string> vs2 = {taskId, "Individual Spectra Files", fileNameWithoutExtension};
        Log("Searching with productMassTolerance: " + std::to_string(initProdTol->getValue()), vs2);
        
        std::vector<std::string> vs3 = {taskId, "Individual Spectra Files", fileNameWithoutExtension};
        ClassicSearchEngine tempVar2(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList,
                                     searchMode, combinedParameters, vs3);
        (&tempVar2)->Run();
#ifdef ORIG
        std::vector<PeptideSpectralMatch*> allPsms = allPsmsArray.Where([&] (std::any b) {
                return b != nullptr;
            }).ToList();
#endif
        std::vector<PeptideSpectralMatch*> allPsms;
        for ( auto b: allPsmsArray ) {
            if ( b != nullptr ) {
                allPsms.push_back(b);
            }
        }            
#ifdef ORIG
        allPsms = allPsms.OrderByDescending([&] (std::any b)  {
                b::Score;
            }).ThenBy([&] (std::any b) {
                    b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
                }).GroupBy([&] (std::any b) {
                        (b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass);
                    })->Select([&] (std::any b){
                            b::First();
                        }).ToList();
#endif
        // OrderByDescending and ThenBY
        std::sort(allPsms.begin(), allPsms.end(), [&] (PeptideSpectralMatch* r, PeptideSpectralMatch* l) {
                if ( r->getScore() > l->getScore() ) return true;
                if ( r->getScore() < l->getScore() ) return false;
                double lval= std::numeric_limits<double>::max(), rval=std::numeric_limits<double>::max();
                if ( l->getPeptideMonisotopicMass().has_value() ) {
                    lval = std::abs(l->getScanPrecursorMass() - l->getPeptideMonisotopicMass().value());
                }
                if ( r->getPeptideMonisotopicMass().has_value() ) {
                    rval = std::abs(r->getScanPrecursorMass() - r->getPeptideMonisotopicMass().value());
                }
                return rval < lval;                    
            });
        //GroupBy tuple
        std::vector<std::vector<PeptideSpectralMatch *>>tvec;
        for ( auto psm : allPsms ) {
            for ( auto t : tvec ) {
                if ( t[0]->getFullFilePath() == psm->getFullFilePath()                      &&
                     t[0]->getScanNumber()   == psm->getScanNumber()                        &&
                     t[0]->getPeptideMonisotopicMass() == psm->getPeptideMonisotopicMass() ) {
                    t.push_back(psm);
                }
                else {
                    std::vector<PeptideSpectralMatch *> *t = new std::vector<PeptideSpectralMatch *>;
                    t->push_back(psm);
                    tvec.push_back(*t);
                }
            }
        }
        // Select first
        allPsms.clear();
        for ( auto t: tvec ) {
            allPsms.push_back(t[0]);
        }
        
        std::vector<std::string> vs4 = {taskId, "Individual Spectra Files", fileNameWithoutExtension};
        FdrAnalysisEngine tempVar3(allPsms, searchMode->getNumNotches(), getCommonParameters(), vs4);
        (&tempVar3)->Run();
        
#ifdef ORIG
        std::vector<PeptideSpectralMatch*> goodIdentifications = allPsms.Where([&] (std::any b) {
                return b::FdrInfo::QValueNotch < 0.001 && !b::IsDecoy && b::FullSequence != nullptr;
            }).ToList();
#endif
        std::vector<PeptideSpectralMatch*> goodIdentifications;
        for ( auto b: allPsms ) {
            if ( b->getFdrInfo()->getQValueNotch() < 0.001 &&
                 !b->getIsDecoy()                     &&
                 !b->getFullSequence().empty() ) {
                goodIdentifications.push_back(b);
            }
        }
        
        if (goodIdentifications.empty())
        {
            std::vector<PeptideSpectralMatch*> vp;
            std::vector<LabeledDataPoint*> vldp1;
            std::vector<LabeledDataPoint*> vldp2;
            return new DataPointAquisitionResults(nullptr, vp, vldp1, vldp2, 0, 0, 0, 0);
        }
        
        std::vector<std::string> vs5 = {taskId, "Individual Spectra Files", fileNameWithoutExtension};
        DataPointAcquisitionEngine tempVar4(goodIdentifications, myMsDataFile, initPrecTol,
                                            getCalibrationParameters()->getMinMS1IsotopicPeaksNeededForConfirmedIdentification(),
                                            getCommonParameters(), vs5);
        DataPointAquisitionResults *currentResult = static_cast<DataPointAquisitionResults*>((&tempVar4)->Run());
        
        return currentResult;
    }
    
    void CalibrationTask::WriteNewExperimentalDesignFile(const std::string &assumedPathToExperDesign, const std::string &outputFolder, std::vector<std::string> &spectraFilesAfterCalibration)
    {
        std::vector<std::string> lines = File::ReadAllLines(assumedPathToExperDesign);
        //std::vector<std::string> newExperimentalDesignOutput = std::vector<std::vector<std::string>>(0) };
        std::vector<std::string> newExperimentalDesignOutput;
        
        for (int i = 1; i < (int)lines.size(); i++)
        {
            auto split = StringHelper::split(lines[i], L'\t');
            std::string oldFileName = Path::GetFileNameWithoutExtension(split[0]);
            std::string newFileName = oldFileName + CalibSuffix;
            std::string newline;
            
            if (std::find(spectraFilesAfterCalibration.begin(), spectraFilesAfterCalibration.end(), newFileName) == spectraFilesAfterCalibration.end())
            {
                // file was not successfully calibrated
                newline = oldFileName + "\t";
            }
            else
            {
                // file was successfully calibrated
                newline = newFileName + "\t";
            }
            
            // add condition, biorep, etc info
            for (int j = 1; j < (int)split.size(); j++)
            {
                newline += split[j] + "\t";
            }
            
            // write the line
            newExperimentalDesignOutput.push_back(newline);
        }
        
        File::WriteAllLines(FileSystem::combine(outputFolder, GlobalVariables::getExperimentalDesignFileName()), newExperimentalDesignOutput);
    }
}

