#include "MetaMorpheusTask.h"
#include "MyTaskResults.h"
#include "FileSpecificParameters.h"
#include "DbForTask.h"
#include "EventArgs/SingleTaskEventArgs.h"

#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/MetaMorpheusException.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../EngineLayer/EventArgs/SingleFileEventArgs.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "../EngineLayer/EventArgs/SingleEngineFinishedEventArgs.h"
#include "../EngineLayer/Indexing/IndexingResults.h"

#include "Chemistry/ClassExtensions.h"
#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
#include "MassSpectrometry/Enums/DissociationType.h"
#include "Proteomics/ProteolyticDigestion/CleavageSpecificity.h"

#include <ctime>
#include <filesystem>
#include <exception>
#include <string>
#include <locale>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <cmath>
#include <limits>

#include <stdio.h>

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

static int MetaMorpheus_offset=0;
#define MetaMorpheus_OffsetIncrement 4
#define DISPLAY_OFFSET(_off) {              \
      if ( MetaMorpheus_offset == 4 )       \
          std::cout << "    ";              \
      else if ( MetaMorpheus_offset == 8 )  \
          std::cout << "        ";          \
      else if ( MetaMorpheus_offset == 12 ) \
          std::cout << "            ";      \
}

#define OFFSET_INCR(_off) {_off+=MetaMorpheus_OffsetIncrement;}
#define OFFSET_DECR(_off) {_off-=MetaMorpheus_OffsetIncrement;}

static std::unordered_set<std::string> EngineProgressKeys;


#ifdef TIMING_INFO
#include <sys/time.h>

double deconvtime=0.0;
double sorttime = 0.0;
double looptime = 0.0;

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
    void MetaMorpheusTask::writeTomlConfig (std::string &filename, std::ofstream &tomlFd )
    {
        if ( !tomlFd.is_open() ) {
            tomlFd.open(filename, std::ios_base::app );
            if ( !tomlFd.is_open() ) {
                std::cout << "Could not open file " << filename << std::endl;
                return;
            }
        }

        CommonParameters *cparams = getCommonParameters();
        toml::Table common_params;
        common_params["MaxThreadsToUsePerFile"] = cparams->getMaxThreadsToUsePerFile();

        auto var1 = cparams->getListOfModsFixed();
        std::string del  = "\t\t";
        std::vector<std::string> vecvar1;
        if ( var1 != nullptr ) {
            for ( auto v : *var1 ) {
                std::string s = std::get<0>(v) + "\t" + std::get<1>(v);
                vecvar1.push_back(s);
            }
            std::string res1 = StringHelper::join(vecvar1, del);
            common_params["ListOfModsFixed"] = res1;
        }

        auto var2 = cparams->getListOfModsVariable();
        std::vector<std::string> vecvar2;
        if ( var2 != nullptr ) {
            for ( auto v : *var2 ) {
                std::string s = std::get<0>(v) + "\t" + std::get<1>(v);
                vecvar2.push_back(s);
            }
            std::string res2 = StringHelper::join(vecvar2, del);            
            common_params["ListOfModsVariable"] = res2;
        }

        common_params["DoPrecursorDeconvolution"] = cparams->getDoPrecursorDeconvolution();
        common_params["UseProvidedPrecursorInfo"] = cparams->getUseProvidedPrecursorInfo();
        common_params["DeconvolutionIntensityRatio"] = cparams->getDeconvolutionIntensityRatio();
        common_params["DeconvolutionMaxAssumedChargeState"] = cparams->getDeconvolutionMaxAssumedChargeState();

        std::string res3 = StringHelper::formatSimple("{0}{1} PPM", "±",
                                                      cparams->getDeconvolutionMassTolerance()->getValue());
        common_params["DeconvolutionMassTolerance"] = res3;
        common_params["TotalPartitions"] = cparams->getTotalPartitions();

        std::string res4 = StringHelper::formatSimple("{0}{1} PPM", "±",
                                                      cparams->getProductMassTolerance()->getValue());
        common_params["ProductMassTolerance"] = res4 ;

        std::string res5 = StringHelper::formatSimple("{0}{1} PPM", "±",
                                                      cparams->getPrecursorMassTolerance()->getValue());
        common_params["PrecursorMassTolerance"] = res5;
        common_params["AddCompIons"] = cparams->getAddCompIons();
        common_params["ScoreCutoff"] = cparams->getScoreCutoff();
        common_params["ReportAllAmbiguity"] = cparams->getReportAllAmbiguity();
        common_params["TopNpeaks"] = cparams->getTopNpeaks();
        common_params["MinRatio"] = cparams->getMinRatio();
        common_params["TrimMs1Peaks"] = cparams->getTrimMs1Peaks();
        common_params["TrimMsMsPeaks"] = cparams->getTrimMsMsPeaks();
        common_params["UseDeltaScore"] = cparams->getUseDeltaScore();
        common_params["CalculateEValue"] =  cparams->getCalculateEValue();
        common_params["QValueOutputFilter"] = cparams->getQValueOutputFilter();

        common_params["DissociationType"] = GetDissocationType::GetDissocationTypeAsString(cparams->getDissociationType());
        common_params["AssumeOrphanPeaksAreZ1Fragments"] = cparams->getAssumeOrphanPeaksAreZ1Fragments();
        common_params["MaxHeterozygousVariants"] = cparams->getMaxHeterozygousVariants();
        common_params["MinVariantDepth"] = cparams->getMinVariantDepth();

        tomlFd << std::endl;
        tomlFd << "[CommonParameters]" << std::endl;
        tomlFd << common_params ;

        toml::Table digest_params;
        DigestionParams *dparams = cparams->getDigestionParams();
        
        digest_params["MaxMissedCleavages"] = dparams->getMaxMissedCleavages();
        digest_params["InitiatorMethionineBehavior"] = ProteolyticDigestion::InitiatorMethionineBehaviorToString(dparams->getInitiatorMethionineBehavior());
        digest_params["MinPeptideLength"] = dparams->getMinPeptideLength();
        digest_params["MaxPeptideLength"] = dparams->getMaxPeptideLength();
        digest_params["MaxModificationIsoforms"] = dparams->getMaxModificationIsoforms();
        digest_params["MaxModsForPeptide"] = dparams->getMaxModsForPeptide();
        digest_params["Protease"] = dparams->getProtease()->ToString();
        digest_params["SearchModeType"] = ProteolyticDigestion::CleavageSpecificityExtension::GetCleavageSpecificityAsString(dparams->getSearchModeType());
        digest_params["FragmentationTerminus"] = Fragmentation:: FragmentationTerminusToString(dparams->getFragmentationTerminus());
        digest_params["SpecificProtease"] = dparams->getSpecificProtease()->ToString();

        tomlFd << std::endl;
        tomlFd << "[CommonParameters.DigestionParams]" << std::endl;
        tomlFd << digest_params;

        tomlFd.close();
        
        return;
    }
    MetaMorpheusTask::MetaMorpheusTask(MyTask taskType)
    {
        this->setTaskType(taskType);
    }
    
    MyTask MetaMorpheusTask::getTaskType() const
    {
        return privateTaskType;
    }
    
    void MetaMorpheusTask::setTaskType(MyTask value)
    {
        privateTaskType = value;
    }
    
    EngineLayer::CommonParameters *MetaMorpheusTask::getCommonParameters() const
    {
        return privateCommonParameters;
    }

    int  MetaMorpheusTask::getVerbose() const
    {
        return privateVerbosityLevel;
    }

    void MetaMorpheusTask::setVerbose( int verbosityLevel )
    {
        privateVerbosityLevel = verbosityLevel;
    }
    
    void MetaMorpheusTask::setCommonParameters(EngineLayer::CommonParameters *value)
    {
        if ( privateCommonParameters != nullptr ) {
            delete privateCommonParameters;
        }
        privateCommonParameters = value;
    }
    
    const std::string MetaMorpheusTask::IndexFolderName = "DatabaseIndex";
    
    std::vector<Ms2ScanWithSpecificMass*> MetaMorpheusTask::GetMs2Scans(MsDataFile *myMSDataFile,
                                                                        const std::string &fullFilePath,
                                                                        EngineLayer::CommonParameters *commonParameters,
                                                                        int firstIndex, int lastIndex )
    {
#ifdef TIMING_INFO
        struct timeval t1, t1e;
        double tgobs=0.0, tgobs2=0.0, tgnef=0.0, tgnef2=0.0;
        int tgobscount=0, tgobscount2=0, tgnefcount=0, tgnefcount2=0;
        gettimeofday (&t1, NULL);
#endif
        std::vector<MsDataScan*> ms2Scans;
        if ( lastIndex == -1 ) {
            // Read entire file.
            for ( auto x : myMSDataFile->GetAllScansList() ) {
                if ( x->getMsnOrder() > 1 ) {
                    ms2Scans.push_back(x);
                }
            }
        }
        else {
            for ( auto x : myMSDataFile->GetMsScansSubset(firstIndex, lastIndex ) ) {
                if ( x->getMsnOrder() > 1 ) {
                    ms2Scans.push_back(x);
                }
            }
        }
        
        std::vector<std::vector<Ms2ScanWithSpecificMass*>> scansWithPrecursors(ms2Scans.size());
        for ( int i = 0; i < (int)ms2Scans.size(); i++ ) {            
            MsDataScan *ms2scan = ms2Scans[i];
            
            std::vector<std::pair<double, int>> precursors;
            if (ms2scan->getOneBasedPrecursorScanNumber().has_value())
            {
#ifdef TIMING_INFO
                struct timeval tx, txe;
                gettimeofday (&tx, NULL);
#endif
                auto precursorSpectrum = myMSDataFile->GetOneBasedScan(ms2scan->getOneBasedPrecursorScanNumber().value());
                
                try
                {
                    ms2scan->RefineSelectedMzAndIntensity(precursorSpectrum->getMassSpectrum());
                }
                catch (const MzLibException &ex)
                {
                    //Warn("Could not get precursor ion for MS2 scan #" + ms2scan->getOneBasedScanNumber() + "; " + ex.Message);
                    std::string s = "Could not get precursor ion for MS2 scan #";
                    s += ms2scan->getOneBasedScanNumber() + "; ";
                    Warn(s);
                    continue;
                }
                if (ms2scan->getSelectedIonMonoisotopicGuessMz().has_value())
                {
                    ms2scan->ComputeMonoisotopicPeakIntensity(precursorSpectrum->getMassSpectrum());
                }
                
                if (commonParameters->getDoPrecursorDeconvolution())
                {
                    for (auto envelope : ms2scan->GetIsolatedMassesAndCharges(precursorSpectrum->getMassSpectrum(), 1,
                                                                           commonParameters->getDeconvolutionMaxAssumedChargeState(),
                                                                           commonParameters->getDeconvolutionMassTolerance()->getValue(),
                                                                           commonParameters->getDeconvolutionIntensityRatio()))
                    {
                        auto monoPeakMz = Chemistry::ClassExtensions::ToMz(envelope->monoisotopicMass, envelope->charge);
                        precursors.push_back(std::make_pair(monoPeakMz, envelope->charge));
                    }
                }
#ifdef TIMING_INFO
                gettimeofday (&txe, NULL);
                tgobs += timediff (tx, txe);
                tgobscount++;
#endif
            }
            
            if (commonParameters->getUseProvidedPrecursorInfo()      &&
                ms2scan->getSelectedIonChargeStateGuess().has_value() )
            {
#ifdef TIMING_INFO
                struct timeval tx, txe;
                gettimeofday (&tx, NULL);
#endif
                auto precursorCharge = ms2scan->getSelectedIonChargeStateGuess().value();
                if (ms2scan->getSelectedIonMonoisotopicGuessMz().has_value())
                {
                    double precursorMZ = ms2scan->getSelectedIonMonoisotopicGuessMz().value();
                    bool found = false;
                    auto tol = commonParameters->getDeconvolutionMassTolerance();
                    for ( auto b: precursors ) {
                        if ( tol->Within(Chemistry::ClassExtensions::ToMass(precursorMZ, precursorCharge),
                                         Chemistry::ClassExtensions::ToMass(b.first, b.second ))) {
                            found = true;
                            break;
                        }
                    }
                    if ( !found )  {
                        precursors.push_back(std::make_pair(precursorMZ, precursorCharge));
                    }
                }
                else
                {
                    double precursorMZ = ms2scan->getSelectedIonMZ().value();
                    bool found = false;
                    auto tol = commonParameters->getDeconvolutionMassTolerance();
                    for ( auto b: precursors ) {
                        if ( tol->Within(Chemistry::ClassExtensions::ToMass(precursorMZ, precursorCharge),
                                         Chemistry::ClassExtensions::ToMass(b.first, b.second ))) {
                            found = true;
                            break;
                        }
                    }
                    if ( !found )       {
                        precursors.push_back(std::make_pair(precursorMZ, precursorCharge));
                    }
                }
#ifdef TIMING_INFO
                gettimeofday (&txe, NULL);
                tgobs2 += timediff (tx, txe);
                tgobscount2++;
#endif
            }
            
#ifdef TIMING_INFO
            struct timeval t2, t2e;
            gettimeofday (&t2, NULL);
#endif
            std::vector<IsotopicEnvelope*> neutralExperimentalFragments = Ms2ScanWithSpecificMass::GetNeutralExperimentalFragments(ms2scan,
                                                                                                                          commonParameters);
#ifdef TIMING_INFO
            gettimeofday (&t2e, NULL);
            tgnef += timediff (t2, t2e);
            tgnefcount++;
            struct timeval t3, t3e;
            gettimeofday (&t3, NULL);
#endif
            
            for (auto precursor : precursors)
            {
                auto  tempVar2 = new Ms2ScanWithSpecificMass(ms2scan, precursor.first, precursor.second, fullFilePath,
                                                             commonParameters, neutralExperimentalFragments);
                scansWithPrecursors[i].push_back(tempVar2);
#ifdef TIMING_INFO
                tgnefcount2++;
#endif
            }
#ifdef TIMING_INFO
            gettimeofday (&t3e, NULL);
            tgnef2 += timediff (t3, t3e);
#endif
        } 
    
        std::vector<Ms2ScanWithSpecificMass*> tmpv;
        for ( auto p: scansWithPrecursors ) {
            for ( auto v: p ) {
                tmpv.push_back(v);
            }
        }

#ifdef TIMING_INFO
        gettimeofday (&t1e, NULL);
        std::cout << " Total time spent in GetMs2Scans : " << timediff( t1, t1e ) << std::endl;
        std::cout << "       time spent in block1      : " << tgobs << " called " << tgobscount << " times" << std::endl;
        std::cout << "       time spent in block2      : " << tgobs2 << " called " << tgobscount2 << " times" << std::endl;
        std::cout << "       time spent in GetNeutralExperimentalFragments : " << tgnef << " called " << tgnefcount << " times" << std::endl;
        std::cout << "            time spent in Deconvolution  : " << deconvtime << std::endl;        
        std::cout << "            time spent in loop           : " << looptime << std::endl;        
        std::cout << "            time spent in sort           : " << sorttime << std::endl;        
        std::cout << "       time spent in Ms2ScanWithSpecificMass constructor : " << tgnef2 << " called " << tgnefcount2 << " times" << std::endl;
#endif
        return tmpv;
    }
    
    EngineLayer::CommonParameters *MetaMorpheusTask::SetAllFileSpecificCommonParams(EngineLayer::CommonParameters *commonParams,
                                                                                    FileSpecificParameters *fileSpecificParams)
    {
        if (fileSpecificParams == nullptr)
        {
            return commonParams;
        }
        
        // set file-specific digestion parameters
        Protease *tempVar = fileSpecificParams->getProtease();
        Protease *protease = (tempVar != nullptr) ? tempVar : commonParams->getDigestionParams()->getProtease();
        std::optional<int> tempVar2 = fileSpecificParams->getMinPeptideLength();
        int minPeptideLength = tempVar2.has_value() ? tempVar2.value() : commonParams->getDigestionParams()->getMinPeptideLength();
        std::optional<int> tempVar3 = fileSpecificParams->getMaxPeptideLength();
        int maxPeptideLength = tempVar3.has_value() ? tempVar3.value() : commonParams->getDigestionParams()->getMaxPeptideLength();
        std::optional<int> tempVar4 = fileSpecificParams->getMaxMissedCleavages();
        int maxMissedCleavages = tempVar4.has_value() ? tempVar4.value() : commonParams->getDigestionParams()->getMaxMissedCleavages();
        std::optional<int> tempVar5 = fileSpecificParams->getMaxModsForPeptide();
        int maxModsForPeptide = tempVar5.has_value() ? tempVar5.value() : commonParams->getDigestionParams()->getMaxModsForPeptide();
        
        DigestionParams *fileSpecificDigestionParams = new DigestionParams(protease->getName(),
                                                                     maxMissedCleavages,
                                                                     minPeptideLength,
                                                                     maxPeptideLength,
                                                                     commonParams->getDigestionParams()->getMaxModificationIsoforms(),
                                                                     commonParams->getDigestionParams()->getInitiatorMethionineBehavior(),
                                                                     maxModsForPeptide,
                                                                     commonParams->getDigestionParams()->getSearchModeType(),
                                                                     commonParams->getDigestionParams()->getFragmentationTerminus());
        
        // set the rest of the file-specific parameters
        Tolerance *tempVar6 = fileSpecificParams->getPrecursorMassTolerance();
        Tolerance *precursorMassTolerance = (tempVar6 != nullptr) ? tempVar6 : commonParams->getPrecursorMassTolerance();
        Tolerance *tempVar7 = fileSpecificParams->getProductMassTolerance();
        Tolerance *productMassTolerance = (tempVar7 != nullptr) ? tempVar7 : commonParams->getProductMassTolerance();
        std::optional<DissociationType*> tempVar8 = fileSpecificParams->getDissociationType();
        DissociationType *tempVar8val = new DissociationType();
        *tempVar8val = commonParams->getDissociationType();
        DissociationType *dissociationType = tempVar8.has_value() ? tempVar8.value() : tempVar8val;
        
        
        EngineLayer::CommonParameters *returnParams = new CommonParameters(commonParams->getTaskDescriptor(),
                                                                           *dissociationType,
                                                                           commonParams->getDoPrecursorDeconvolution(),
                                                                           commonParams->getUseProvidedPrecursorInfo(),
                                                                           commonParams->getDeconvolutionIntensityRatio(),
                                                                           commonParams->getDeconvolutionMaxAssumedChargeState(),
                                                                           commonParams->getReportAllAmbiguity(),
                                                                           commonParams->getAddCompIons(),
                                                                           commonParams->getTotalPartitions(),
                                                                           commonParams->getScoreCutoff(),
                                                                           commonParams->getTopNpeaks(),
                                                                           commonParams->getMinRatio(),
                                                                           commonParams->getTrimMs1Peaks(),
                                                                           commonParams->getTrimMsMsPeaks(),
                                                                           commonParams->getUseDeltaScore(),
                                                                           commonParams->getCalculateEValue(),
                                                                           productMassTolerance,
                                                                           precursorMassTolerance,
                                                                           commonParams->getDeconvolutionMassTolerance(),
                                                                           commonParams->getMaxThreadsToUsePerFile(),
                                                                           fileSpecificDigestionParams,
                                                                           commonParams->getListOfModsVariable(),
                                                                           commonParams->getListOfModsFixed(),
                                                                           commonParams->getQValueOutputFilter(),
                                                                           commonParams->getAssumeOrphanPeaksAreZ1Fragments());

        //delete fileSpecificDigestionParams;
        return returnParams;
    }

    MyTaskResults *MetaMorpheusTask::RunTask(const std::string &output_folder,
                                             std::vector<DbForTask*> &currentProteinDbFilenameList,
                                             std::vector<std::string> &currentRawDataFilepathList,
                                             const std::string &displayName)
    {
        StartingSingleTask(displayName, privateVerbosityLevel);
        
        GlobalVariables::GlobalVariables_init();

        int rank;
        MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
        
        std::filesystem::path output_directory = output_folder;
        std::string output_dir = output_directory.parent_path().string() + "/Task Settings";
        if ( rank == 0 ) {
            if ( !std::filesystem::exists(output_dir ) ){
                std::filesystem::create_directory(output_dir);
            }       
        
            std::string output_path = output_dir + "/config.toml";
            this->writeTomlConfig(output_path, tomlFile);
            std::vector<std::string> vtd = {displayName};
            FinishedWritingFile(output_path, vtd, privateVerbosityLevel );
        }
        MPI_Barrier ( MPI_COMM_WORLD );
        
        try
        {
            clock_t begin = clock();
            
            std::vector<FileSpecificParameters*> fileSettingsList(currentRawDataFilepathList.size());
            for (int i = 0; i < (int)currentRawDataFilepathList.size(); i++)
            {
                std::string rawFilePath = currentRawDataFilepathList[i];
                std::filesystem::path rdpath = rawFilePath;
                std::string directory = rdpath.parent_path();
                std::string basename = rawFilePath.substr(rawFilePath.find_last_of("/"));
                std::string fname = basename.substr(1, basename.find_last_of(".")-1);
                std::string fileSpecificTomlPath = directory + "/" + fname + ".toml";
                if (std::filesystem::exists(fileSpecificTomlPath))
                {
                    Toml trw;
                    toml::Value toml_value = trw.tomlReadFile(fileSpecificTomlPath);
                    
                    toml::Value* fileParameters = trw.getValue(toml_value, "CommonParameters");
                    toml::Table fileSpecificSettings = fileParameters->as<toml::Table>();              
                    
                    try
                    {
                        fileSettingsList[i] = new FileSpecificParameters(fileSpecificSettings);
                    }
                    catch (const MetaMorpheusException &e)
                    {
                        // file-specific toml has already been validated in the GUI when the spectra files were added, so...
                        // probably the only time you can get here is if the user modifies the file-specific parameter
                        // file in the middle of a run...
                        std::string sw = "Problem parsing the file-specific toml " + fileSpecificTomlPath; 
                        Warn(sw);
                    }
                }
            }
            
            RunSpecific(output_folder, currentProteinDbFilenameList, currentRawDataFilepathList, displayName, fileSettingsList);
            
            clock_t end = clock();
            myTaskResults->Time = (end - begin)/CLOCKS_PER_SEC;
            std::string resultsFileName = output_folder + "/results.txt";

            std::ofstream file(resultsFileName);
            
            file << "Software: " <<  GlobalVariables::getMetaMorpheusVersion() << std::endl;
            file << myTaskResults->ToString() << std::endl;

            std::vector<std::string> svec = {displayName};
            FinishedWritingFile(resultsFileName, svec,  privateVerbosityLevel );
            FinishedSingleTask(displayName,  privateVerbosityLevel);
        }
        catch (const std::runtime_error &e)
        {
            std::string resultsFileName = output_folder + "/results.txt";
            //e.Data->Add("folder", output_folder);
            std::ofstream file(resultsFileName);
            std::string ver = "1.0.0.0";
            if (  GlobalVariables::getMetaMorpheusVersion() == ver ) {
                file << "HPCMorpheus: Not a release version"  << std::endl;
            }
            else {
                file << "HPCMorpheus: version "  <<  GlobalVariables::getMetaMorpheusVersion() << std::endl;
            }
            //file.WriteLine(SystemInfo::CompleteSystemInfo());
            //OS, OS Version, .Net Version, RAM, processor count, MSFileReader .dll versions X3
            //file << "e: " <<  e;
            file << "e.Message: " << e.what();
            //file << "e.InnerException: " <<  e.InnerException;
            //file << "e.Source: " << e.Source;
            //file << "e.StackTrace: " << e.StackTrace;
            //file << "e.TargetSite: " << e.TargetSite;
        } 
       
        
        std::string proseFilePath = output_folder + "/prose.txt";
        std::ofstream file(proseFilePath);
        file << "The data analysis was performed using HPCMorpheus version" <<
            GlobalVariables::getMetaMorpheusVersion() <<
            ", available at " <<   "https://github.com/smith-chem-wisc/MetaMorpheus.";
        file << ProseCreatedWhileRunning->toString();
        //file << SystemInfo::SystemProse()->Replace("\r\n", "") << " ";

        std::string tasktype = "";
        MyTask task = getTaskType();
        if ( task == MyTask::Search ) 
            tasktype = "Search";
        else if ( task == MyTask::Gptmd ) 
            tasktype = "Gptmd";
        else if ( task == MyTask::Calibrate ) 
            tasktype = "Calibrate";
        else if ( task == MyTask::XLSearch ) 
            tasktype == "XLSearch";

        int time_min = (int) (myTaskResults->Time/60);
        int time_sec = myTaskResults->Time - (time_min * 60);
        file << "The total time to perform the " << tasktype << " task on " <<
            std::to_string(currentRawDataFilepathList.size()) <<
            " spectra file(s) was " << time_min << ":" << time_sec <<
            " minutes." << std::endl << std::endl;
        file << "Published works using MetaMorpheus software are encouraged to cite: " <<
            "Solntsev, S. K.; Shortreed, M. R.; Frey, B. L.; Smith, L. M. " <<
            "Enhanced Global Post-translational Modification Discovery with MetaMorpheus. " <<
            "Journal of Proteome Research. 2018, 17 (5), 1844-1851." << std::endl << std::endl;
        
        file <<"Spectra files: " << std::endl;
        
        std::string sjoint1;
        for ( auto b: currentRawDataFilepathList ) {
            sjoint1 += '\t' + b;
        }
        file << sjoint1 << std::endl;

        file << "Databases:" << std::endl;

        std::string sjoint2="";
        for ( auto b: currentProteinDbFilenameList ) {
            if ( b->getIsContaminant() ) {
                sjoint2 += '\t' +  "Contaminant "  + b->getFilePath();
            }
            else {
                sjoint2 += '\t' + b->getFilePath();
            }
        }
        file << sjoint2 << std::endl;
        
        std::vector<std::string> svec2 = {displayName};
        FinishedWritingFile(proseFilePath, svec2,  privateVerbosityLevel);
        
        return myTaskResults;
    }
    
    std::vector<Protein*> MetaMorpheusTask::LoadProteins(const std::string &taskId,
                                                         std::vector<DbForTask*> &dbFilenameList,
                                                         bool searchTarget,
                                                         DecoyType decoyType,
                                                         std::vector<std::string> &localizeableModificationTypes,
                                                         EngineLayer::CommonParameters *commonParameters)
    {
        std::vector<std::string> svec = {taskId};
        Status("Loading proteins...", svec,  privateVerbosityLevel);
        int emptyProteinEntries = 0;
        std::vector<Protein*> proteinList;
        for (auto db : dbFilenameList)
        {
            int emptyProteinEntriesForThisDb = 0;
            std::unordered_map<std::string, Modification*> unknownModifications;
            auto dbProteinList = LoadProteinDb(db->getFilePath(), searchTarget, decoyType,
                                               localizeableModificationTypes, db->getIsContaminant(),
                                               unknownModifications, emptyProteinEntriesForThisDb,
                                               commonParameters);
            proteinList.insert(proteinList.end(), dbProteinList.begin(), dbProteinList.end() );
            emptyProteinEntries += emptyProteinEntriesForThisDb;
        }
        if (proteinList.empty())
        {
            Warn("Warning: No protein entries were found in the database");
        }
        else if (emptyProteinEntries > 0)
        {
            Warn("Warning: " + std::to_string(emptyProteinEntries) + " empty protein entries ignored");
        }
        return proteinList;
    }
    
    std::vector<Protein*> MetaMorpheusTask::LoadProteinDb(const std::string &fileName,
                                                          bool generateTargets,
                                                          DecoyType decoyType,
                                                          std::vector<std::string> &localizeableModificationTypes,
                                                          bool isContaminant,
                                                          std::unordered_map<std::string, Modification*> &um,
                                                          int &emptyEntriesCount,
                                                          EngineLayer::CommonParameters *commonParameters)
    {
        std::vector<std::string> dbErrors;
        std::vector<Protein*> proteinList;
        
        std::string theExtension = fileName.substr(fileName.find_last_of("."));
        std::transform(theExtension.begin(), theExtension.end(), theExtension.begin(), [] (unsigned char c) {
                return std::tolower(c); });
        bool compressed = StringHelper::endsWith(theExtension, "gz");

        // allows for .bgz and .tgz, too which are used on occasion
        std::string fname = fileName.substr(0, fileName.find_last_of("."));
        std::transform(fname.begin(), fname.end(), fname.begin(), [] (unsigned char c) {
                return std::tolower(c); });

        theExtension = compressed ? fname.substr(fname.find_last_of(".")) : theExtension;
        
        if (theExtension == ".fasta" || theExtension == ".fa")
        {
            um.clear();
            proteinList = ProteinDbLoader::LoadProteinFasta(fileName, generateTargets,
                                                            decoyType, isContaminant,
                                                            ProteinDbLoader::UniprotAccessionRegex,
                                                            ProteinDbLoader::UniprotFullNameRegex,
                                                            ProteinDbLoader::UniprotFullNameRegex,
                                                            ProteinDbLoader::UniprotGeneNameRegex,
                                                            ProteinDbLoader::UniprotOrganismRegex,
                                                            dbErrors,
                                                            commonParameters->getMaxThreadsToUsePerFile());
        }
        else
        {
            std::vector<std::string> modTypesToExclude;
            for ( auto b: GlobalVariables::getAllModTypesKnown() ) {
                if (std::find(localizeableModificationTypes.begin(), localizeableModificationTypes.end(), b) ==
                    localizeableModificationTypes.end()) {
                    modTypesToExclude.push_back(b);
                }
            }

            auto tmpmods = GlobalVariables::getAllModsKnown();
            proteinList = ProteinDbLoader::LoadProteinXML(fileName, generateTargets, decoyType,
                                                          tmpmods,
                                                          isContaminant,
                                                          modTypesToExclude, um,
                                                          commonParameters->getMaxThreadsToUsePerFile(),
                                                          commonParameters->getMaxHeterozygousVariants(),
                                                          commonParameters->getMinVariantDepth());
        }

        emptyEntriesCount=0;
        for ( auto p: proteinList ) {
            if ( p->getBaseSequence().length() == 0 ) {
                emptyEntriesCount++;
            }
        }
        
        std::vector<Proteomics::Protein*> tmpvec;
        for ( auto p: proteinList ) {
            if ( p->getBaseSequence().length() > 0 ) {
                tmpvec.push_back(p);
            }
        }

        return tmpvec;
    }
    
    void MetaMorpheusTask::LoadModifications(const std::string &taskId,
                                             std::vector<Modification*> &variableModifications,
                                             std::vector<Modification*> &fixedModifications,
                                             std::vector<std::string> &localizableModificationTypes)
    {
        // load modifications
        Status("Loading modifications...", taskId,  privateVerbosityLevel);

        variableModifications.clear();
        for ( auto b: GlobalVariables::getAllModsKnown() ) {
            auto tmp = getCommonParameters()->getListOfModsVariable();
            std::tuple<std::string, std::string> elem = std::make_tuple( b->getModificationType(),
                                                                         b->getIdWithMotif());
            if ( std::find(tmp->begin(), tmp->end(), elem ) != tmp->end() ) {
                variableModifications.push_back(b);
            }
        }
        
        fixedModifications.clear();
        for ( auto b: GlobalVariables::getAllModsKnown() ) {
            auto tmp = getCommonParameters()->getListOfModsFixed();
            std::tuple<std::string, std::string> elem = std::make_tuple( b->getModificationType(),
                                                                         b->getIdWithMotif());
            if ( std::find(tmp->begin(), tmp->end(), elem ) != tmp->end() ) {
                fixedModifications.push_back(b);
            }
        }
        
        //localizableModificationTypes = GlobalVariables::getAllModTypesKnown().ToList();
        localizableModificationTypes.clear();
        auto tmpmods = GlobalVariables::getAllModTypesKnown();
        for ( auto p = tmpmods.begin(); p != tmpmods.end(); p++ ) {
            localizableModificationTypes.push_back(*p);
        }
        
        std::vector<std::string> recognizedVariable;
        for ( auto p: variableModifications ) {
            recognizedVariable.push_back(p->getIdWithMotif() );
        }

        std::vector<std::string> recognizedFixed;
        for ( auto p: fixedModifications ) {
            recognizedFixed.push_back(p->getIdWithMotif() );
        }

        std::vector<std::string> unknownMods;
        for ( auto p: *(getCommonParameters()->getListOfModsVariable()) ) {
            if ( std::find(recognizedVariable.begin(),recognizedVariable.end(), std::get<1>(p)) ==
                 recognizedVariable.end() ) {
                unknownMods.push_back(std::get<1>(p));
            }
        }
        for ( auto p: *(getCommonParameters()->getListOfModsFixed()) )  {
            if ( std::find(recognizedFixed.begin(),recognizedFixed.end(), std::get<1>(p)) ==
                 recognizedFixed.end() ) {
                unknownMods.push_back(std::get<1>(p));
            }
        }
        
        for (auto unrecognizedMod : unknownMods)
        {
            Warn("Unrecognized mod " + unrecognizedMod + "; are you using an old .toml?");
        }
    }

    void MetaMorpheusTask::WritePsmsToTsv(std::vector<PeptideSpectralMatch*> &psms,
                                          const std::string &filePath,
                                          std::unordered_map<std::string, int> *modstoWritePruned)
    {
        std::ofstream output(filePath);
        output << PeptideSpectralMatch::GetTabSeparatedHeader() << std::endl;
        for (auto psm : psms)
        {
            output << psm->ToString(modstoWritePruned) << std::endl;
        }
    }

    void MetaMorpheusTask::StartingSingleTask(const std::string &taskName, int verbosityLevel)
    {
        if ( verbosityLevel >= 1 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);
            std::cout << "Starting task: " << taskName << std::endl;
        }
        OFFSET_INCR(MetaMorpheus_offset);
    }

    void MetaMorpheusTask::FinishedSingleTask(const std::string &taskName, int verbosityLevel)
    {
        OFFSET_DECR(MetaMorpheus_offset);
        if ( verbosityLevel >= 1 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);                        
            std::cout << "Finished task: " << taskName << std::endl;
        }
    }
    void MetaMorpheusTask::StartingSingleEngine(std::vector<std::string> &nestedIDs, int verbosityLevel)
    {
        if ( verbosityLevel >= 1 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);
            std::cout << "Starting engine: ";
            for ( auto p : nestedIDs ) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
        OFFSET_INCR(MetaMorpheus_offset);
    }

    void MetaMorpheusTask::FinishedSingleEngine(std::vector<std::string> &nestedIDs, MetaMorpheusEngineResults *myResults, int verbosityLevel)
    {

        int time_hour = (int) (myResults->Time/3600);
        int time_min = (int) ((myResults->Time - (time_hour*3600))/60);
        int time_sec = myResults->Time - ((time_min * 60) + (time_hour *3600));
        char timestr[64];
        sprintf(timestr, "%02d:%02d:%02d", time_hour, time_min, time_sec);
        if ( verbosityLevel >= 1 ) {
            std::cout << std::endl;
            DISPLAY_OFFSET(MetaMorpheus_offset);
            std::cout << "Engine results: Time to run : " << timestr << std::endl;
        }

        OFFSET_DECR(MetaMorpheus_offset);
        if ( verbosityLevel >= 1 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);
            std::cout << "Finished engine: ";
            for ( auto p : nestedIDs ) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
        // ToDO print myResults
    }
    
    void MetaMorpheusTask::ReportProgress(ProgressEventArgs *v, int verbosityLevel)
    {
        if ( verbosityLevel >= 2 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);            
            std::cout << v->V << " " ;
            for ( auto p : v->NestedIDs ) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
    }

    void MetaMorpheusTask::ReportEngineProgress(std::string id, int value, int verbosityLevel)
    {
        if ( EngineProgressKeys.find(id) != EngineProgressKeys.end() ) {
            if ( verbosityLevel >= 1 ) {
                std::cout << value << " " << std::flush;
            }
        }
        else {
            EngineProgressKeys.insert(id);
            if ( verbosityLevel >= 1 ) {                
                std::cout << std::endl;
                DISPLAY_OFFSET(MetaMorpheus_offset);            
                std::cout << id << " " << value << " " << std::flush;
            }
        }
    }
    
    
    void MetaMorpheusTask::FinishedWritingFile(const std::string &v, std::vector<std::string> &nestedIDs, int verbosityLevel)
    {
        if ( verbosityLevel >= 1 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);
            std::cout << "FinishedWritingFile: " << v << " " ;
            for ( auto p : nestedIDs ) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
    }
    
    void MetaMorpheusTask::StartingDataFile(const std::string &v, std::vector<std::string> &nestedIDs, int verbosityLevel)
    {
        if ( verbosityLevel >= 4) {
            DISPLAY_OFFSET(MetaMorpheus_offset);
            std::cout << "StartingDataFile " << v << " " ;
            for ( auto p : nestedIDs ) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
    }
    
    void MetaMorpheusTask::FinishedDataFile(const std::string &v, std::vector<std::string> &nestedIDs, int verbosityLevel)
    {
        if ( verbosityLevel >= 4 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);
            std::cout << "FinishedDataFile " << v << " " ;
            for ( auto p : nestedIDs ) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
    }
    
    void MetaMorpheusTask::Status(const std::string &v, const std::string &id, int verbosityLevel)
    {
        if ( verbosityLevel >= 4 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);
            std::cout << v << " " << id << std::endl;
        }
    }
    
    void MetaMorpheusTask::Status(const std::string &v, std::vector<std::string> &nestedIds, int verbosityLevel)
    {
        if ( verbosityLevel >= 4 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);
            std::cout << v << " " ;
            for ( auto p : nestedIds ) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
    }
    
    void MetaMorpheusTask::Warn(const std::string &v)
    {
        std::cout << "Warn : " << v << std::endl;
    }
    
    void MetaMorpheusTask::Log(const std::string &v, std::vector<std::string> &nestedIds, int verbosityLevel)
    {
        if ( verbosityLevel >= 5 ) {
            DISPLAY_OFFSET(MetaMorpheus_offset);  
            std::cout << v << " " ;
            for ( auto p : nestedIds ) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
    }
    
    //void MetaMorpheusTask::NewCollection(const std::string &displayName, std::vector<std::string> &nestedIds)
    //{
    //    StringEventArgs tempVar(displayName, nestedIds);
    //}
    
    std::vector<std::string> MetaMorpheusTask::GetModsTypesFromString(const std::string &value)
    {
        return StringHelper::split(value, '\t');
    }
    
    bool MetaMorpheusTask::SameSettings(const std::string &pathToOldParamsFile, IndexingEngine *indexEngine)
    {
        std::ifstream reader(pathToOldParamsFile);
        std::string content ( (std::istreambuf_iterator<char>(reader) ),
                              (std::istreambuf_iterator<char>() ) );
        if ( content == indexEngine->ToString())
        {
            return true;
        }
        return false;
    }
    
    void MetaMorpheusTask::WritePeptideIndex(std::vector<PeptideWithSetModifications*> &peptideIndex,
                                             std::string &peptideIndexFile)
    {
        PeptideWithSetModifications::Serialize (peptideIndexFile, peptideIndex );
    }
    
    void MetaMorpheusTask::WriteFragmentIndexSerializer(std::vector<std::vector<int>> &fragmentIndex,
                                                        const std::string &fragmentIndexFile)
    {
        std::ofstream output (fragmentIndexFile );
        if ( output.is_open() ) {
            // We are just writing sparse data, since the initial vector can be huge, but most elements are not populated.
            output << fragmentIndex.size() << std::endl;
            for ( unsigned long i = 0;i < fragmentIndex.size() ; i++ ) {
                if ( fragmentIndex[i].size() != 0 ) {
                    output << i << std::endl;
                    for ( auto e : fragmentIndex[i] ) {
                        output << e << " ";
                    }
                    output << std::endl;
                }
            }
            output.close();
        }
        else {
            std::cout << "WriteFragmentIndexSerializer : Could not create file " << fragmentIndexFile << std::endl;
        }
    }

    void MetaMorpheusTask::ReadFragmentIndexDeserializer(std::vector<std::vector<int>> &fragmentIndex,
                                                         const std::string &fragmentIndexFile)
    {
        std::ifstream input (fragmentIndexFile );
        if ( input.is_open() ){
            std::string line;
            std::getline (input, line );
            unsigned long int vecsize = std::stoi ( line );

             try
             {
                 fragmentIndex.resize(vecsize);
             }
             catch (int e1)
             {
                 throw MetaMorpheusException("ReadFragmentIndexDeserializer : Max fragment mass too large for indexing engine");
             }

             while ( std::getline ( input, line ) ) {
                 int index = std::stoi ( line);
                 std::getline (input, line );

                 std::vector<std::string> values = StringHelper::split (line, ' ');
                 for (auto value : values ) {
                     int ivalue = std::stoi ( value );
                     fragmentIndex[index].push_back(ivalue);
                 }
             }
             
            input.close();
        }
        else {
            std::cout << "ReadFragmentIndexDeserializer: Could not open file " << fragmentIndexFile << std::endl;
        }
    }
    
    std::string MetaMorpheusTask::GetExistingFolderWithIndices(IndexingEngine *indexEngine,
                                                               std::vector<DbForTask*> &dbFilenameList)
    {
        for (auto database : dbFilenameList)
        {
            std::filesystem::path dp = database->getFilePath();
            std::string baseDir = dp.parent_path();
            std::string indexDirectory = baseDir + "/" + IndexFolderName;
            
            if (! std::filesystem::exists(indexDirectory) )
            {
                return "";
            }
            
            // all directories in the same directory as the protein database
            std::vector<std::string> directories;
            for ( auto p : std::filesystem::recursive_directory_iterator( indexDirectory) ){
                if ( std::filesystem::is_directory(p) ) {
                    directories.push_back(p.path().string() );
                }
            }
            
            
            // look in each subdirectory to find indexes folder
            for (auto possibleFolder : directories)
            {
                std::string result = CheckFiles(indexEngine, possibleFolder);
                
                if (result != "")
                {
                    return result;
                }
            }
        }
        
        return "";
    }

    std::string MetaMorpheusTask::CheckFiles(IndexingEngine *indexEngine, std::string &folder)
    {
        std::string ret;
        if ( std::filesystem::exists(folder + "/indexEngine.params")    &&
             std::filesystem::exists(folder + "/peptideIndex.ind")      &&
             std::filesystem::exists(folder + "/fragmentIndex.ind")     &&
             (std::filesystem::exists(folder + "/precursorIndex.ind") ||
              !indexEngine->GeneratePrecursorIndex)                                    &&
             SameSettings(folder + "/indexEngine.params", indexEngine) )
        {
            return folder;
        }
        return ret;
    }
    
    void MetaMorpheusTask::WriteIndexEngineParams(IndexingEngine *indexEngine, const std::string &fileName)
    {
        std::ofstream output(fileName);
        output << indexEngine->ToString() ;
        output.close();
    }
    
    std::string MetaMorpheusTask::GenerateOutputFolderForIndices(std::vector<DbForTask*> &dbFilenameList)
    {
        std::filesystem::path dp = dbFilenameList.front()->getFilePath();
        std::string pathToIndexes = dp.parent_path().string() +  "/" + IndexFolderName;
        if (!std::filesystem::exists(pathToIndexes))
        {
            std::filesystem::create_directory(pathToIndexes);
        }
        char dates[100];
        time_t curr_time;
        tm *curr_tm;

        time(&curr_time);
        curr_tm = localtime(&curr_time);
        strftime(dates, 100, "%Y-%m-%d-%H-%M-%S", curr_tm);
        std::string date_string(dates);
        std::string folder = pathToIndexes + "/" + date_string;
        std::filesystem::create_directory(folder);
        return folder;
    }
    
    void MetaMorpheusTask::GenerateIndexes(IndexingEngine *indexEngine,
                                           std::vector<DbForTask*> &dbFilenameList,
                                           std::vector<PeptideWithSetModifications*> &peptideIndex,
                                           std::vector<std::vector<int>> &fragmentIndex,
                                           std::vector<std::vector<int>> &precursorIndex,
                                           std::vector<Protein*> &allKnownProteins,
                                           std::vector<Modification*> &allKnownModifications,
                                           const std::string &taskId,
                                           MPI_Comm comm )
    {
        std::vector<std::string> svec1 = {taskId};        
        Status("Running Index Engine...", svec1,  privateVerbosityLevel);

        int rank;
        MPI_Comm_rank ( comm, &rank );
        if ( rank == 0 ) {
            std::string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);
            
            if (pathToFolderWithIndices == "")
            {
                auto output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
                Status("Writing params...", svec1,  privateVerbosityLevel);
                
                std::string paramsFile = output_folderForIndices + "/indexEngine.params";
                WriteIndexEngineParams(indexEngine, paramsFile);
                FinishedWritingFile(paramsFile, svec1,  privateVerbosityLevel );
                
                Status("Running Index Engine...", svec1,  privateVerbosityLevel);
                auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
                peptideIndex = indexResults->getPeptideIndex();
                fragmentIndex = indexResults->getFragmentIndex();
                precursorIndex = indexResults->getPrecursorIndex();
                
                Status("Writing peptide index...", svec1,  privateVerbosityLevel);
                std::string peptideIndexFile = output_folderForIndices + "/peptideIndex.ind";
                WritePeptideIndex(peptideIndex, peptideIndexFile);
                FinishedWritingFile(peptideIndexFile, svec1,  privateVerbosityLevel);
                
                Status("Writing fragment index...", svec1,  privateVerbosityLevel);
                std::string fragmentIndexFile = output_folderForIndices + "/fragmentIndex.ind";
                WriteFragmentIndexSerializer(fragmentIndex, fragmentIndexFile);
                FinishedWritingFile(fragmentIndexFile, svec1,  privateVerbosityLevel );
                
                if (indexEngine->GeneratePrecursorIndex)
                {
                    Status("Writing precursor index...", svec1,  privateVerbosityLevel );
                    std::string precursorIndexFile = output_folderForIndices + "/precursorIndex.ind";
                    WriteFragmentIndexSerializer(precursorIndex, precursorIndexFile);
                    FinishedWritingFile(precursorIndexFile, svec1,  privateVerbosityLevel );
                }
            }
        }

        // This Barrier is required to ensure that the files are written by rank 0
        // if they didn't already exist
        MPI_Barrier ( comm);
        
        std::string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);
        Status("Reading fragment index...", svec1,  privateVerbosityLevel );
        std::string file = pathToFolderWithIndices + "/fragmentIndex.ind";
        ReadFragmentIndexDeserializer(fragmentIndex, file );
        
        if (indexEngine->GeneratePrecursorIndex)
        {
            Status("Reading precursor index...", svec1,  privateVerbosityLevel );
            file = pathToFolderWithIndices + "/precursorIndex.ind";
            ReadFragmentIndexDeserializer(precursorIndex, file );
        }
        
        Status("Reading peptide index...", svec1,  privateVerbosityLevel);
        file = pathToFolderWithIndices + "/peptideIndex.ind";
        PeptideWithSetModifications::Deserialize(file, peptideIndex );
        
        // populate dictionaries of known proteins for deserialization
        std::unordered_map<std::string, Protein*> proteinDictionary;
        
        for (auto protein : allKnownProteins)
        {
            if (proteinDictionary.find(protein->getAccession()) == proteinDictionary.end())
            {
                proteinDictionary.emplace(protein->getAccession(), protein);
            }
            else if (proteinDictionary[protein->getAccession()]->getBaseSequence() != protein->getBaseSequence())
            {
                throw MetaMorpheusException(StringHelper::formatSimple("The protein database contained multiple proteins with "
                                                                       "accession {0} ! This is not allowed for index-based searches "
                                                                       "(modern, non-specific, crosslink searches)",
                                                                       protein->getAccession()));
            }
        }
        
        // get non-serialized information for the peptides (proteins, mod info)
        auto tmp = GlobalVariables::getAllModsKnownDictionary();
        for (auto peptide : peptideIndex)
        {
            peptide->SetNonSerializedPeptideInfo( tmp, proteinDictionary);
        }    
    }

    // This next routine is just a temporary hack until we find a way to get this information
    // from MsDataFile
    int MetaMorpheusTask::getNumScans ( std::string &filename )
    {
        int size=10;
        if ( filename.find("BSADSSO200-1-08032018_Slot1-12_01_520modified.mgf") != std::string::npos ){
            size = 67711;            
        }
        else if ( filename.find("RibosomeA-10182018_Slot2-01_01_628modified.mgf") != std::string::npos ){
            size = 41642;            
        }
        else if ( filename.find("Shaun-Exp-AP-10122019.mgf") != std::string::npos ){
            size = 23474;            
        }
        else if ( filename.find("B170110_02_Lumos_PR_IN_190_mito-DSS_18A15") != std::string::npos ){
            size = 21614;            
        }
        else if ( filename.find("B170110_03_Lumos_PR_IN_190_mito-DSS_18A16") != std::string::npos ){
            size = 34659;            
        }
        else if ( filename.find("B170110_04_Lumos_PR_IN_190_mito-DSS_18A17") != std::string::npos ){
            size = 32574;            
        }
        else if ( filename.find("B170110_06_Lumos_PR_IN_190_mito-DSS_18A19") != std::string::npos ){
            size = 16169;            
        }
        else {
            std::cout << "MetaMorpheusTask::getNumScans(): file " << filename << " unknown. Returning default value of 10 Spectras\n";
        }
        
        return size;
    }

    void MetaMorpheusTask::DataFilePartitioning ( std::vector<std::string> &allFiles,                   // IN
                                                  MPI_Comm comm,                                        // IN
                                                  std::vector<std::string> &myFile,                     // OUT
                                                  std::vector<std::tuple<int, int>> &myFirstLastIndex ) // OUT

    {
        int rank, numProcs;
        MPI_Comm_rank ( comm, &rank );
        MPI_Comm_size ( comm, &numProcs );
        
        int numFiles = (int) allFiles.size();
        
        if ( numFiles > numProcs ) {
            // currently not handled. Interfaces are however ready to deal with this scenario if
            // we have to later on.
            std::cout << "Number of files exceeds number of process. This scenario is currently " <<
                "not supported. Returning.\n";
            return;
        }
        
        // Generic partitioning is based on the number of Scans per file.   
        // Step 1: determine ideal avg. num. Scans per proc.
        std::vector<int> numScansPerFile (numFiles);
        int totalNumScans = 0;
        for ( int i =0; i < numFiles; i++  ) {
            numScansPerFile[i] = getNumScans ( allFiles[i] );
            totalNumScans     += numScansPerFile[i];
        }
        
        int avgNumScans = std::round(totalNumScans / numProcs );
        
        // Step 2: determine number of Processes assigned to each file.
        //         total sum has to match numProcs
        std::vector<int> avgScansPerFile(numFiles);
        std::vector<int> procsPerFile(numFiles);
        int totalProcs = 0;
        for ( int i=0; i < numFiles; i++ ) {
            procsPerFile[i]    = numScansPerFile[i] / avgNumScans;
            if ( procsPerFile[i] == 0 ) {
                procsPerFile[i] = 1;
            }
            avgScansPerFile[i] = numScansPerFile[i] / procsPerFile[i];
            totalProcs        += procsPerFile[i];
        }
        
        if ( totalProcs < numProcs ) {
            // Assign more procs to the files with largers avg Scans per File
            while ( totalProcs != numProcs ) {
                // Find file with largest avg. and assign an additional process to this file.
                int max_val   = avgScansPerFile[0];
                int max_index = 0;
                for ( int j=1; j< numFiles; j++ ) {
                    if ( avgScansPerFile[j] > max_val ) {
                        max_val   = avgScansPerFile[j];
                        max_index = j;
                    }
                }
                
                procsPerFile[max_index]++;
                avgScansPerFile[max_index] = numScansPerFile[max_index]/procsPerFile[max_index];
                
                totalProcs++;
            }
        }
        else if ( totalProcs > numProcs ) {
            while ( totalProcs != numProcs ) {
                // Find file with smallest avg. and remove a proc from this file.
                // Make sure to not remove a proc from a file that has only one 
                // proc assigned.
                int min_val   = std::numeric_limits<int>::max();
                int min_index = -1;
                
                for ( int j=1; j < numFiles; j++ ) {
                    if ( avgScansPerFile[j] < min_val &&
                         procsPerFile[j] > 1 ) {
                        min_val   = avgScansPerFile[j];
                        min_index = j;
                    }
                }
                
                procsPerFile[min_index]--;
                avgScansPerFile[min_index] = numScansPerFile[min_index]/procsPerFile[min_index];
                
                totalProcs--;           
            }
        }
        
        // Step3: make the assignments for each proc
        //       
        // 3a: which file is assigned to a proc
        std::vector<std::string> filePerProc(numProcs);
        int currentFile=0, procCount=0;
        for ( int i=0; i<numProcs; i++ ) {
            filePerProc[i] = allFiles[currentFile];
            procCount++;
            if ( procCount == procsPerFile[currentFile] ) {
                currentFile++;
                procCount=0;
            }
        }
        
        // 3b: which ranks are assigned to a file
        std::vector<std::vector<int>> ranksPerFile(numFiles);
        int currentRank=0;
        for ( int i=0; i<numFiles; i++ ) {
            for ( int k =0; k < procsPerFile[i]; k++ ) {
                ranksPerFile[i].push_back(currentRank);
                currentRank++;
            }        
        }
        
        // 3c: first and last indices for each rank.
        std::vector<std::tuple<int, int>> indecesPerProc;
        for ( int i =0; i < numFiles; i++ ) {
            int firstIndex=0, lastIndex=0;
            for ( int j=0; j < procsPerFile[i]; j++ ) {
                firstIndex  = lastIndex;
                lastIndex  +=  avgScansPerFile[i];
                if ( j == procsPerFile[i]-1 ) {
                    lastIndex = numScansPerFile[i];
                }
                
                indecesPerProc.push_back(std::make_tuple(firstIndex, lastIndex));
            }
        }
        
#ifdef DEBUG
        //Print out information for debugging
        for (int i=0; i < numProcs; i++ ) {
            std::cout << "Rank: " << i << " file: " << filePerProc[i] << " firstIndex: " <<
                std::get<0>(indecesPerProc[i]) << " lastIndex: " << std::get<1>(indecesPerProc[i])
                      << std::endl;
        }
#endif
        // Step 4: extract information for my rank.
        myFile.push_back(filePerProc[rank]);
        myFirstLastIndex.push_back(indecesPerProc[rank]);
        
        return;
    }
}
