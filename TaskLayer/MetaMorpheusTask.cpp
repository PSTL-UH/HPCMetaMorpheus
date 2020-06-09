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

//For XML serialization / deserialization
#include "include/cereal/types/memory.hpp"
#include "include/cereal/archives/xml.hpp"
#include <include/cereal/types/vector.hpp>

#include <ctime>
#include <experimental/filesystem>
#include <exception>
#include <string>
#include <locale>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <typeinfo>

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{

//The TomlSettings tomlConfig seems to serve as the settings required to read and process or write toml content later on.  
//I'm having trouble determining exactly what settings it is providing though.  It clearly has Tolerance, PpmTolerance, AbsoluteTolerance, and Protease information
//but I havent figured out how this translates to providing read or write settings yet.
#ifdef ORIG
    TomlSettings *const MetaMorpheusTask::tomlConfig = TomlSettings::Create([&] (std::any cfg)   {
            cfg::ConfigureType<Tolerance*>([&] (std::any type) {
                    type::WithConversionFor<TomlString*>([&] (std::any convert) {
                            convert::FromToml([&] (std::any tmlString) {
                                    Tolerance::ParseToleranceString(tmlString->Value);
                                });
                        });
                }).ConfigureType<PpmTolerance*>([&] (std::any type){
                        type::WithConversionFor<TomlString*>([&] (std::any convert){
                                convert::ToToml([&] (std::any custom) {
                                        custom.ToString();
                                    });
                            });
                    }).ConfigureType<AbsoluteTolerance*>([&] (std::any type) {
                            type::WithConversionFor<TomlString*>([&] (std::any convert) {
                                    convert::ToToml([&] (std::any custom){
                                            custom.ToString();
                                        });
                                });
                        }).ConfigureType<Protease*>([&] (std::any type) {
                                type::WithConversionFor<TomlString*>([&] (std::any convert)	{
                                        convert::ToToml([&] (std::any custom) {
                                                custom.ToString();
                                            }).FromToml([&] (std::any tmlString){
                                                    ProteaseDictionary::Dictionary[tmlString->Value];
                                                });
                                    });
                            }).ConfigureType<std::vector<std::string>>([&] (std::any type)	{
                                    type::WithConversionFor<TomlString*>([&] (std::any convert)	{
                                            convert::ToToml([&] (std::any custom) {
                                                    std::string::Join("\t", custom);
                                                }).FromToml([&] (std::any tmlString)	{
                                                        GetModsTypesFromString(tmlString->Value);
                                                    });
                                        });
                                }).ConfigureType<std::vector<(std::string, std::string)*>>([&] (std::any type)	{
                                        type::WithConversionFor<TomlString*>([&] (std::any convert) {
                                                convert::ToToml([&] (std::any custom)	{
                                                        std::string::Join("\t\t", custom->Select([&] (std::any b) {
                                                                    return b::Item1 + "\t" + b::Item2;
                                                                }));
                                                    }).FromToml([&] (std::any tmlString){
                                                            GetModsFromString(tmlString->Value);
                                                        });
                                            });
                                    }); 
        });
#endif

    // toml::Table MetaMorpheusTask::tomlConfig;
    toml::Value MetaMorpheusTask::tomlConfig;

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
    
    void MetaMorpheusTask::setCommonParameters(EngineLayer::CommonParameters *value)
    {
        privateCommonParameters = value;
    }
    
    const std::string MetaMorpheusTask::IndexFolderName = "DatabaseIndex";
    
    std::vector<Ms2ScanWithSpecificMass*> MetaMorpheusTask::GetMs2Scans(MsDataFile *myMSDataFile,
                                                                        const std::string &fullFilePath,
                                                                        EngineLayer::CommonParameters *commonParameters)
    {
#ifdef ORIG
        auto ms2Scans = myMSDataFile->GetAllScansList().Where([&] (std::any x)		{
                return x::MsnOrder > 1;
            })->ToArray();
#endif
        std::vector<MsDataScan*> ms2Scans;
        for ( auto x : myMSDataFile->GetAllScansList() ) {
            if ( x->getMsnOrder() > 1 ) {
                ms2Scans.push_back(x);
            }
        }
                

        std::vector<std::vector<Ms2ScanWithSpecificMass*>> scansWithPrecursors(ms2Scans.size());
        
#ifdef ORIG
        //ParallelOptions *tempVar = new ParallelOptions();
        //tempVar->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
        //Parallel::ForEach(Partitioner::Create(0, ms2Scans.size()), tempVar, [&] (partitionRange, loopState)       {
        //        for (int i = partitionRange::Item1; i < partitionRange::Item2; i++)
#endif
        for ( int i = 0; i < (int)ms2Scans.size(); i++ ) {
            if (GlobalVariables::getStopLoops())
            {
                break;
            }
            
            MsDataScan *ms2scan = ms2Scans[i];
            
            std::vector<std::pair<double, int>> precursors;
            if (ms2scan->getOneBasedPrecursorScanNumber().has_value())
            {
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
                    //Warn(s);
                    std::cout << s << std::endl;
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
            }
            
            if (commonParameters->getUseProvidedPrecursorInfo()      &&
                ms2scan->getSelectedIonChargeStateGuess().has_value() )
            {
                auto precursorCharge = ms2scan->getSelectedIonChargeStateGuess().value();
                if (ms2scan->getSelectedIonMonoisotopicGuessMz().has_value())
                {
                    double precursorMZ = ms2scan->getSelectedIonMonoisotopicGuessMz().value();
#ifdef ORIG
                    //if (!precursors.Any([&] (std::any b)   {
                    //            commonParameters->getDeconvolutionMassTolerance()->Within(precursorMZ.ToMass(precursorCharge),
                    //                                                                      b::Item1->ToMass(b::Item2));
                    //        }))
#endif
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
#ifdef ORIG
                    if (!precursors.Any([&] (std::any b){
                                commonParameters->getDeconvolutionMassTolerance()->Within(precursorMZ.ToMass(precursorCharge),
                                                                                          b::Item1->ToMass(b::Item2));
                            }))
#endif
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
            }
            
            scansWithPrecursors[i] = std::vector<Ms2ScanWithSpecificMass*>();
            std::vector<IsotopicEnvelope*> neutralExperimentalFragments = Ms2ScanWithSpecificMass::GetNeutralExperimentalFragments(ms2scan, commonParameters);
            
            for (auto precursor : precursors)
            {
                auto  tempVar2 = new Ms2ScanWithSpecificMass(ms2scan, precursor.first, precursor.second, fullFilePath,
                                                             commonParameters, neutralExperimentalFragments);
                scansWithPrecursors[i].push_back(tempVar2);
            }
        } 
    

#ifdef ORIG
        return scansWithPrecursors.SelectMany([&] (std::any p)	{
                return p;
            });
#endif
        std::vector<Ms2ScanWithSpecificMass*> tmpv;
        for ( auto p: scansWithPrecursors ) {
            for ( auto v: p ) {
                tmpv.push_back(v);
            }
        }
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

        //C# TO C++ CONVERTER TODO TASK: A 'delete returnParams' statement was not added since returnParams
        //was used in a 'return' or 'throw' statement.
        delete fileSpecificDigestionParams;
        return returnParams;
    }

    MyTaskResults *MetaMorpheusTask::RunTask(const std::string &output_folder,
                                             std::vector<DbForTask*> &currentProteinDbFilenameList,
                                             std::vector<std::string> &currentRawDataFilepathList,
                                             const std::string &displayName)
    {
        StartingSingleTask(displayName);
        
#ifdef ORIG
        auto tomlFileName = FileSystem::combine(Directory::GetParent(output_folder)->ToString(), "Task Settings", displayName + "config.toml");
        Toml::WriteFile(this, tomlFileName, tomlConfig);
        FinishedWritingFile(tomlFileName, std::vector<std::string> {displayName});
#endif

        std::experimental::filesystem::path output_directory = output_folder;
        std::string output_path = output_directory.parent_path().string() + "Task Settings" + displayName + "config.toml";
        Toml trw;
        trw.tomlWriteNewFile(output_path, tomlConfig);

        
        //MetaMorpheusEngine::FinishedSingleEngineHandler->addListener("SingleEngineHandlerInTask",
        //                                                             SingleEngineHandlerInTask);
        try
        {
            clock_t begin = clock();
            
            std::vector<FileSpecificParameters*> fileSettingsList(currentRawDataFilepathList.size());
            for (int i = 0; i < (int)currentRawDataFilepathList.size(); i++)
            {
                if (GlobalVariables::getStopLoops())
                {
                    break;
                }
                std::string rawFilePath = currentRawDataFilepathList[i];
                std::experimental::filesystem::path rdpath = rawFilePath;
                std::string directory = rdpath.parent_path();
                std::string fname = rawFilePath.substr(rawFilePath.find_last_of("."));
                std::string fileSpecificTomlPath = directory + "/" + fname + ".toml";
                if (std::experimental::filesystem::exists(fileSpecificTomlPath))
                {
#ifdef ORIG
                    //In the Nett package the second parameter for ReadFile are the settings used to process the toml content
                    //public static T ReadFile<T>(string filePath, TomlSettings settings)
                    TomlTable *fileSpecificSettings = Toml::ReadFile(fileSpecificTomlPath, tomlConfig);
#endif
                    toml::Value toml_value = trw.tomlReadFile(fileSpecificTomlPath);
                    
                    //original
                    // toml::Value* fileParameters = trw.getValue(toml_value, tomlConfig);

                    //need string of header we are looking for in toml file.  is it "CommonParameters"?
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
                        //Warn("Problem parsing the file-specific toml " + FileSystem::getFileName(fileSpecificTomlPath) + "; "
                        //     + e->what() +  "; is the toml from an older version of MetaMorpheus?");
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
            
            file << "MetaMorpheus: version " <<  GlobalVariables::getMetaMorpheusVersion() << std::endl;
            file << myTaskResults->ToString() << std::endl;

            std::vector<std::string> svec = {displayName};
            FinishedWritingFile(resultsFileName, svec );
            FinishedSingleTask(displayName);
        }
        catch (const std::runtime_error &e)
        {
            //MetaMorpheusEngine::FinishedSingleEngineHandler->removeListener("SingleEngineHandlerInTask");
            std::string resultsFileName = output_folder + "/results.txt";
            //e.Data->Add("folder", output_folder);
            std::ofstream file(resultsFileName);
            std::string ver = "1.0.0.0";
            if (  GlobalVariables::getMetaMorpheusVersion() == ver ) {
                file << "MetaMorpheus: Not a release version"  << std::endl;
            }
            else {
                file << "MetaMorpheus: version "  <<  GlobalVariables::getMetaMorpheusVersion() << std::endl;
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
        file << "The data analysis was performed using MetaMorpheus version " <<
            GlobalVariables::getMetaMorpheusVersion() <<
            ", available at " <<   "https://github.com/smith-chem-wisc/MetaMorpheus.";
        file << ProseCreatedWhileRunning->toString();
        //file << SystemInfo::SystemProse()->Replace("\r\n", "") << " ";
        file << "The total time to perform the " << static_cast<int>(getTaskType()) << " task on " <<
            std::to_string(currentRawDataFilepathList.size()) <<
            " spectra file(s) was " << myTaskResults->Time <<
            " minutes." << std::endl;
        file << "\n Published works using MetaMorpheus software are encouraged to cite: " <<
            "Solntsev, S. K.; Shortreed, M. R.; Frey, B. L.; Smith, L. M. " <<
            "Enhanced Global Post-translational Modification Discovery with MetaMorpheus. " <<
            "Journal of Proteome Research. 2018, 17 (5), 1844-1851." << std::endl;
        
        file <<"\n Spectra files: " << std::endl;
#ifdef ORIG
        //file << std::string::Join("\r\n", currentRawDataFilepathList.Select([&] (std::any b)    {
        //            return '\t' + b;
        //            }))) << std::endl;
#endif
        
        std::string sjoint1;
        for ( auto b: currentRawDataFilepathList ) {
            sjoint1 += '\t' + b;
        }
        file << sjoint1 << std::endl;

        file << "Databases:" << std::endl;
#ifdef ORIG    
        //file << std::string::Join("\r\n", currentProteinDbFilenameList.Select([&] (std::any b)		{
        //        return '\t' + (b::IsContaminant ? "Contaminant " : "") + b::FilePath;
        //        })) << std::endl;
#endif
        std::string sjoint2;
        for ( auto b: currentProteinDbFilenameList ) {
            sjoint2 += '\t' + (b->getIsContaminant() ? "Contaminant " : "") + b->getFilePath();
        }
        std::vector<std::string> svec2 = {displayName};
        FinishedWritingFile(proseFilePath, svec2);
        
        //MetaMorpheusEngine::FinishedSingleEngineHandler->removeListener("SingleEngineHandlerInTask");
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
        Status("Loading proteins...", svec);
        int emptyProteinEntries = 0;
        std::vector<Protein*> proteinList;
        for (auto db : dbFilenameList)
        {
            int emptyProteinEntriesForThisDb = 0;
            std::unordered_map<std::string, Modification*> unknownModifications;
            proteinList = LoadProteinDb(db->getFilePath(), searchTarget, decoyType,
                                        localizeableModificationTypes, db->getIsContaminant(),
                                        unknownModifications, emptyProteinEntriesForThisDb,
                                        commonParameters);
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
#ifdef ORIG
            std::vector<std::string> modTypesToExclude = GlobalVariables::getAllModTypesKnown().Where([&] (
                                                                                             std::any b)  {
                    std::find(localizeableModificationTypes.begin(), localizeableModificationTypes.end(), b) ==
                    localizeableModificationTypes.end();
                }).ToList();
#endif
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

#ifdef ORIG
        emptyEntriesCount = proteinList.size()([&] (std::any p)	{
                return p::BaseSequence->Length == 0;
            });
#endif
        emptyEntriesCount=0;
        for ( auto p: proteinList ) {
            if ( p->getBaseSequence().length() == 0 ) {
                emptyEntriesCount++;
            }
        }
        
#ifdef ORIG
        return proteinList.Where([&] (std::any p){
                return p::BaseSequence->Length > 0;
            }).ToList();
#endif
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
        Status("Loading modifications...", taskId);
#ifdef ORIG
        variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)    {
                getCommonParameters()->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
            }).ToList();
#endif
        variableModifications.clear();
        for ( auto b: GlobalVariables::getAllModsKnown() ) {
            auto tmp = getCommonParameters()->getListOfModsVariable();
            std::tuple<std::string, std::string> elem = std::make_tuple( b->getModificationType(),
                                                                         b->getIdWithMotif());
            if ( std::find(tmp->begin(), tmp->end(), elem ) != tmp->end() ) {
                variableModifications.push_back(b);
            }
        }
        
#ifdef ORIG
        fixedModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)	{
                getCommonParameters()->ListOfModsFixed->Contains((b::ModificationType, b::IdWithMotif));
            }).ToList();
#endif
        fixedModifications.clear();
        for ( auto b: GlobalVariables::getAllModsKnown() ) {
            auto tmp = getCommonParameters()->getListOfModsFixed();
            std::tuple<std::string, std::string> elem = std::make_tuple( b->getModificationType(),
                                                                         b->getIdWithMotif());
            if ( std::find(tmp->begin(), tmp->end(), elem ) != tmp->end() ) {
                variableModifications.push_back(b);
            }
        }
        
        //localizableModificationTypes = GlobalVariables::getAllModTypesKnown().ToList();
        localizableModificationTypes.clear();
        auto tmpmods = GlobalVariables::getAllModTypesKnown();
        for ( auto p = tmpmods.begin(); p != tmpmods.end(); p++ ) {
            localizableModificationTypes.push_back(*p);
        }
        
#ifdef ORIG
        auto recognizedVariable = variableModifications.Select([&] (std::any p)      {
                p::IdWithMotif;
            });
#endif
        std::vector<std::string> recognizedVariable;
        for ( auto p: variableModifications ) {
            recognizedVariable.push_back(p->getIdWithMotif() );
        }

#ifdef ORIG
        auto recognizedFixed = fixedModifications.Select([&] (std::any p){
                p::IdWithMotif;
            });
#endif
        std::vector<std::string> recognizedFixed;
        for ( auto p: fixedModifications ) {
            recognizedFixed.push_back(p->getIdWithMotif() );
        }

#ifdef ORIG
        auto unknownMods = getCommonParameters()->ListOfModsVariable->Select([&] (std::any p){
                p::Item2;
            }).Except(recognizedVariable).ToList();
        unknownMods.AddRange(getCommonParameters()->ListOfModsFixed->Select([&] (std::any p){
                    p::Item2;
		}).Except(recognizedFixed));
#endif
        std::vector<std::string> unknownMods;
        for ( auto p: *(getCommonParameters()->getListOfModsVariable()) ) {
            if ( std::find(recognizedVariable.begin(),recognizedVariable.end(), std::get<1>(p)) !=
                 recognizedVariable.end() ) {
                unknownMods.push_back(std::get<1>(p));
            }
        }
        for ( auto p: *(getCommonParameters()->getListOfModsFixed()) )  {
            if ( std::find(recognizedFixed.begin(),recognizedFixed.end(), std::get<1>(p)) !=
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
    
    void MetaMorpheusTask::ReportProgress(ProgressEventArgs *v)
    {
        if ( OutProgressHandler != nullptr )
            OutProgressHandler->Invoke(*v);
    }
    
    void MetaMorpheusTask::FinishedWritingFile(const std::string &path, std::vector<std::string> &nestedIDs)
    {
        SingleFileEventArgs tempVar(path, nestedIDs);
        if ( FinishedWritingFileHandler != nullptr )
            FinishedWritingFileHandler->Invoke(tempVar);
    }
    
    void MetaMorpheusTask::StartingDataFile(const std::string &v, std::vector<std::string> &nestedIDs)
    {
        StringEventArgs tempVar(v, nestedIDs);
        if ( StartingDataFileHandler != nullptr )
            StartingDataFileHandler->Invoke(tempVar);
    }
    
    void MetaMorpheusTask::FinishedDataFile(const std::string &v, std::vector<std::string> &nestedIDs)
    {
        StringEventArgs tempVar(v, nestedIDs);
        if ( FinishedDataFileHandler != nullptr )
            FinishedDataFileHandler->Invoke(tempVar);
    }
    
    void MetaMorpheusTask::Status(const std::string &v, const std::string &id)
    {
        std::vector<std::string> s = {id};
        StringEventArgs tempVar(v, s);
        if ( OutLabelStatusHandler != nullptr )
            OutLabelStatusHandler->Invoke(tempVar);
    }
    
    void MetaMorpheusTask::Status(const std::string &v, std::vector<std::string> &nestedIds)
    {
        StringEventArgs tempVar(v, nestedIds);
        if ( OutLabelStatusHandler != nullptr )
            OutLabelStatusHandler->Invoke(tempVar);
    }
    
    void MetaMorpheusTask::Warn(const std::string &v)
    {
        StringEventArgs tempVar(v, std::vector<std::string>() );
        if ( WarnHandler != nullptr )
            WarnHandler->Invoke(tempVar);
    }
    
    void MetaMorpheusTask::Log(const std::string &v, std::vector<std::string> &nestedIds)
    {
        StringEventArgs tempVar(v, nestedIds);
        if ( LogHandler != nullptr )
            LogHandler->Invoke(tempVar);
    }
    
    void MetaMorpheusTask::NewCollection(const std::string &displayName, std::vector<std::string> &nestedIds)
    {
        StringEventArgs tempVar(displayName, nestedIds);
        if ( NewCollectionHandler != nullptr )
            NewCollectionHandler->Invoke(tempVar);
    }
    
    std::vector<std::string> MetaMorpheusTask::GetModsTypesFromString(const std::string &value)
    {
        //return value.Split({"\t"}, StringSplitOptions::RemoveEmptyEntries).ToList();
        return StringHelper::split(value, '\t');
    }
    
    //void MetaMorpheusTask::SingleEngineHandlerInTask(std::any sender, SingleEngineFinishedEventArgs *e)
    void MetaMorpheusTask::SingleEngineHandlerInTask(SingleEngineFinishedEventArgs e)
    {
        myTaskResults->AddResultText(e.ToString());
    }
    
    void MetaMorpheusTask::FinishedSingleTask(const std::string &displayName)
    {
        SingleTaskEventArgs tempVar(displayName);
        if ( FinishedSingleTaskHandler != nullptr )
            FinishedSingleTaskHandler->Invoke(tempVar);
    }
    
    void MetaMorpheusTask::StartingSingleTask(const std::string &displayName)
    {
        SingleTaskEventArgs tempVar(displayName);
        if ( StartingSingleTaskHandler != nullptr )
            StartingSingleTaskHandler->Invoke(tempVar);
    }

#ifdef NOT_NOW
    std::vector<std::type_info> MetaMorpheusTask::GetSubclassesAndItself(std::type_info type)
    {
        std::vector<type> tmp;
        return tmp;
    }
#endif
    
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
                                             const std::string &peptideIndexFile)
    {
        //auto messageTypes = GetSubclassesAndItself(std::vector<PeptideWithSetModifications*>::typeid);
        //auto messageTypes = GetSubclassesAndItself(typeid(std::vector<PeptideWithSetModifications*>));

        //original
        // auto ser = new NetSerializer::Serializer(messageTypes);        
        // auto file = File::Create(peptideIndexFile);
        // ser->Serialize(file, peptideIndex);

        std::ofstream os(peptideIndexFile);
        cereal::XMLOutputArchive archive( os );

        //Cereal only works with smart pointers, so all raw pointers in
        //std::vector<PeptideWithSetModifications*> peptideIndex
        //must be converted to unique pointers with std::make_unique.  This requires traversing the
        //vector and performing the
        //transformation for each entry.
        std::vector<std::unique_ptr<PeptideWithSetModifications>> unique_vector;

        for (int i=0; i < (int)peptideIndex.size(); i++){

            //create unique pointer from raw pointer
            auto p = std::make_unique<PeptideWithSetModifications>(*peptideIndex[i]);
                
            //push unique pointer into vector of unique pointers
            unique_vector.push_back(std::move(p));
        }

        //finally serialize the vector of unique pointers to an XML file
        //this will write the info contained in peptideIndex to the XML file 
        //specified in the peptideIndexFile variable.
#ifdef DO_A_BIT_LATER
        archive(unique_vector);
#else
#pragma message("Warning: a CEREAL statement in Line 876 in MetaMorpheusTask.cpp still needs fixing.")
#endif
        //----------------------------------------------------------------
    }
    
    void MetaMorpheusTask::WriteFragmentIndexNetSerializer(std::vector<std::vector<int>> &fragmentIndex,
                                                           const std::string &fragmentIndexFile)
    {
        //auto messageTypes = GetSubclassesAndItself(std::vector<std::vector<int>>::typeid);
        
        //original
        // auto ser = new NetSerializer::Serializer(messageTypes);        
        // auto file = File::Create(fragmentIndexFile);
        // ser->Serialize(file, fragmentIndex);

        std::ofstream os(fragmentIndexFile);
        cereal::XMLOutputArchive archive( os );

        //Cereal only works with smart pointers, so all raw pointers in
        //std::vector<std::vector<int>> fragmentIndex
        //must be converted to unique pointers with std::make_unique.  This requires traversing the
        //vector and performing the
        //transformation for each entry.
        std::vector<std::vector<std::unique_ptr<int>>> unique_vector;

        for (int i=0; i < (int)fragmentIndex.size(); i++){
            std::vector<std::unique_ptr<int>> unique_vec_row;
            for (int j=0; j < (int)fragmentIndex[i].size(); j++){
                //create unique pointer from raw pointer
                auto p = std::make_unique<int>(fragmentIndex[i][j]);
                    
                //push unique pointer into vector of unique pointers
                unique_vec_row.push_back(std::move(p));
            }
            unique_vector.push_back(std::move(unique_vec_row));
        }

        //finally serialize the vector of vectors of unique pointers to an XML file
        //this will write the info contained in fragmentIndex to the XML file 
        //specified in the indexPath variable.
        archive(fragmentIndex);
        //----------------------------------------------------------------
    }
    
    std::string MetaMorpheusTask::GetExistingFolderWithIndices(IndexingEngine *indexEngine,
                                                               std::vector<DbForTask*> &dbFilenameList)
    {
        for (auto database : dbFilenameList)
        {
            std::experimental::filesystem::path dp = database->getFilePath();
            std::string baseDir = dp.parent_path();
            std::string indexDirectory = baseDir + "/" + IndexFolderName;
            
            if (! std::experimental::filesystem::exists(indexDirectory) )
            {
                return "";
            }
            
            // all directories in the same directory as the protein database
            std::vector<std::string> directories;
            for ( auto p : std::experimental::filesystem::recursive_directory_iterator( indexDirectory) ){
                if ( std::experimental::filesystem::is_directory(p) ) {
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
        if ( std::experimental::filesystem::exists(folder + "/indexEngine.params")    &&
             std::experimental::filesystem::exists(folder + "/peptideIndex.ind")      &&
             std::experimental::filesystem::exists(folder + "/fragmentIndex.ind")     &&
             (std::experimental::filesystem::exists(folder + "/precursorIndex.ind") ||
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
        output << indexEngine ;
    }
    
    std::string MetaMorpheusTask::GenerateOutputFolderForIndices(std::vector<DbForTask*> &dbFilenameList)
    {
        std::string pathToIndexes = dbFilenameList.front()->getFilePath() +  IndexFolderName;
        if (!std::experimental::filesystem::exists(pathToIndexes))
        {
            std::experimental::filesystem::create_directory(pathToIndexes);
        }
        char dates[100];
        time_t curr_time;
        tm *curr_tm;

        time(&curr_time);
        curr_tm = localtime(&curr_time);
        strftime(dates, 100, "%Y-%m-%d-%H-%M-%S", curr_tm);
        std::string date_string(dates);
        std::string folder = pathToIndexes + date_string;
        std::experimental::filesystem::create_directory(folder);
        return folder;
    }
    
    void MetaMorpheusTask::GenerateIndexes(IndexingEngine *indexEngine,
                                           std::vector<DbForTask*> &dbFilenameList,
                                           std::vector<PeptideWithSetModifications*> &peptideIndex,
                                           std::vector<std::vector<int>> &fragmentIndex,
                                           std::vector<std::vector<int>> &precursorIndex,
                                           std::vector<Protein*> &allKnownProteins,
                                           std::vector<Modification*> &allKnownModifications,
                                           const std::string &taskId)
    {
        std::string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);
        std::vector<std::string> svec1 = {taskId};
            
        if (pathToFolderWithIndices == "")
        {
            auto output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
            Status("Writing params...", svec1);

            std::string paramsFile = output_folderForIndices + "/indexEngine.params";
            WriteIndexEngineParams(indexEngine, paramsFile);
            FinishedWritingFile(paramsFile, svec1 );
            
            Status("Running Index Engine...", svec1);
            auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
            peptideIndex = indexResults->getPeptideIndex();
            fragmentIndex = indexResults->getFragmentIndex();
            precursorIndex = indexResults->getPrecursorIndex();
            
            Status("Writing peptide index...", svec1);
            std::string peptideIndexFile = output_folderForIndices + "/peptideIndex.ind";
            WritePeptideIndex(peptideIndex, peptideIndexFile);
            FinishedWritingFile(peptideIndexFile, svec1);
            
            Status("Writing fragment index...", svec1);
            std::string fragmentIndexFile = output_folderForIndices + "/fragmentIndex.ind";
            WriteFragmentIndexNetSerializer(fragmentIndex, fragmentIndexFile);
            FinishedWritingFile(fragmentIndexFile, svec1 );
            
            if (indexEngine->GeneratePrecursorIndex)
            {
                Status("Writing precursor index...", svec1 );
                std::string precursorIndexFile = output_folderForIndices + "/precursorIndex.ind";
                WriteFragmentIndexNetSerializer(precursorIndex, precursorIndexFile);
                FinishedWritingFile(precursorIndexFile, svec1 );
            }
        }
        else
        {
            Status("Reading peptide index...", svec1 );
            //auto messageTypes = GetSubclassesAndItself(std::vector<PeptideWithSetModifications*>::typeid);


            // auto ser = new NetSerializer::Serializer(messageTypes);
            // auto file = File::OpenRead(pathToFolderWithIndices + "/peptideIndex.ind");
            // peptideIndex = static_cast<std::vector<PeptideWithSetModifications*>>(ser->Deserialize(file));

            std::string file = pathToFolderWithIndices + "/peptideIndex.ind";

            //------------------------------------------------------
            //DESERIALIZE FILE INDEXFILE TO OBJECT

            //open ifstream for indexPath
            std::ifstream is(file.c_str());
            cereal::XMLInputArchive archive_read(is);

            //Cereal only has functionality for smart pointers.  Must deserialize data from file to vector
            //of vecotrs of unique pointers
            std::vector<std::unique_ptr<PeptideWithSetModifications>> unique_vec;

            //deserialize xml file to unique_vec structure
            archive_read(unique_vec);

            //now must convert unique pointers to raw pointers.  This requires traversing the structure
            std::vector<PeptideWithSetModifications*> raw_vec;

            for (int i=0; i < (int)unique_vec.size(); i++){

                std::unordered_map<std::string, Modification*> allModsOneIsNterminus;
                PeptideWithSetModifications* p = new PeptideWithSetModifications(unique_vec[i]->getFullSequence(), allModsOneIsNterminus);
                raw_vec.push_back(p);
            }

            //set peptideIndex to equal the vector of vectors of raw pointers.
            peptideIndex = raw_vec;
            //------------------------------------------------------
            is.close();


            
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
                    throw MetaMorpheusException(StringHelper::formatSimple("The protein database contained multiple proteins with accession {0} ! This is not allowed for index-based searches (modern, non-specific, crosslink searches)", protein->getAccession()));
                }
            }
            
            // get non-serialized information for the peptides (proteins, mod info)
            for (auto peptide : peptideIndex)
            {
                auto tmp = GlobalVariables::getAllModsKnownDictionary();
                peptide->SetNonSerializedPeptideInfo( tmp, proteinDictionary);
            }
            
            Status("Reading fragment index...", svec1 );
            //messageTypes = GetSubclassesAndItself(std::vector<std::vector<int>>::typeid);

            // ser = new NetSerializer::Serializer(messageTypes);
            // auto file = File::OpenRead(pathToFolderWithIndices + "/fragmentIndex.ind");
            // fragmentIndex = static_cast<std::vector<std::vector<int>>>(ser->Deserialize(file));


            file = pathToFolderWithIndices + "/fragmentIndex.ind";
            //------------------------------------------------------
            //DESERIALIZE FILE INDEXFILE TO OBJECT

            //open ifstream for indexPath
            std::ifstream is2(file.c_str());
            cereal::XMLInputArchive archive_read_2(is2);

            //Cereal only has functionality for smart pointers.  Must deserialize data from file to vector
            //of vecotrs of unique pointers
            std::vector<std::vector<std::unique_ptr<int>>> unique_vector;

            //deserialize xml file to unique_vec structure
            archive_read_2(unique_vector);

            //now must convert unique pointers to raw pointers.  This requires traversing the structure
            std::vector<std::vector<int>> raw_vector;

            for (int i=0; i < (int)unique_vector.size(); i++){
                std::vector<int> raw_vec_row;
                
                for (int j=0;j< (int)unique_vector[i].size();j++){
                    int* p = unique_vector[i][j].get();
                    raw_vec_row.push_back(*p);
                }

                raw_vector.push_back(raw_vec_row);
            }

            //set _indexedpeaks to equal the vector of vectors of raw pointers.
            fragmentIndex = raw_vector;
            // return raw_vec;
            //------------------------------------------------------
            is2.close();            
            
            if (indexEngine->GeneratePrecursorIndex)
            {
                Status("Reading precursor index...", svec1 );
                //messageTypes = GetSubclassesAndItself(std::vector<std::vector<int>>::typeid);


                // ser = new NetSerializer::Serializer(messageTypes);
                // auto file = File::OpenRead(pathToFolderWithIndices + "/precursorIndex.ind");
                // precursorIndex = static_cast<std::vector<std::vector<int>>>(ser->Deserialize(file));

                file = pathToFolderWithIndices + "/precursorIndex.ind";
                //------------------------------------------------------
                //DESERIALIZE FILE INDEXFILE TO OBJECT

                //open ifstream for indexPath
                std::ifstream is3(file.c_str());
                cereal::XMLInputArchive archive_read_3(is3);

                //Cereal only has functionality for smart pointers.  Must deserialize data from file to vector
                //of vecotrs of unique pointers
                std::vector<std::vector<std::unique_ptr<int>>> unique_vector2;

                //deserialize xml file to unique_vec structure
                archive_read_3(unique_vector2);

                //now must convert unique pointers to raw pointers.  This requires traversing the structure
                std::vector<std::vector<int>> raw_vector2;

                for (int i=0; i < (int)unique_vector2.size(); i++){
                    std::vector<int> raw_vector_row;
                    
                    for (int j=0;j< (int)unique_vector2[i].size();j++){
                        int* p = unique_vector2[i][j].get();
                        raw_vector_row.push_back(*p);
                    }

                    raw_vector2.push_back(raw_vector_row);
                }

                //set _indexedpeaks to equal the vector of vectors of raw pointers.
                precursorIndex = raw_vector2;
                // return raw_vec;
                //------------------------------------------------------
                is3.close(); 
            }
        }
    }
}
