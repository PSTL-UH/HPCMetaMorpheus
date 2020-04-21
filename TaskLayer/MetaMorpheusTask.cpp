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


#include "Chemistry/ClassExtensions.h"
#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"

#include <ctime>
#include <experimental/filesystem>
#include <exception>
#include <string>
#include <locale>
#include <algorithm>

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
                Ms2ScanWithSpecificMass tempVar2(ms2scan, precursor.first, precursor.second, fullFilePath,
                                                 commonParameters, neutralExperimentalFragments);
                scansWithPrecursors[i].push_back(&tempVar2);
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
        //Dr. Gabriel, maybe we can save tomlConfig as a tomlTable rather than the TomlSettings type used above?
        //This might simplify using it here.
        std::experimental::filesystem::path output_directory = output_folder;
        std::string output_path = output_directory.parent_path().string() + "Task Settings" + displayName + "config.toml";
        Toml trw;
        trw.tomlWriteNewFile(tomlFileName, tomlConfig);

        
        MetaMorpheusEngine::FinishedSingleEngineHandler->addListener("SingleEngineHandlerInTask", [&] (SingleEngineFinishedEventArgs* e) {
                SingleEngineHandlerInTask( e);});
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
                    
                    //What is the header or key we are looking for here in the toml file?
                    //In a couple of config files I was looking at it was "CommonParameters", I'm not 
                    //entirely sure that's correct here though.
                    toml::Value* fileParameters = trw.getValue(toml_value, tomlConfig);
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
            MetaMorpheusEngine::FinishedSingleEngineHandler->removeListener("SingleEngineHandlerInTask");
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
        
        MetaMorpheusEngine::FinishedSingleEngineHandler->removeListener("SingleEngineHandlerInTask");
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
            
            proteinList = ProteinDbLoader::LoadProteinXML(fileName, generateTargets, decoyType,
                                                          GlobalVariables::getAllModsKnown(), isContaminant,
                                                          modTypesToExclude, um,
                                                          commonParameters->getMaxThreadsToUsePerFile(),
                                                          commonParameters->getMaxHeterozygousVariants(),
                                                          commonParameters->getMinVariantDepth());
        }

        emptyEntriesCount = proteinList.size()([&] (std::any p)	{
                return p::BaseSequence->Length == 0;
            });

        return proteinList.Where([&] (std::any p){
                return p::BaseSequence->Length > 0;
            }).ToList();
    }
    
    void MetaMorpheusTask::LoadModifications(const std::string &taskId, std::vector<Modification*> &variableModifications,
                                             std::vector<Modification*> &fixedModifications,
                                             std::vector<std::string> &localizableModificationTypes)
    {
        // load modifications
        Status("Loading modifications...", taskId);
        variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)    {
                getCommonParameters()->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
            }).ToList();
        fixedModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)	{
                getCommonParameters()->ListOfModsFixed->Contains((b::ModificationType, b::IdWithMotif));
            }).ToList();
        localizableModificationTypes = GlobalVariables::getAllModTypesKnown().ToList();
        
        auto recognizedVariable = variableModifications.Select([&] (std::any p)      {
                p::IdWithMotif;
            });
        auto recognizedFixed = fixedModifications.Select([&] (std::any p){
                p::IdWithMotif;
            });
        auto unknownMods = getCommonParameters()->ListOfModsVariable->Select([&] (std::any p){
                p::Item2;
            }).Except(recognizedVariable).ToList();
        unknownMods.AddRange(getCommonParameters()->ListOfModsFixed->Select([&] (std::any p){
                    p::Item2;
		}).Except(recognizedFixed));
        for (auto unrecognizedMod : unknownMods)
        {
            Warn("Unrecognized mod " + unrecognizedMod + "; are you using an old .toml?");
        }
    }

    void MetaMorpheusTask::WritePsmsToTsv(std::vector<PeptideSpectralMatch*> &psms, const std::string &filePath,
                                          IReadOnlyDictionary<std::string, int> *modstoWritePruned)
    {
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(filePath))
        {
            StreamWriter output = StreamWriter(filePath);
            output.WriteLine(PeptideSpectralMatch::GetTabSeparatedHeader());
            for (auto psm : psms)
            {
                output.WriteLine(psm->ToString(modstoWritePruned));
            }
        }
    }
    
    void MetaMorpheusTask::ReportProgress(ProgressEventArgs *v)
    {
        OutProgressHandler +== nullptr ? nullptr : OutProgressHandler::Invoke(this, v);
    }
    
    void MetaMorpheusTask::FinishedWritingFile(const std::string &path, std::vector<std::string> &nestedIDs)
    {
        SingleFileEventArgs tempVar(path, nestedIDs);
        FinishedWritingFileHandler +== nullptr ? nullptr : FinishedWritingFileHandler::Invoke(this, &tempVar);
    }
    
    void MetaMorpheusTask::StartingDataFile(const std::string &v, std::vector<std::string> &nestedIDs)
    {
        StringEventArgs tempVar(v, nestedIDs);
        StartingDataFileHandler +== nullptr ? nullptr : StartingDataFileHandler::Invoke(this, &tempVar);
    }
    
    void MetaMorpheusTask::FinishedDataFile(const std::string &v, std::vector<std::string> &nestedIDs)
    {
        StringEventArgs tempVar(v, nestedIDs);
        FinishedDataFileHandler +== nullptr ? nullptr : FinishedDataFileHandler::Invoke(this, &tempVar);
    }
    
    void MetaMorpheusTask::Status(const std::string &v, const std::string &id)
    {
        StringEventArgs tempVar(v, new std::vector<std::string> {id});
        OutLabelStatusHandler +== nullptr ? nullptr : OutLabelStatusHandler::Invoke(this, &tempVar);
    }
    
    void MetaMorpheusTask::Status(const std::string &v, std::vector<std::string> &nestedIds)
    {
        StringEventArgs tempVar(v, nestedIds);
        OutLabelStatusHandler +== nullptr ? nullptr : OutLabelStatusHandler::Invoke(this, &tempVar);
    }
    
    void MetaMorpheusTask::Warn(const std::string &v)
    {
        StringEventArgs tempVar(v, nullptr);
        WarnHandler +== nullptr ? nullptr : WarnHandler::Invoke(nullptr, &tempVar);
    }
    
    void MetaMorpheusTask::Log(const std::string &v, std::vector<std::string> &nestedIds)
    {
        StringEventArgs tempVar(v, nestedIds);
        LogHandler +== nullptr ? nullptr : LogHandler::Invoke(this, &tempVar);
    }
    
    void MetaMorpheusTask::NewCollection(const std::string &displayName, std::vector<std::string> &nestedIds)
    {
        StringEventArgs tempVar(displayName, nestedIds);
        NewCollectionHandler +== nullptr ? nullptr : NewCollectionHandler::Invoke(this, &tempVar);
    }
    
    std::vector<std::string> MetaMorpheusTask::GetModsTypesFromString(const std::string &value)
    {
        return value.Split({"\t"}, StringSplitOptions::RemoveEmptyEntries).ToList();
    }
    
    //void MetaMorpheusTask::SingleEngineHandlerInTask(std::any sender, SingleEngineFinishedEventArgs *e)
    void MetaMorpheusTask::SingleEngineHandlerInTask(SingleEngineFinishedEventArgs *e)
    {
        myTaskResults->AddResultText(e->ToString());
    }
    
    void MetaMorpheusTask::FinishedSingleTask(const std::string &displayName)
    {
        SingleTaskEventArgs tempVar(displayName);
        FinishedSingleTaskHandler +== nullptr ? nullptr : FinishedSingleTaskHandler::Invoke(this, &tempVar);
    }
    
    void MetaMorpheusTask::StartingSingleTask(const std::string &displayName)
    {
        SingleTaskEventArgs tempVar(displayName);
        StartingSingleTaskHander +== nullptr ? nullptr : StartingSingleTaskHander::Invoke(this, &tempVar);
    }
    
    std::vector<std::type_info> MetaMorpheusTask::GetSubclassesAndItself(std::type_info type)
    {
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
        yield return type;
    }
    
    bool MetaMorpheusTask::SameSettings(const std::string &pathToOldParamsFile, IndexingEngine *indexEngine)
    {
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamReader reader = new StreamReader(pathToOldParamsFile))
        {
            StreamReader reader = StreamReader(pathToOldParamsFile);
            if (reader.ReadToEnd() == indexEngine->ToString())
            {
                return true;
            }
        }
        return false;
    }
    
    void MetaMorpheusTask::WritePeptideIndex(std::vector<PeptideWithSetModifications*> &peptideIndex, const std::string &peptideIndexFile)
    {
        auto messageTypes = GetSubclassesAndItself(std::vector<PeptideWithSetModifications*>::typeid);
        auto ser = new NetSerializer::Serializer(messageTypes);
        
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.Create(peptideIndexFile))
        {
            auto file = File::Create(peptideIndexFile);
            ser->Serialize(file, peptideIndex);
        }
        
        delete ser;
    }
    
    void MetaMorpheusTask::WriteFragmentIndexNetSerializer(std::vector<std::vector<int>&> &fragmentIndex, const std::string &fragmentIndexFile)
    {
        auto messageTypes = GetSubclassesAndItself(std::vector<std::vector<int>>::typeid);
        auto ser = new NetSerializer::Serializer(messageTypes);
        
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.Create(fragmentIndexFile))
        {
            auto file = File::Create(fragmentIndexFile);
            ser->Serialize(file, fragmentIndex);
        }
        
        delete ser;
    }
    
    std::string MetaMorpheusTask::GetExistingFolderWithIndices(IndexingEngine *indexEngine, std::vector<DbForTask*> &dbFilenameList)
    {
        for (auto database : dbFilenameList)
        {
            std::string baseDir = FileSystem::getDirectoryName(database->getFilePath());
            DirectoryInfo *indexDirectory = new DirectoryInfo(FileSystem::combine(baseDir, IndexFolderName));
            
            if (!FileSystem::directoryExists(indexDirectory->FullName))
            {
                delete indexDirectory;
                return "";
            }
            
            // all directories in the same directory as the protein database
            std::vector<DirectoryInfo*> directories = indexDirectory->GetDirectories();
            
            // look in each subdirectory to find indexes folder
            for (auto possibleFolder : directories)
            {
                std::string result = CheckFiles(indexEngine, possibleFolder);
                
                if (result != "")
                {
                    delete indexDirectory;
                    return result;
                }
            }
            
            delete indexDirectory;
        }
        
        return "";
    }
    
    std::string MetaMorpheusTask::CheckFiles(IndexingEngine *indexEngine, DirectoryInfo *folder)
    {
        if (FileSystem::fileExists(FileSystem::combine(folder->FullName, "indexEngine.params"))    &&
            FileSystem::fileExists(FileSystem::combine(folder->FullName, "peptideIndex.ind"))      &&
            FileSystem::fileExists(FileSystem::combine(folder->FullName, "fragmentIndex.ind"))     &&
            (FileSystem::fileExists(FileSystem::combine(folder->FullName, "precursorIndex.ind")) ||
             !indexEngine->GeneratePrecursorIndex)                                                 &&
            SameSettings(FileSystem::combine(folder->FullName, "indexEngine.params"), indexEngine))
        {
            return folder->FullName;
        }
        return "";
    }
    
    void MetaMorpheusTask::WriteIndexEngineParams(IndexingEngine *indexEngine, const std::string &fileName)
    {
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(fileName))
        {
            StreamWriter output = StreamWriter(fileName);
            output.Write(indexEngine);
        }
    }
    
    std::string MetaMorpheusTask::GenerateOutputFolderForIndices(std::vector<DbForTask*> &dbFilenameList)
    {
        auto pathToIndexes = FileSystem::combine(FileSystem::getDirectoryName(dbFilenameList.front().FilePath), IndexFolderName);
        if (!FileSystem::fileExists(pathToIndexes))
        {
            FileSystem::createDirectory(pathToIndexes);
        }
        auto folder = FileSystem::combine(pathToIndexes, DateTime::Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo::InvariantCulture));
        FileSystem::createDirectory(folder);
        return folder;
    }
    
    void MetaMorpheusTask::GenerateIndexes(IndexingEngine *indexEngine, std::vector<DbForTask*> &dbFilenameList,
                                           std::vector<PeptideWithSetModifications*> &peptideIndex,
                                           std::vector<std::vector<int>> &fragmentIndex,
                                           std::vector<std::vector<int>> &precursorIndex,
                                           std::vector<Protein*> &allKnownProteins,
                                           std::vector<Modification*> &allKnownModifications,
                                           const std::string &taskId)
    {
        std::string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);
        if (pathToFolderWithIndices == "")
        {
            auto output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
            Status("Writing params...", std::vector<std::string> {taskId});
            auto paramsFile = FileSystem::combine(output_folderForIndices, "indexEngine.params");
            WriteIndexEngineParams(indexEngine, paramsFile);
            FinishedWritingFile(paramsFile, std::vector<std::string> {taskId});
            
            Status("Running Index Engine...", std::vector<std::string> {taskId});
            auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
            peptideIndex = indexResults->getPeptideIndex();
            fragmentIndex = indexResults->getFragmentIndex();
            precursorIndex = indexResults->getPrecursorIndex();
            
            Status("Writing peptide index...", std::vector<std::string> {taskId});
            auto peptideIndexFile = FileSystem::combine(output_folderForIndices, "peptideIndex.ind");
            WritePeptideIndex(peptideIndex, peptideIndexFile);
            FinishedWritingFile(peptideIndexFile, std::vector<std::string> {taskId});
            
            Status("Writing fragment index...", std::vector<std::string> {taskId});
            auto fragmentIndexFile = FileSystem::combine(output_folderForIndices, "fragmentIndex.ind");
            WriteFragmentIndexNetSerializer(fragmentIndex, fragmentIndexFile);
            FinishedWritingFile(fragmentIndexFile, std::vector<std::string> {taskId});
            
            if (indexEngine->GeneratePrecursorIndex)
            {
                Status("Writing precursor index...", std::vector<std::string> {taskId});
                auto precursorIndexFile = FileSystem::combine(output_folderForIndices, "precursorIndex.ind");
                WriteFragmentIndexNetSerializer(precursorIndex, precursorIndexFile);
                FinishedWritingFile(precursorIndexFile, std::vector<std::string> {taskId});
            }
        }
        else
        {
            Status("Reading peptide index...", std::vector<std::string> {taskId});
            auto messageTypes = GetSubclassesAndItself(std::vector<PeptideWithSetModifications*>::typeid);
            auto ser = new NetSerializer::Serializer(messageTypes);
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "peptideIndex.ind")))
            {
                auto file = File::OpenRead(FileSystem::combine(pathToFolderWithIndices, "peptideIndex.ind"));
                peptideIndex = static_cast<std::vector<PeptideWithSetModifications*>>(ser->Deserialize(file));
            }
            
            // populate dictionaries of known proteins for deserialization
            std::unordered_map<std::string, Protein*> proteinDictionary;
            
            for (auto protein : allKnownProteins)
            {
                if (proteinDictionary.find(protein->Accession) == proteinDictionary.end())
                {
                    proteinDictionary.emplace(protein->Accession, protein);
                }
                else if (proteinDictionary[protein->Accession]->BaseSequence != protein->BaseSequence)
                {
                    delete ser;
                    throw MetaMorpheusException(StringHelper::formatSimple("The protein database contained multiple proteins with accession {0} ! This is not allowed for index-based searches (modern, non-specific, crosslink searches)", protein->Accession));
                }
            }
            
            // get non-serialized information for the peptides (proteins, mod info)
            for (auto peptide : peptideIndex)
            {
                peptide->SetNonSerializedPeptideInfo(GlobalVariables::getAllModsKnownDictionary(), proteinDictionary);
            }
            
            Status("Reading fragment index...", std::vector<std::string> {taskId});
            messageTypes = GetSubclassesAndItself(std::vector<std::vector<int>>::typeid);
            ser = new NetSerializer::Serializer(messageTypes);
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "fragmentIndex.ind")))
            {
                auto file = File::OpenRead(FileSystem::combine(pathToFolderWithIndices, "fragmentIndex.ind"));
                fragmentIndex = static_cast<std::vector<std::vector<int>>>(ser->Deserialize(file));
            }
            
            if (indexEngine->GeneratePrecursorIndex)
            {
                Status("Reading precursor index...", std::vector<std::string> {taskId});
                messageTypes = GetSubclassesAndItself(std::vector<std::vector<int>>::typeid);
                ser = new NetSerializer::Serializer(messageTypes);
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "precursorIndex.ind")))
                {
                    auto file = File::OpenRead(FileSystem::combine(pathToFolderWithIndices, "precursorIndex.ind"));
                    precursorIndex = static_cast<std::vector<std::vector<int>>>(ser->Deserialize(file));
                }
            }
            
            delete ser;
        }
    }
}
