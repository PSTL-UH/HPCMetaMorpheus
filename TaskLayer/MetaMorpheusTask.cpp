﻿#include "MetaMorpheusTask.h"
#include "../EngineLayer/CommonParameters.h"
#include "MyTaskResults.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "FileSpecificParameters.h"
#include "DbForTask.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/EventArgs/SingleEngineFinishedEventArgs.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/MetaMorpheusException.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../EngineLayer/EventArgs/SingleFileEventArgs.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "EventArgs/SingleTaskEventArgs.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Nett;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{

    //Dr. Gabriel, I havent sorted out what to do with this yet.  It looks to me like they are creating and setting the tomlConfig key-value pairs
    //from the Tolerance, PpmTolerance, AbsoluteTolerance, and Protease configurations.  Its a little confusing because not all of these 
    //are member variables present in MetaMorpheusTask.h or FileSpecificParameters.h.  It looks like for each of these, they are converting
    //the value to a string and then performing some operation with that string.
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
        auto ms2Scans = myMSDataFile->GetAllScansList().Where([&] (std::any x)		{
                return x::MsnOrder > 1;
            })->ToArray();
        std::vector<std::vector<Ms2ScanWithSpecificMass*>> scansWithPrecursors(ms2Scans.size());

        ParallelOptions *tempVar = new ParallelOptions();
        tempVar->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
        Parallel::ForEach(Partitioner::Create(0, ms2Scans.size()), tempVar, [&] (partitionRange, loopState)       {
                for (int i = partitionRange::Item1; i < partitionRange::Item2; i++)
                {
                    if (GlobalVariables::getStopLoops())
                    {
                        break;
                    }
                    
                    MsDataScan *ms2scan = ms2Scans[i];
                    
                    std::vector<(double, int)*> precursors;
                    if (ms2scan->OneBasedPrecursorScanNumber.HasValue)
                    {
                        auto precursorSpectrum = myMSDataFile->GetOneBasedScan(ms2scan->OneBasedPrecursorScanNumber->Value);
                        
                        try
                        {
                            ms2scan->RefineSelectedMzAndIntensity(precursorSpectrum->MassSpectrum);
                        }
                        catch (const MzLibException &ex)
                        {
                            Warn("Could not get precursor ion for MS2 scan #" + ms2scan->OneBasedScanNumber + "; " + ex->Message);
                            continue;
                        }
                        
                        if (ms2scan->SelectedIonMonoisotopicGuessMz.HasValue)
                        {
                            ms2scan->ComputeMonoisotopicPeakIntensity(precursorSpectrum->MassSpectrum);
                        }
                        
                        if (commonParameters->getDoPrecursorDeconvolution())
                        {
                            for (auto envelope : ms2scan->GetIsolatedMassesAndCharges(precursorSpectrum->MassSpectrum, 1,
                                                                                      commonParameters->getDeconvolutionMaxAssumedChargeState(),
                                                                                      commonParameters->getDeconvolutionMassTolerance()->Value,
                                                                                      commonParameters->getDeconvolutionIntensityRatio()))
                            {
                                auto monoPeakMz = envelope->monoisotopicMass.ToMz(envelope->charge);
                                precursors.Add((monoPeakMz, envelope->charge));
                            }
                        }
                    }
                    
                    if (commonParameters->getUseProvidedPrecursorInfo() && ms2scan->SelectedIonChargeStateGuess.HasValue)
                    {
                        auto precursorCharge = ms2scan->SelectedIonChargeStateGuess->Value;
                        if (ms2scan->SelectedIonMonoisotopicGuessMz.HasValue)
                        {
                            double precursorMZ = ms2scan->SelectedIonMonoisotopicGuessMz->Value;
                            if (!precursors.Any([&] (std::any b)   {
                                        commonParameters->getDeconvolutionMassTolerance()->Within(precursorMZ.ToMass(precursorCharge), b::Item1->ToMass(b::Item2));
                                    }))
                            {
                                precursors.Add((precursorMZ, precursorCharge));
                            }
                        }
                        else
                        {
                            double precursorMZ = ms2scan->SelectedIonMZ->Value;
                            if (!precursors.Any([&] (std::any b){
                                        commonParameters->getDeconvolutionMassTolerance()->Within(precursorMZ.ToMass(precursorCharge), b::Item1->ToMass(b::Item2));
                                    }))
                            {
                                precursors.Add((precursorMZ, precursorCharge));
                            }
                        }
                    }
                    
                    scansWithPrecursors[i] = std::vector<Ms2ScanWithSpecificMass*>();
                    std::vector<IsotopicEnvelope*> neutralExperimentalFragments = Ms2ScanWithSpecificMass::GetNeutralExperimentalFragments(ms2scan, commonParameters);
                    
                    for (auto precursor : precursors)
                    {
                        Ms2ScanWithSpecificMass tempVar2(ms2scan, precursor->Item1, precursor->Item2, fullFilePath, commonParameters, neutralExperimentalFragments);
                        scansWithPrecursors[i].Add(&tempVar2);
                    }
                }
            });
        
        return scansWithPrecursors.SelectMany([&] (std::any p)	{
                return p;
            });
    }
    
    EngineLayer::CommonParameters *MetaMorpheusTask::SetAllFileSpecificCommonParams(EngineLayer::CommonParameters *commonParams, FileSpecificParameters *fileSpecificParams)
    {
        if (fileSpecificParams == nullptr)
        {
            return commonParams;
        }
        
        // set file-specific digestion parameters
        Protease tempVar = fileSpecificParams.getProtease();
        Protease *protease = (tempVar != nullptr) ? tempVar : commonParams->getDigestionParams()->Protease;
        std::optional<int> tempVar2 = fileSpecificParams.getMinPeptideLength();
        int minPeptideLength = tempVar2 ? tempVar2 : commonParams->getDigestionParams()->MinPeptideLength;
        std::optional<int> tempVar3 = fileSpecificParams.getMaxPeptideLength();
        int maxPeptideLength = tempVar3 ? tempVar3 : commonParams->getDigestionParams()->MaxPeptideLength;
        std::optional<int> tempVar4 = fileSpecificParams.getMaxMissedCleavages();
        int maxMissedCleavages = tempVar4 ? tempVar4 : commonParams->getDigestionParams()->MaxMissedCleavages;
        std::optional<int> tempVar5 = fileSpecificParams.getMaxModsForPeptide();
        int maxModsForPeptide = tempVar5 ? tempVar5 : commonParams->getDigestionParams()->MaxModsForPeptide;

        DigestionParams *fileSpecificDigestionParams = new DigestionParams(protease: protease->Name,
                                                                           maxMissedCleavages: maxMissedCleavages,
                                                                           minPeptideLength: minPeptideLength,
                                                                           maxPeptideLength: maxPeptideLength,
                                                                           maxModsForPeptides: maxModsForPeptide,
                                                                           maxModificationIsoforms: commonParams->getDigestionParams()->MaxModificationIsoforms,
                                                                           initiatorMethionineBehavior: commonParams->getDigestionParams()->InitiatorMethionineBehavior,
                                                                           fragmentationTerminus: commonParams->getDigestionParams()->FragmentationTerminus,
                                                                           searchModeType: commonParams->getDigestionParams()->SearchModeType);
        
        // set the rest of the file-specific parameters
        Tolerance tempVar6 = fileSpecificParams.getPrecursorMassTolerance();
        Tolerance *precursorMassTolerance = (tempVar6 != nullptr) ? tempVar6 : commonParams->getPrecursorMassTolerance();
        Tolerance tempVar7 = fileSpecificParams.getProductMassTolerance();
        Tolerance *productMassTolerance = (tempVar7 != nullptr) ? tempVar7 : commonParams->getProductMassTolerance();
        std::optional<DissociationType*> tempVar8 = fileSpecificParams.getDissociationType();
        DissociationType *dissociationType = tempVar8 ? tempVar8 : commonParams->getDissociationType();
        
        
        EngineLayer::CommonParameters *returnParams = new CommonParameters(commonParams->getTaskDescriptor(),
                                                                           dissociationType,
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
                                                                           commonParams->ListOfModsVariable,
                                                                           commonParams->ListOfModsFixed,
                                                                           commonParams->getQValueOutputFilter(),
                                                                           commonParams->getAssumeOrphanPeaksAreZ1Fragments());

        //C# TO C++ CONVERTER TODO TASK: A 'delete returnParams' statement was not added since returnParams was used in a 'return' or 'throw' statement.
        delete fileSpecificDigestionParams;
        return returnParams;
    }

    MyTaskResults *MetaMorpheusTask::RunTask(const std::string &output_folder, std::vector<DbForTask*> &currentProteinDbFilenameList,
                                             std::vector<std::string> &currentRawDataFilepathList, const std::string &displayName)
    {
        StartingSingleTask(displayName);
        
#ifdef ORIG
        auto tomlFileName = FileSystem::combine(Directory::GetParent(output_folder)->ToString(), "Task Settings", displayName + "config.toml");
        Toml::WriteFile(this, tomlFileName, tomlConfig);
        FinishedWritingFile(tomlFileName, std::vector<std::string> {displayName});
#endif
        //looks like the file is being written to the parent directory of the output_folder
        //we need the output folder as type std::experimental::filesystem::path rather than std::string
        //I made a new function for writing new toml files that does not require setting key-value pairs, it just outputs the Toml value to the new file.
        //However I havent sorted out how to create the tomlConfig yet.
        std::experimental::filesystem::path output_directory = output_folder;
        std::string output_path = output_directory.parent_path().string() + "Task Settings" + displayName + "config.toml";
        Toml trw;
        trw.tomlWriteNewFile(tomlFileName, tomlConfig);

        
        MetaMorpheusEngine::FinishedSingleEngineHandler->addListener("SingleEngineHandlerInTask", [&] (std::any sender, SingleEngineFinishedEventArgs* e) {
                SingleEngineHandlerInTask(sender, e);});
        try
        {
            auto stopWatch = new Stopwatch();
            stopWatch->Start();
            
            std::vector<FileSpecificParameters*> fileSettingsList(currentRawDataFilepathList.size());
            for (int i = 0; i < currentRawDataFilepathList.size(); i++)
            {
                if (GlobalVariables::getStopLoops())
                {
                    break;
                }
                std::string rawFilePath = currentRawDataFilepathList[i];
                std::string directory = Directory::GetParent(rawFilePath)->ToString();
                std::string fileSpecificTomlPath = FileSystem::combine(directory, Path::GetFileNameWithoutExtension(rawFilePath)) + ".toml";
                if (FileSystem::fileExists(fileSpecificTomlPath))
                {
#ifdef ORIG
                    //Dr. Gabriel, Im not quite sure why they are passing the tomlConfig to the ReadFile function here
                    //when it looks like theyre setting fileSpecificSettings.  We do need to find a table header to make a toml table
                    //maybe this is the reason.  But I'm not sure what the header would be
                    TomlTable *fileSpecificSettings = Toml::ReadFile(fileSpecificTomlPath, tomlConfig);
#endif
                    toml::Value toml_value = trw.tomlReadFile(fileSpecificTomlPath);

                    //here we need a header to search for to construct the toml table.  I assume its in tomlConfig, but Im not sure what it is,
                    //so as it is this now wont work yet.
                    toml::Value* fileParameters = trw.getValue(toml_value, tomlConfig);
                    toml::Table fileSpecificSettings = fileParameters->as<toml::Table>();


                    try
                    {
                        fileSettingsList[i] = new FileSpecificParameters(fileSpecificSettings);
                    }
                    catch (const MetaMorpheusException &e)
                    {
                        // file-specific toml has already been validated in the GUI when the spectra files were added, so...
                        // probably the only time you can get here is if the user modifies the file-specific parameter file in the middle of a run...
                        Warn("Problem parsing the file-specific toml " + FileSystem::getFileName(fileSpecificTomlPath) + "; " + e->what() +
                             "; is the toml from an older version of MetaMorpheus?");
                    }
                }
            }
            
            RunSpecific(output_folder, currentProteinDbFilenameList, currentRawDataFilepathList, displayName, fileSettingsList);
            stopWatch->Stop();
            myTaskResults->Time = stopWatch->Elapsed;
            auto resultsFileName = FileSystem::combine(output_folder, "results.txt");
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                StreamWriter file = StreamWriter(resultsFileName);
                file.WriteLine("MetaMorpheus: version " + GlobalVariables::getMetaMorpheusVersion());
                file.Write(myTaskResults->ToString());
            }
            FinishedWritingFile(resultsFileName, std::vector<std::string> {displayName});
            FinishedSingleTask(displayName);
            
            delete stopWatch;
        }
        catch (const std::runtime_error &e)
        {
            MetaMorpheusEngine::FinishedSingleEngineHandler->removeListener("SingleEngineHandlerInTask");
            auto resultsFileName = FileSystem::combine(output_folder, "results.txt");
            e.Data->Add("folder", output_folder);
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                StreamWriter file = StreamWriter(resultsFileName);
                file.WriteLine(GlobalVariables::getMetaMorpheusVersion() == "1.0.0.0" ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + GlobalVariables::getMetaMorpheusVersion());
                file.WriteLine(SystemInfo::CompleteSystemInfo()); //OS, OS Version, .Net Version, RAM, processor count, MSFileReader .dll versions X3
                file.Write("e: " + e);
                file.Write("e.Message: " + e.what());
                file.Write("e.InnerException: " + e.InnerException);
                file.Write("e.Source: " + e.Source);
                file.Write("e.StackTrace: " + e.StackTrace);
                file.Write("e.TargetSite: " + e.TargetSite);
            }
            throw;
        }
        
        {
            auto proseFilePath = FileSystem::combine(output_folder, "prose.txt");
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter file = new StreamWriter(proseFilePath))
            {
                StreamWriter file = StreamWriter(proseFilePath);
                file.Write("The data analysis was performed using MetaMorpheus version " + GlobalVariables::getMetaMorpheusVersion() + ", available at " +
                           "https://github.com/smith-chem-wisc/MetaMorpheus.");
                file.Write(ProseCreatedWhileRunning->toString());
                file.Write(SystemInfo::SystemProse()->Replace("\r\n", "") + " ");
                file.WriteLine("The total time to perform the " + getTaskType() + " task on " + std::to_string(currentRawDataFilepathList.size()) +
                               " spectra file(s) was " + std::string::Format("{0:0.00}", myTaskResults->Time.TotalMinutes) + " minutes.");
                file.WriteLine();
                file.WriteLine("Published works using MetaMorpheus software are encouraged to cite: Solntsev, S. K.; Shortreed, M. R.; Frey, B. L.; Smith, L. M. Enhanced Global Post-translational Modification Discovery with MetaMorpheus. Journal of Proteome Research. 2018, 17 (5), 1844-1851.");

                file.WriteLine();
                file.WriteLine("Spectra files: ");
                file.WriteLine(std::string::Join("\r\n", currentRawDataFilepathList.Select([&] (std::any b)    {
                                return '\t' + b;
                            })));
                file.WriteLine("Databases:");
                file.Write(std::string::Join("\r\n", currentProteinDbFilenameList.Select([&] (std::any b)		{
                                return '\t' + (b::IsContaminant ? "Contaminant " : "") + b::FilePath;
                            })));
            }
            FinishedWritingFile(proseFilePath, std::vector<std::string> {displayName});
        }
        
        MetaMorpheusEngine::FinishedSingleEngineHandler->removeListener("SingleEngineHandlerInTask");
        return myTaskResults;
    }
    
    std::vector<Protein*> MetaMorpheusTask::LoadProteins(const std::string &taskId, std::vector<DbForTask*> &dbFilenameList, bool searchTarget,
                                                         DecoyType decoyType, std::vector<std::string> &localizeableModificationTypes,
                                                         EngineLayer::CommonParameters *commonParameters)
    {
        Status("Loading proteins...", std::vector<std::string> {taskId});
        int emptyProteinEntries = 0;
        std::vector<Protein*> proteinList;
        for (auto db : dbFilenameList)
        {
            int emptyProteinEntriesForThisDb = 0;
            Dictionary<std::string, Modification*> unknownModifications;
            auto dbProteinList = LoadProteinDb(db->getFilePath(), searchTarget, decoyType, localizeableModificationTypes, db->getIsContaminant(),
                                               unknownModifications, emptyProteinEntriesForThisDb, commonParameters);
            proteinList = proteinList.Concat(dbProteinList)->ToList();
            emptyProteinEntries += emptyProteinEntriesForThisDb;
        }
        if (!proteinList.Any())
        {
            Warn("Warning: No protein entries were found in the database");
        }
        else if (emptyProteinEntries > 0)
        {
            Warn("Warning: " + std::to_string(emptyProteinEntries) + " empty protein entries ignored");
        }
        return proteinList;
    }
    
    std::vector<Protein*> MetaMorpheusTask::LoadProteinDb(const std::string &fileName, bool generateTargets, DecoyType decoyType,
                                                          std::vector<std::string> &localizeableModificationTypes, bool isContaminant,
                                                          std::unordered_map<std::string, Modification*> &um, int &emptyEntriesCount,
                                                          EngineLayer::CommonParameters *commonParameters)
    {
        std::vector<std::string> dbErrors;
        std::vector<Protein*> proteinList;
        
//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
        std::string theExtension = Path::GetExtension(fileName).ToLowerInvariant();
        bool compressed = StringHelper::endsWith(theExtension, "gz"); // allows for .bgz and .tgz, too which are used on occasion
//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
        theExtension = compressed ? Path::GetExtension(Path::GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;
        
        if (theExtension == ".fasta" || theExtension == ".fa")
        {
            um.clear();
            proteinList = ProteinDbLoader::LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant,
                                                            ProteinDbLoader::UniprotAccessionRegex, ProteinDbLoader::UniprotFullNameRegex,
                                                            ProteinDbLoader::UniprotFullNameRegex, ProteinDbLoader::UniprotGeneNameRegex,
                                                            ProteinDbLoader::UniprotOrganismRegex, dbErrors,
                                                            commonParameters->getMaxThreadsToUsePerFile());
        }
        else
        {
            std::vector<std::string> modTypesToExclude = GlobalVariables::getAllModTypesKnown().Where([&] (std::any b)  {
                    std::find(localizeableModificationTypes.begin(), localizeableModificationTypes.end(), b) == localizeableModificationTypes.end();
                }).ToList();
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
                                             std::vector<Modification*> &fixedModifications, std::vector<std::string> &localizableModificationTypes)
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
    
    void MetaMorpheusTask::SingleEngineHandlerInTask(std::any sender, SingleEngineFinishedEventArgs *e)
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
