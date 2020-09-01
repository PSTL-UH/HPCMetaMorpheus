#include "SearchTask.h"
#include "SearchParameters.h"
#include "../../EngineLayer/CommonParameters.h"
#include "../../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/SingleAbsoluteAroundZeroMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/DotMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/IntervalMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/OpenMassDiffAcceptor.h"
#include "../../EngineLayer/MetaMorpheusException.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"
#include "../MyTaskResults.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../../EngineLayer/ProteinScoringAndFdr/FdrCategory.h"
#include "../MyFileManager.h"
#include "../../EngineLayer/GlobalVariables.h"
#include "../../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "PostSearchAnalysisParameters.h"
#include "PostSearchAnalysisTask.h"
#include "UsefulProteomicsDatabases/DecoyType.h"

#include "stringhelper.h"


namespace TaskLayer
{

    SearchTask::SearchTask() : MetaMorpheusTask(MyTask::Search)
    {
        auto tempVar = new EngineLayer::CommonParameters();
        setCommonParameters(tempVar);
        
        auto tempVar2 = new TaskLayer::SearchParameters();
        setSearchParameters(tempVar2);
    }
    
    SearchTask::SearchTask(std::string tomlFile ) : MetaMorpheusTask(MyTask::Search)
    {
        //Check ending that it is a toml file.

        Toml trw;
        toml::Value toml_value = trw.tomlReadFile(tomlFile);
        toml::Value* fileParameters = trw.getValue(toml_value, "SearchParameters");
        toml::Table tomlTable = fileParameters->as<toml::Table>();

        auto sParams = new TaskLayer::SearchParameters();

        // parse toml file and set the values
        for (auto const& keyValuePair : tomlTable)
        {
            if ( keyValuePair.first == "DisposeOfFileWhenDone" ) {
                auto val = keyValuePair.second.as<bool>();
                sParams->setDisposeOfFileWhenDone (val);
            }
            else if ( keyValuePair.first == "DoParsimony" ) {
                auto val = keyValuePair.second.as<bool>();
                sParams->setDoParsimony (val);
            }
            else if ( keyValuePair.first == "ModPeptidesAreDifferent" ) {
                auto val = keyValuePair.second.as<bool>();
                sParams->setModPeptidesAreDifferent (val);
            }
            else if ( keyValuePair.first == "NoOneHitWonders" ) {
                auto val = keyValuePair.second.as<bool>();
                sParams->setNoOneHitWonders(val);
            }
            else if ( keyValuePair.first == "MatchBetweenRuns" ) {
                auto val = keyValuePair.second.as<bool>();
                sParams->setMatchBetweenRuns(val);
            }
            else if ( keyValuePair.first == "Normalize" ) {
                auto val = keyValuePair.second.as<bool>();
                sParams->setNormalize(val);
            }
            else if ( keyValuePair.first == "QuantifyPpmTol" ) {
                auto val = keyValuePair.second.as<double>();
                sParams->setQuantifyPpmTol(val);
            }
            else if ( keyValuePair.first == "DoHistogramAnalysis" ) {
                auto val = keyValuePair.second.as<bool>();
                sParams->setDoHistogramAnalysis(val);
            }
            else if ( keyValuePair.first == "SearchTarget" ) {
                auto val = keyValuePair.second.as<bool>();
                sParams->setSearchTarget(val);
            }
            else if (keyValuePair.first == "DecoyType") {
                auto val = keyValuePair.second.as<std::string>();
                sParams->setDecoyType(DecoyTypeFromString(val ));
            }
            else if (keyValuePair.first == "MassDiffAcceptorType") {
                auto val = keyValuePair.second.as<std::string>();
                sParams->setMassDiffAcceptorType(TaskLayer::MassDiffAcceptorTypeFromString(val));
            }
            else if (keyValuePair.first == "WritePrunedDatabase") {
                auto val = keyValuePair.second.as<bool>();
                sParams->setWritePrunedDatabase(val);
            }
            else if (keyValuePair.first == "KeepAllUniprotMods") {
                auto val = keyValuePair.second.as<bool>();
                sParams->setKeepAllUniprotMods(val);
            }
            else if (keyValuePair.first == "DoLocalizationAnalysis") {
                auto val = keyValuePair.second.as<bool>();
                sParams->setDoLocalizationAnalysis(val);
            }
            else if (keyValuePair.first == "DoQuantification") {
                auto val = keyValuePair.second.as<bool>();
                sParams->setDoQuantification(val);
            }
            else if (keyValuePair.first == "SearchType") {
                auto val = keyValuePair.second.as<std::string>();
                sParams->setSearchType(TaskLayer::SearchTypeFromString(val));
            }
            else if (keyValuePair.first == "LocalFdrCategories") {
                std::vector<FdrCategory> vec;
                auto val = keyValuePair.second.as<toml::Array>();
                for ( auto v : val ) {
                    vec.push_back (EngineLayer::FdrCategoryFromString(v.as<std::string>()));
                }
                sParams->setLocalFdrCategories(vec);
            }
            else if (keyValuePair.first == "MaxFragmentSize") {
                auto val = keyValuePair.second.as<double>();
                sParams->setMaxFragmentSize(val);
            }
            else if (keyValuePair.first == "HistogramBinTolInDaltons") {
                auto val = keyValuePair.second.as<double>();
                sParams->setHistogramBinTolInDaltons(val);
            }
            else if (keyValuePair.first == "MaximumMassThatFragmentIonScoreIsDoubled") {
                auto val = keyValuePair.second.as<double>();
                sParams->setMaximumMassThatFragmentIonScoreIsDoubled(val);
            }
            else if (keyValuePair.first == "WriteMzId") {
                auto val = keyValuePair.second.as<bool>();
                sParams->setWriteMzId(val);
            }
            else if (keyValuePair.first == "WritePepXml") {
                auto val = keyValuePair.second.as<bool>();
                sParams->setWritePepXml(val);
            }
            else if (keyValuePair.first == "WriteDecoys") {
                auto val = keyValuePair.second.as<bool>();
                sParams->setWriteDecoys(val);
            }
            else if (keyValuePair.first == "WriteContaminants") {
                auto val = keyValuePair.second.as<bool>();
                sParams->setWriteContaminants(val);
            }
        }
        
        toml::Value* fileParameters2 = trw.getValue(toml_value, "SearchParameters.ModsToWriteSelection");
        toml::Table tomlTable2 = fileParameters2->as<toml::Table>();
        std::unordered_map<std::string, int> modsmap;
        for (auto b = tomlTable2.begin(); b != tomlTable2.end(); b++ ) {
            std::string key = (*b).first;
            int value = (*b).second.as<int>();
            modsmap.insert ({key,value});
        }
        sParams->setModsToWriteSelection(modsmap);
        setSearchParameters (sParams);
        
        std::string taskDescr = "";
        auto tempVar = new EngineLayer::CommonParameters( tomlFile, taskDescr );
        setCommonParameters(tempVar);
    }
    
    void  SearchTask::writeTomlConfig(std::string &filename, std::ofstream &tomlFd )
    {
        if ( !tomlFd.is_open() ) {
            tomlFd.open(filename );
            if ( !tomlFd.is_open() ) {
                std::cout << "Could not open file " << filename << std::endl;
                return;
            }
        }
        
        toml::Value v;
        std::string key = "TaskType", value = "Search";
        v.set ( key, value);
        tomlFd << v;

        SearchParameters *sparams = getSearchParameters();
        toml::Table search_params;

        search_params["DisposeOfFileWhenDone"] = sparams->getDisposeOfFileWhenDone();
        search_params["DoParsimony"] = sparams->getDoParsimony();
        search_params["ModPeptidesAreDifferent"] = sparams->getModPeptidesAreDifferent();
        search_params["NoOneHitWonders"] = sparams->getNoOneHitWonders();
        search_params["MatchBetweenRuns"] = sparams->getMatchBetweenRuns();
        search_params["Normalize"] = sparams->getNormalize();
        search_params["QuantifyPpmTol"] = sparams->getQuantifyPpmTol();
        search_params["DoHistogramAnalysis"] = sparams->getDoHistogramAnalysis();
        search_params["SearchTarget"] = sparams->getSearchTarget();
        auto var = sparams->getDecoyType();
        search_params["DecoyType"] =  UsefulProteomicsDatabases::DecoyTypeToString(var);
        search_params["MassDiffAcceptorType"] = TaskLayer::MassDiffAcceptorTypeToString(sparams->getMassDiffAcceptorType());
        search_params["WritePrunedDatabase"] = sparams->getWritePrunedDatabase();
        search_params["KeepAllUniprotMods"] = sparams->getKeepAllUniprotMods();
        search_params["DoLocalizationAnalysis"] = sparams->getDoLocalizationAnalysis();
        search_params["DoQuantification"] = sparams->getDoQuantification();
        search_params["SearchType"] = TaskLayer::SearchTypeToString(sparams->getSearchType());

        toml::Array sarr;
        auto var2 = sparams->getLocalFdrCategories();
        for ( auto t: var2 ) {
            sarr.push_back(EngineLayer::FdrCategoryToString(t));
        }
        search_params["LocalFdrCategories"] = sarr;
        search_params["MaxFragmentSize"] = sparams->getMaxFragmentSize();
        search_params["HistogramBinTolInDaltons"] = sparams->getHistogramBinTolInDaltons();
        search_params["MaximumMassThatFragmentIonScoreIsDoubled"] = sparams->getMaximumMassThatFragmentIonScoreIsDoubled();
        search_params["WriteMzId"] = sparams->getWriteMzId();
        search_params["WritePepXml"] = sparams->getWritePepXml();
        search_params["WriteDecoys"] = sparams->getWriteDecoys();
        search_params["WriteContaminants"] = sparams->getWriteContaminants();
        
        tomlFd << std::endl;
        tomlFd << "[SearchParameters]" << std::endl;
        tomlFd << search_params;        
        
        std::unordered_map<std::string, int> modsmap = sparams->getModsToWriteSelection();
        toml::Value v2;
        for ( auto p: modsmap ) {
            std::string key = "\'" + p.first + "\'";
            int value = p.second;
            v2.set ( key, value);
        }
        tomlFd << std::endl;
        tomlFd << "[SearchParameters.ModsToWriteSelection]" << std::endl;
        tomlFd << v2;
        
        MetaMorpheusTask::writeTomlConfig(filename, tomlFd );
        if ( tomlFd.is_open() ) {
            tomlFd.close();
        }
        return;
    }
    TaskLayer::SearchParameters *SearchTask::getSearchParameters() const
    {
        return privateSearchParameters;
    }
    
    void SearchTask::setSearchParameters(TaskLayer::SearchParameters *value)
    {
        if ( privateSearchParameters != nullptr ) {
            delete privateSearchParameters;
        }
        privateSearchParameters = value;
    }
    
    MassDiffAcceptor *SearchTask::GetMassDiffAcceptor(Tolerance *precursorMassTolerance,
                                                      MassDiffAcceptorType massDiffAcceptorType,
                                                      const std::string &customMdac)
    {
        if ( massDiffAcceptorType ==  MassDiffAcceptorType::Exact  )
        {
            if (dynamic_cast<PpmTolerance*>(precursorMassTolerance) != nullptr)
            {
                return new SinglePpmAroundZeroSearchMode(precursorMassTolerance->getValue());
            }
            else
            {
                return new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance->getValue());
            }
        }
        else if ( massDiffAcceptorType == MassDiffAcceptorType::OneMM )
        {
            std::vector<double> vec = {0, 1.0029};
            return new DotMassDiffAcceptor("1mm", vec, precursorMassTolerance);
        }
        
        else if ( massDiffAcceptorType ==MassDiffAcceptorType::TwoMM )
        {
            std::vector<double> vec1 = {0, 1.0029, 2.0052};
            return new DotMassDiffAcceptor("2mm", vec1, precursorMassTolerance);
        }    
        else if ( massDiffAcceptorType == MassDiffAcceptorType::ThreeMM ) 
        {
            std::vector<double> vec2 = {0, 1.0029, 2.0052, 3.0077};
            return new DotMassDiffAcceptor("3mm", vec2, precursorMassTolerance);
        }
        else if( massDiffAcceptorType == MassDiffAcceptorType::ModOpen )
        {
            std::vector<DoubleRange*> v1 = {new DoubleRange(-187, std::numeric_limits<double>::infinity())};
            return new IntervalMassDiffAcceptor("-187andUp", v1);
        }    
        else if (massDiffAcceptorType == MassDiffAcceptorType::Open ) 
        {
            return new OpenSearchMode();
        }       
        else if ( massDiffAcceptorType == MassDiffAcceptorType::Custom ) 
        {
            return ParseSearchMode(customMdac);
        }

        //default:
        throw MetaMorpheusException("Unknown MassDiffAcceptorType");
        return nullptr;
    }
    
    MyTaskResults *SearchTask::RunSpecific(const std::string &OutputFolder,
                                           std::vector<DbForTask*> &dbFilenameList,
                                           std::vector<std::string> &currentRawFileList,
                                           const std::string &taskId,
                                           std::vector<FileSpecificParameters*> &fileSettingsList)
    {
        // disable quantification if a .mgf is being used
        bool found = false;
        for ( auto x: currentRawFileList ) {
            int pos = x.find_last_of(".")+1;
            std::string extension = x.substr(pos, x.length() );
            if ( extension == "mgf" ) {
                found = true;
                break;
            }
        }                    
        
        if ( found )     
        {
            getSearchParameters()->setDoQuantification(false);
        }
        
        std::vector<Modification*> variableModifications;
        std::vector<Modification*> fixedModifications;
        std::vector<std::string> localizeableModificationTypes;
        LoadModifications(taskId, variableModifications, fixedModifications, localizeableModificationTypes);
        
        // load proteins
        std::vector<Protein*> proteinList = LoadProteins(taskId, dbFilenameList,
                                                         getSearchParameters()->getSearchTarget(),
                                                         getSearchParameters()->getDecoyType(),
                                                         localizeableModificationTypes,
                                                         getCommonParameters());
        
        // write prose settings
        ProseCreatedWhileRunning->append("The following search settings were used: ");
        std::string s1 = "protease = ";
        ProseCreatedWhileRunning->append( s1 +
                                          getCommonParameters()->getDigestionParams()->getProtease()->ToString() + "; ");
        std::string s2 = "maximum missed cleavages = ";
        ProseCreatedWhileRunning->append( s2 +
                                          std::to_string(getCommonParameters()->getDigestionParams()->getMaxMissedCleavages()) + "; ");
        std::string s3 = "minimum peptide length = ";
        ProseCreatedWhileRunning->append( s3 +
                                          std::to_string(getCommonParameters()->getDigestionParams()->getMinPeptideLength()) + "; ");
        std::string s4 = "maximum peptide length = unspecified; ";
        std::string s5 = "maximum peptide length = ";
        ProseCreatedWhileRunning->append(getCommonParameters()->getDigestionParams()->getMaxPeptideLength() ==
                                         std::numeric_limits<int>::max() ?
                                         s4 : s5 +
                                         std::to_string(getCommonParameters()->getDigestionParams()->getMaxPeptideLength()) + "; ");
        std::string s6 = "initiator methionine behavior = ";
        ProseCreatedWhileRunning->append( s6 +
                                          InitiatorMethionineBehaviorToString(getCommonParameters()->getDigestionParams()->getInitiatorMethionineBehavior()) + "; ");
        std::vector<std::string> vs;
        for (auto p: fixedModifications ) {
            vs.push_back(p->getIdWithMotif());
        }
        std::string comma = ", ";
        ProseCreatedWhileRunning->append("fixed modifications = " + StringHelper::join(vs, comma) + "; ");
        
        vs.clear();
        for (auto p: variableModifications ) {
            vs.push_back(p->getIdWithMotif());
        }
        std::string s7 = "variable modifications = ";
        ProseCreatedWhileRunning->append(s7 + StringHelper::join(vs, comma) + "; ");

        std::string s8 = "max mods per peptide = ";
        ProseCreatedWhileRunning->append( s8 +
                      std::to_string(getCommonParameters()->getDigestionParams()->getMaxModsForPeptide()) + "; ");
        std::string s9 = "max modification isoforms = ";
        ProseCreatedWhileRunning->append( s9 +
                      std::to_string(getCommonParameters()->getDigestionParams()->getMaxModificationIsoforms()) + "; ");
        std::string s10 = "precursor mass tolerance = ";
        ProseCreatedWhileRunning->append( s10 +
                                          std::to_string(getCommonParameters()->getPrecursorMassTolerance()->getValue()) + "; ");
        std::string s11 = "product mass tolerance = ";
        ProseCreatedWhileRunning->append( s11 +
                                          std::to_string(getCommonParameters()->getProductMassTolerance()->getValue()) + "; ");
        std::string s12 = "report PSM ambiguity = ";
        ProseCreatedWhileRunning->append(s12 +
                      StringHelper::toString(getCommonParameters()->getReportAllAmbiguity()) + ". ");
        int pL_iD_count=0, pL_iC_count=0;
        for ( auto p: proteinList ) {
            if ( p->getIsDecoy() ){
                pL_iD_count++;
            }
            if ( p->getIsContaminant() ) {
                pL_iC_count++;
            }
        }
        ProseCreatedWhileRunning->append("The combined search database contained " + std::to_string(pL_iD_count) +
                                         " non-decoy protein entries including " + std::to_string(pL_iC_count) +
                                         " contaminant sequences. ");
        // start the search task
        myTaskResults = new MyTaskResults(this);
        std::vector<PeptideSpectralMatch*> allPsms;
        
        //generate an array to store category specific fdr values (for speedy semi/nonspecific searches)
        //+1 because it starts at zero
        //int numFdrCategories = static_cast<int>(Enum::GetValues(FdrCategory::typeid)->Cast<FdrCategory>().Last() + 1); 

        //EDGAR: setting for now to the max value;
        int numFdrCategories = 3;
        std::vector<std::vector<PeptideSpectralMatch*>> allCategorySpecificPsms(numFdrCategories);
        for (int i = 0; i < numFdrCategories; i++)
        {
            allCategorySpecificPsms[i] = std::vector<PeptideSpectralMatch*>();
        }
        
        FlashLfqResults *flashLfqResults = nullptr;
        
        MyFileManager *myFileManager = new MyFileManager(getSearchParameters()->getDisposeOfFileWhenDone());
        
        std::vector<CommonParameters *> fileSpecificCommonParams;
        for ( auto b: fileSettingsList ) {
            fileSpecificCommonParams.push_back(SetAllFileSpecificCommonParams(getCommonParameters(), b));
        }
        
        int completedFiles = 0;
        //std::any indexLock = std::any();
        //std::any psmLock = std::any();
        
        Status("Searching files...", taskId, getVerbose());
        std::vector<std::string> vsx = {taskId, "Individual Spectra Files"};
        Status("Searching files...", vsx, getVerbose() );
        
        std::unordered_map<std::string, std::vector<int>> numMs2SpectraPerFile;
        for (int spectraFileIndex = 0; spectraFileIndex < (int)currentRawFileList.size(); spectraFileIndex++)
        {
            auto origDataFile = currentRawFileList[spectraFileIndex];
            
            // mark the file as in-progress
            std::vector<std::string> vs = {taskId, "Individual Spectra Files", origDataFile};
            StartingDataFile(origDataFile, vs, getVerbose() );
            
            EngineLayer::CommonParameters *combinedParams = SetAllFileSpecificCommonParams(getCommonParameters(),
                                                                             fileSettingsList[spectraFileIndex]);
            
            MassDiffAcceptor *massDiffAcceptor = GetMassDiffAcceptor(combinedParams->getPrecursorMassTolerance(),
                                                                     getSearchParameters()->getMassDiffAcceptorType(),
                                                                     getSearchParameters()->getCustomMdac());
            
            std::vector<std::string > thisId = {taskId, "Individual Spectra Files", origDataFile};
            //NewCollection(FileSystem::getFileName(origDataFile), thisId);
            Status("Loading spectra file...", thisId, getVerbose());
            MsDataFile *myMsDataFile = myFileManager->LoadFile(origDataFile,
                                                               std::make_optional(combinedParams->getTopNpeaks()),
                                                               std::make_optional(combinedParams->getMinRatio()),
                                                               combinedParams->getTrimMs1Peaks(),
                                                               combinedParams->getTrimMsMsPeaks(), combinedParams);
            Status("Getting ms2 scans...", thisId, getVerbose());
            std::vector<Ms2ScanWithSpecificMass*> arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile,
                                                                                            origDataFile,
                                                                                            combinedParams);
            std::sort(arrayOfMs2ScansSortedByMass.begin(), arrayOfMs2ScansSortedByMass.end(), [&]
                      (Ms2ScanWithSpecificMass *l, Ms2ScanWithSpecificMass *r) {
                          return l->getPrecursorMass() < r->getPrecursorMass();
                      });

            std::string fname = origDataFile.substr(0, origDataFile.find_last_of(".")); 
            int count=0;
            for ( auto p: myMsDataFile->GetAllScansList() ) {
                if ( p->getMsnOrder() == 2 ) {
                    count++;
                }
            }
            std::vector<int> v2= {count, (int)arrayOfMs2ScansSortedByMass.size() };
            numMs2SpectraPerFile.emplace(fname, v2);
            
            myFileManager->DoneWithFile(origDataFile);
            
            std::vector<PeptideSpectralMatch*> fileSpecificPsms(arrayOfMs2ScansSortedByMass.size());
            
            // modern search
            if (getSearchParameters()->getSearchType() == SearchType::Modern)
            {
                for (int currentPartition = 0; currentPartition < combinedParams->getTotalPartitions(); currentPartition++)
                {
                    std::vector<PeptideWithSetModifications*> peptideIndex;
                    int start = currentPartition * proteinList.size() / combinedParams->getTotalPartitions();
                    int count = ((currentPartition + 1) * proteinList.size() / combinedParams->getTotalPartitions()) -
                        (currentPartition * proteinList.size() / combinedParams->getTotalPartitions());

                    std::vector<Protein*> proteinListSubset;
                    for (auto i = 0; i < count; i++ ) {
                        proteinListSubset.push_back( proteinList[start+i] );
                    }

                    std::vector<std::string> vs2 = {taskId};
                    Status("Getting fragment dictionary...", vs2, getVerbose() );

                    std::vector<std::string> dbfnames;
                    for ( auto p: dbFilenameList ) {
                        dbfnames.push_back(p->getFilePath() );
                    }
                    auto indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications,
                                                          currentPartition, getSearchParameters()->getDecoyType(),
                                                          combinedParams, getSearchParameters()->getMaxFragmentSize(), false,
                                                          dbfnames, vs2, getVerbose());
                    
                    std::vector<std::vector<int>> fragmentIndex;
                    std::vector<std::vector<int>> precursorIndex;

                    //std::lock_guard<std::mutex> lock(indexLock);
                    auto tempmods = GlobalVariables::getAllModsKnown();
                    GenerateIndexes(indexEngine, dbFilenameList, peptideIndex, fragmentIndex, precursorIndex,
                                    proteinList, tempmods, taskId);
                    
                    Status("Searching files...", taskId, getVerbose());
                    
                    auto tempVar = new ModernSearchEngine (fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex,
                                               fragmentIndex, currentPartition, combinedParams, massDiffAcceptor,
                                               getSearchParameters()->getMaximumMassThatFragmentIonScoreIsDoubled(),
                                               thisId, getVerbose());
                    tempVar->Run();
                    
                    ProgressEventArgs tempVar2(100, "Done with search " + std::to_string(currentPartition + 1) + "/" +
                                               std::to_string(combinedParams->getTotalPartitions()) + "!", thisId);
                    ReportProgress(&tempVar2, getVerbose());
                    
                    delete indexEngine;
                    delete tempVar;
                }
            }
            // nonspecific search
            else if (getSearchParameters()->getSearchType() == SearchType::NonSpecific)
            {
                //generate an array of all possible locals
                std::vector<std::vector<PeptideSpectralMatch*>> fileSpecificPsmsSeparatedByFdrCategory(numFdrCategories); 
                for (int i = 0; i < numFdrCategories; i++) //only add if we're using for FDR, else ignore it as null.
                {
                    fileSpecificPsmsSeparatedByFdrCategory[i] = std::vector<PeptideSpectralMatch*>(arrayOfMs2ScansSortedByMass.size());
                }
                
                std::vector<EngineLayer::CommonParameters*> paramsToUse = {combinedParams};
                if (combinedParams->getDigestionParams()->getSearchModeType() == CleavageSpecificity::Semi) 
                {
                    //if semi, we need to do both N and C to hit everything
                    paramsToUse.clear();
                    std::vector<FragmentationTerminus> terminiToUse = {FragmentationTerminus::N, FragmentationTerminus::C};
                    for (auto terminus : terminiToUse) //set both termini
                    {
                        //paramsToUse.push_back(combinedParams->CloneWithNewTerminus(std::make_optional(terminus)));
                        auto c = new CommonParameters (combinedParams, std::make_optional(terminus), std::nullopt);
                        paramsToUse.push_back( c );
                    }
                }
                for (auto paramToUse : paramsToUse)
                {
                    for (int currentPartition = 0; currentPartition < paramToUse->getTotalPartitions(); currentPartition++)
                    {
                        std::vector<PeptideWithSetModifications*> peptideIndex;
                        
                        int start = currentPartition * proteinList.size() / paramToUse->getTotalPartitions();
                        int count = ((currentPartition + 1) * proteinList.size() / paramToUse->getTotalPartitions()) -
                            (currentPartition * proteinList.size() / paramToUse->getTotalPartitions());
                        std::vector<Protein*> proteinListSubset;
                        for ( auto i=0; i<count; i++ ) {
                            proteinListSubset.push_back(proteinList[start+i]);
                        }
                        
                        std::vector<std::vector<int>> fragmentIndex(1);
                        std::vector<std::vector<int>> precursorIndex(1);
                        
                        std::vector<std::string> vs3 = {taskId};
                        Status("Getting fragment dictionary...", vs3, getVerbose());
                        std::vector<std::string> dbfnames;
                        for ( auto p: dbFilenameList ) {
                            dbfnames.push_back(p->getFilePath() );
                        }

                        auto indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications,
                                                              currentPartition, getSearchParameters()->getDecoyType(),
                                                              paramToUse, getSearchParameters()->getMaxFragmentSize(), true,
                                                              dbfnames, vs3, getVerbose());
                        //std::lock_guard<std::mutex> lock(indexLock);
                        auto tempmods = GlobalVariables::getAllModsKnown();
                        GenerateIndexes(indexEngine, dbFilenameList, peptideIndex, fragmentIndex, precursorIndex,
                                        proteinList, tempmods, taskId);
                        
                        Status("Searching files...", taskId, getVerbose());
                        
                        auto tempVar3 = new NonSpecificEnzymeSearchEngine (fileSpecificPsmsSeparatedByFdrCategory,
                                                               arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex,
                                                               precursorIndex, currentPartition, paramToUse, massDiffAcceptor,
                                                               getSearchParameters()->getMaximumMassThatFragmentIonScoreIsDoubled(),
                                                               thisId, getVerbose() );
                        tempVar3->Run();
                        
                        ProgressEventArgs tempVar4(100, "Done with search " + std::to_string(currentPartition + 1) + "/" +
                                                   std::to_string(paramToUse->getTotalPartitions()) + "!", thisId);
                        ReportProgress(&tempVar4, getVerbose());
                         
                        delete indexEngine;
                        delete tempVar3;
                    }
                }
                {
                    //std::lock_guard<std::mutex> lock(psmLock);
                    for (int i = 0; i < (int)allCategorySpecificPsms.size(); i++)
                    {
                        if (allCategorySpecificPsms[i].size() > 0)
                        {
                            allCategorySpecificPsms[i].insert(allCategorySpecificPsms[i].end(),
                                                              fileSpecificPsmsSeparatedByFdrCategory[i].begin(),
                                                              fileSpecificPsmsSeparatedByFdrCategory[i].end());
                        }
                    }
                }
            }
            // classic search
            else
            {
                Status("Starting search...", thisId, getVerbose());
                auto tempVar5 = new ClassicSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, variableModifications,
                                                        fixedModifications, proteinList, massDiffAcceptor, combinedParams,
                                                        thisId, getVerbose());
                tempVar5->Run();
                
                ProgressEventArgs tempVar6(100, "Done with search!", thisId);
                ReportProgress(&tempVar6, getVerbose());
                delete tempVar5;
            }
            
            //std::lock_guard<std::mutex> lock(psmLock);
            allPsms.insert(allPsms.end(), fileSpecificPsms.begin(), fileSpecificPsms.end());
            
            completedFiles++;
            std::vector<std::string> vsy= {taskId, "Individual Spectra Files", origDataFile};
            FinishedDataFile(origDataFile, vsy, getVerbose() );
            std::vector<std::string> vs4 = {taskId, "Individual Spectra Files"};
            ProgressEventArgs tempVar7(completedFiles / currentRawFileList.size(), "Searching...", vs4);
            ReportProgress(&tempVar7, getVerbose());
        }

        std::vector<std::string> vs5 = {taskId, "Individual Spectra Files"};
        ProgressEventArgs tempVar8(100, "Done with all searches!", vs5);
        ReportProgress(&tempVar8, getVerbose());
        
        int numNotches = GetNumNotches(getSearchParameters()->getMassDiffAcceptorType(),
                                       getSearchParameters()->getCustomMdac());
        //resolve category specific fdrs (for speedy semi and nonspecific
        if (getSearchParameters()->getSearchType() == SearchType::NonSpecific)
        {
            allPsms = NonSpecificEnzymeSearchEngine::ResolveFdrCategorySpecificPsms(allCategorySpecificPsms,
                                                                                    numNotches, taskId, getCommonParameters());
        }
        
        PostSearchAnalysisParameters *parameters = new PostSearchAnalysisParameters();
        parameters->setSearchTaskResults(myTaskResults);
        parameters->setSearchTaskId(taskId);
        parameters->setSearchParameters(getSearchParameters());
        parameters->setProteinList(proteinList);
        parameters->setAllPsms(allPsms);
        parameters->setFixedModifications(fixedModifications);
        parameters->setVariableModifications(variableModifications);

        std::unordered_set<DigestionParams*> ptv;
        for ( auto p = fileSpecificCommonParams.begin(); p != fileSpecificCommonParams.end(); p++ ) {
            ptv.insert((*p)->getDigestionParams() );
        }
        parameters->setListOfDigestionParams(ptv);
        parameters->setCurrentRawFileList(currentRawFileList);
        parameters->setMyFileManager(myFileManager);
        parameters->setNumNotches(numNotches);
        parameters->setOutputFolder(OutputFolder);
        parameters->setIndividualResultsOutputFolder(OutputFolder + "/Individual File Results");
        parameters->setFlashLfqResults(flashLfqResults);
        parameters->setFileSettingsList(fileSettingsList);
        parameters->setNumMs2SpectraPerFile(numMs2SpectraPerFile);
        parameters->setDatabaseFilenameList(dbFilenameList);
        PostSearchAnalysisTask *postProcessing = new PostSearchAnalysisTask();
        postProcessing->setParameters(parameters);
        postProcessing->setCommonParameters(getCommonParameters());
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete parameters' statement was not added since parameters
        //was assigned to another object. Handle memory management manually.
        //C# TO C++ CONVERTER TODO TASK: A 'delete myFileManager' statement was not added since myFileManager
        //was assigned to another object. Handle memory management manually.
        return postProcessing->Run();
    }
    
    int SearchTask::GetNumNotches(MassDiffAcceptorType massDiffAcceptorType, const std::string &customMdac)
    {
        switch (massDiffAcceptorType)
        {
            case MassDiffAcceptorType::Exact:
                return 1;
            case MassDiffAcceptorType::OneMM:
                return 2;
            case MassDiffAcceptorType::TwoMM:
                return 3;
            case MassDiffAcceptorType::ThreeMM:
                return 4;
            case MassDiffAcceptorType::ModOpen:
                return 1;
            case MassDiffAcceptorType::Open:
                return 1;
            case MassDiffAcceptorType::Custom:
                return ParseSearchMode(customMdac)->getNumNotches();
                
            default:
                throw MetaMorpheusException("Unknown mass difference acceptor type");
        }
    }
    
    MassDiffAcceptor *SearchTask::ParseSearchMode(const std::string &text)
    {
        MassDiffAcceptor *massDiffAcceptor = nullptr;
        
        try
        {
            auto split = StringHelper::split(text, ' ');
            
            if (split[1] == "dot")
            {
#ifdef ORIG
                std::vector<double> massShifts = Array::ConvertAll(StringHelper::split(split[4], ','),
                                                                   double::Parse);
#endif
                auto stringvec = StringHelper::split(split[4], ',');
                std::vector<double> massShifts;
                for ( auto  p: stringvec ) {
                    massShifts.push_back(std::stod (p));
                }
                
                std::string nestring = StringHelper::replace(split[2], "�", "");
                double toleranceValue = std::stod(nestring);
                std::for_each(split[3].begin(), split[3].end(), [] (char &c)  { c = std::toupper(c);});
                if (split[3] == "PPM")
                {
                    auto  tempVar = new PpmTolerance (toleranceValue);
                    massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, tempVar);
                }
                else if (split[3] == "DA")
                {
                    auto  tempVar2 = new AbsoluteTolerance(toleranceValue);
                    massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, tempVar2);
                }
                
            }
            else if (split[1] == "interval")
            {
#ifdef ORIG
                std::vector<DoubleRange*> doubleRanges = Array::ConvertAll(StringHelper::split(split[2], ';'), [&] (std::any b)  {
                        new DoubleRange(std::stod(b->Trim(std::vector<char_t> {'[', ']'})->Split(',')[0], CultureInfo::InvariantCulture),
                                        std::stod(b->Trim(std::vector<char_t> {'[', ']'})->Split(',')[1], CultureInfo::InvariantCulture));
                    });
#endif
                std::vector<DoubleRange*> doubleRanges;
                std::vector<std::string> splitstring = StringHelper::split(split[2], ';');
                for ( auto b: splitstring ) {
                    StringHelper::replace(b, "[", "" );
                    StringHelper::replace(b, "]", "" );
                    std::vector<std::string> v = StringHelper::split(b, ',');
                    doubleRanges.push_back(new DoubleRange(std::stod(v[0]), std::stod(v[1])));
                }
                
                massDiffAcceptor = new IntervalMassDiffAcceptor(split[0], doubleRanges);
                
            }
            else if (split[1] == "OpenSearch")
            {
                massDiffAcceptor = new OpenSearchMode();
                
            }
            else if (split[1] == "daltonsAroundZero")
            {
                massDiffAcceptor = new SingleAbsoluteAroundZeroSearchMode(std::stod(split[2]));
                
            }
            else if (split[1] == "ppmAroundZero")
            {
                massDiffAcceptor = new SinglePpmAroundZeroSearchMode(std::stod(split[2]));
                
            }
            else
            {
                delete massDiffAcceptor;
                throw MetaMorpheusException("Unrecognized search mode type: " + split[1]);
            }
        }
        catch (const std::runtime_error &e)
        {
            delete massDiffAcceptor;
            std::string s = "Could not parse search mode string: ";
            throw MetaMorpheusException( s + e.what());
        }
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete massDiffAcceptor' statement was not added since
        //massDiffAcceptor was used in a 'return' or 'throw' statement.
        return massDiffAcceptor;
    }
}
