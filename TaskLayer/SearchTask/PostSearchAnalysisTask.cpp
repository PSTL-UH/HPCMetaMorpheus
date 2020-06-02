#include "PostSearchAnalysisTask.h"
#include "PostSearchAnalysisParameters.h"
#include "../../EngineLayer/ProteinParsimony/ProteinGroup.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../MyTaskResults.h"
#include "../../EngineLayer/GlobalVariables.h"
#include "../../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"
#include "../../EngineLayer/CommonParameters.h"
#include "../../EngineLayer/ProteinParsimony/ProteinParsimonyEngine.h"
#include "../../EngineLayer/ProteinParsimony/ProteinParsimonyResults.h"
#include "../../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.h"
#include "../../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrResults.h"
#include "../../EngineLayer/MetaMorpheusException.h"
#include "../../EngineLayer/FdrAnalysis/FdrAnalysisEngine.h"
#include "../../EngineLayer/ModificationAnalysis/ModificationAnalysisEngine.h"
#include "MzIdentMLWriter.h"
#include "../PepXMLWriter.h"

#include <fstream>
#include <iostream>
#include <algorithm>

#include <experimental/filesystem>
#include "Group.h"

using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::HistogramAnalysis;
using namespace EngineLayer::Localization;
using namespace EngineLayer::ModificationAnalysis;
using namespace FlashLFQ;
using namespace MassSpectrometry;
//using namespace MathNet::Numerics::Distributions;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{

    PostSearchAnalysisParameters *PostSearchAnalysisTask::getParameters() const
    {
        return privateParameters;
    }
    
    void PostSearchAnalysisTask::setParameters(PostSearchAnalysisParameters *value)
    {
        privateParameters = value;
    }
    
    std::vector<EngineLayer::ProteinGroup*> PostSearchAnalysisTask::getProteinGroups() const
    {
        return privateProteinGroups;
    }
    
    void PostSearchAnalysisTask::setProteinGroups(const std::vector<EngineLayer::ProteinGroup*> &value)
    {
        privateProteinGroups = value;
    }
    
    std::vector<std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>>*> PostSearchAnalysisTask::getPsmsGroupedByFile() const
    {
        return privatePsmsGroupedByFile;
    }
    
    
    void PostSearchAnalysisTask::setPsmsGroupedByFile(const std::vector<std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>>*> &value)    {
        privatePsmsGroupedByFile = value;
    }
    
    PostSearchAnalysisTask::PostSearchAnalysisTask() : MetaMorpheusTask(MyTask::Search)
    {
    }
    
    MyTaskResults *PostSearchAnalysisTask::Run()
    {
        // Stop loop if canceled
        if (GlobalVariables::getStopLoops())
        {
            return getParameters()->getSearchTaskResults();
        }
        
        if ( getParameters()->getSearchParameters()->getMassDiffAcceptorType() == MassDiffAcceptorType::ModOpen ||
             getParameters()->getSearchParameters()->getMassDiffAcceptorType() == MassDiffAcceptorType::Open ||
             getParameters()->getSearchParameters()->getMassDiffAcceptorType() == MassDiffAcceptorType::Custom)
        {
            // This only makes sense if there is a mass difference that you want to localize.
            //No use for exact and missed monoisotopic mass searches.
            getParameters()->getSearchParameters()->setDoLocalizationAnalysis(true);
        }
        else
        {
            getParameters()->getSearchParameters()->setDoLocalizationAnalysis(false);
        }
        
        //update all psms with peptide info
        //if it hasn't been done already
        if (getParameters()->getSearchParameters()->getSearchType() != SearchType::NonSpecific) 
        {
#ifdef ORIG
            getParameters()->setAllPsms(getParameters()->getAllPsms().Where([&] (std::any psm){
                        return psm != nullptr;
                    }).ToList());
#endif
            std::vector<PeptideSpectralMatch*> tmppsms;
            for ( auto psm: getParameters()->getAllPsms() ){
                if ( psm != nullptr ) {
                    tmppsms.push_back(psm);
                }
            }
            getParameters()->setAllPsms(tmppsms);
            
            for (auto psm : tmppsms )	{
                psm->ResolveAllAmbiguities();
            }
            
#ifdef ORIG            
            getParameters()->setAllPsms(getParameters()->getAllPsms().OrderByDescending([&] (std::any b) {
                        b::Score;
                    }).ThenBy([&] (std::any b)	{
                            b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
                        }).GroupBy([&] (std::any b) {
				(b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass);
                            })->Select([&] (std::any b)  {
                                    b::front();
                                }).ToList());
#endif
            //
            std::sort(tmppsms.begin(), tmppsms.end(), [&] (PeptideSpectralMatch *r, PeptideSpectralMatch *l) {
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
            for ( auto psm : tmppsms ) {
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
            std::vector<PeptideSpectralMatch*> scoreSorted;
            for ( auto t: tvec ) {
                scoreSorted.push_back(t[0]);
            }
            getParameters()->setAllPsms(scoreSorted);
            
            CalculatePsmFdr();
        }
        
        DoMassDifferenceLocalizationAnalysis();
        ProteinAnalysis();
        QuantificationAnalysis();
        
        std::vector<std::string> svec;
        svec.push_back(getParameters()->getSearchTaskId());
        svec.push_back("Individual Spectra Files");
        ProgressEventArgs tempVar(100, "Done!", svec);
        ReportProgress(&tempVar);
        
        HistogramAnalysis();
        WritePsmResults();
        WriteProteinResults();
        WriteQuantificationResults();
        WritePrunedDatabase();
        
        return getParameters()->getSearchTaskResults();
    }
    
    MyTaskResults *PostSearchAnalysisTask::RunSpecific(const std::string &OutputFolder,
                                                       std::vector<DbForTask*> &dbFilenameList,
                                                       std::vector<std::string> &currentRawFileList,
                                                       const std::string &taskId,
                                                       std::vector<FileSpecificParameters*> &fileSettingsList)
    {
        return nullptr;
    }
    
    void PostSearchAnalysisTask::CalculatePsmFdr()
    {
        // TODO: because FDR is done before parsimony, if a PSM matches to a target and a decoy protein,
        // there may be conflicts between how it's handled in parsimony and the FDR engine here
        // for example, here it may be treated as a decoy PSM, where as in parsimony it will be determined
        // by the parsimony algorithm which is agnostic of target/decoy assignments
        // this could cause weird PSM FDR issues
        
        Status("Estimating PSM FDR...", getParameters()->getSearchTaskId());
        int massDiffAcceptorNumNotches = getParameters()->getNumNotches();
        std::vector<std::string> svec;
        svec.push_back(getParameters()->getSearchTaskId());
        auto tmppsms = getParameters()->getAllPsms();
        FdrAnalysisEngine tempVar( tmppsms, massDiffAcceptorNumNotches, getCommonParameters(), svec);
        (&tempVar)->Run();
        
        // sort by q-value because of group FDR stuff
        // e.g. multiprotease FDR, non/semi-specific protease, etc
#ifdef ORIG
        getParameters()->setAllPsms(getParameters()->getAllPsms().OrderBy([&] (std::any p) {
                    p::FdrInfo::QValue;
		}).ThenByDescending([&] (std::any p) {
			p::Score;
                    }).ThenBy([&] (std::any p) {
                            p::FdrInfo::CumulativeTarget;
                        }).ToList());
#endif
        tmppsms.clear();
        tmppsms = getParameters()->getAllPsms();
        std::sort(tmppsms.begin(), tmppsms.end(), [&] (PeptideSpectralMatch *r, PeptideSpectralMatch *l) {
                if ( r->getFdrInfo()->getQValue() <  l->getFdrInfo()->getQValue() ) return true;
                if ( r->getFdrInfo()->getQValue() >  l->getFdrInfo()->getQValue() ) return false;
                
                if ( r->getScore()  > l->getScore() ) return true;
                if ( r->getScore()  < l->getScore() ) return false;
                
                if ( r->getFdrInfo()->getCumulativeTarget() <  l->getFdrInfo()->getCumulativeTarget() ) return true;
                return false;
            });
        getParameters()->setAllPsms(tmppsms);
        
        Status("Done estimating PSM FDR!", getParameters()->getSearchTaskId());
    }
    
    void PostSearchAnalysisTask::ProteinAnalysis()
    {
        if (!getParameters()->getSearchParameters()->getDoParsimony())
        {
            return;
        }
        
        Status("Constructing protein groups...", getParameters()->getSearchTaskId());
        
        // run parsimony
        auto tmppsms = getParameters()->getAllPsms();
        std::vector<std::string> svec = {getParameters()->getSearchTaskId()};
        ProteinParsimonyEngine tempVar(tmppsms, getParameters()->getSearchParameters()->getModPeptidesAreDifferent(),
                                       getCommonParameters(), svec );
        ProteinParsimonyResults *proteinAnalysisResults = static_cast<ProteinParsimonyResults*>((&tempVar)->Run());
        
        // score protein groups and calculate FDR
        auto tmp = proteinAnalysisResults->getProteinGroups();
        tmppsms = getParameters()->getAllPsms();
        ProteinScoringAndFdrEngine tempVar2(tmp, tmppsms,
                                            getParameters()->getSearchParameters()->getNoOneHitWonders(),
                                            getParameters()->getSearchParameters()->getModPeptidesAreDifferent(),
                                            true, getCommonParameters(), svec);
        ProteinScoringAndFdrResults *proteinScoringAndFdrResults = static_cast<ProteinScoringAndFdrResults*>((&tempVar2)->Run());
        
        setProteinGroups(proteinScoringAndFdrResults->SortedAndScoredProteinGroups);
        
        for (auto psm : getParameters()->getAllPsms())
        {
            psm->ResolveAllAmbiguities();
        }
        
        Status("Done constructing protein groups!", getParameters()->getSearchTaskId());
    }
    
    void PostSearchAnalysisTask::DoMassDifferenceLocalizationAnalysis()
    {
        if (getParameters()->getSearchParameters()->getDoLocalizationAnalysis())
        {
            Status("Running mass-difference localization analysis...", getParameters()->getSearchTaskId());
            for (int spectraFileIndex = 0; spectraFileIndex < (int)getParameters()->getCurrentRawFileList().size();
                 spectraFileIndex++)
            {
                EngineLayer::CommonParameters *combinedParams = SetAllFileSpecificCommonParams(getCommonParameters(),
                                                                                               getParameters()->getFileSettingsList()[spectraFileIndex]);
                
                auto origDataFile = getParameters()->getCurrentRawFileList()[spectraFileIndex];
                std::vector<std::string> svec = {getParameters()->getSearchTaskId(), "Individual Spectra Files",
                                                 origDataFile};
                Status("Running mass-difference localization analysis...", svec);
                MsDataFile *myMsDataFile = getParameters()->getMyFileManager()->LoadFile(origDataFile,
                                                                                         std::make_optional(combinedParams->getTopNpeaks()),
                                                                                         std::make_optional(combinedParams->getMinRatio()),
                                                                                         combinedParams->getTrimMs1Peaks(),
                                                                                         combinedParams->getTrimMsMsPeaks(), combinedParams);
                
                std::vector<std::string> svec2 = {getParameters()->getSearchTaskId(), "Individual Spectra Files",
                                                  origDataFile};
#ifdef ORIG
                LocalizationEngine tempVar(getParameters()->getAllPsms().Where([&] (std::any b) {
                            b::FullFilePath->Equals(origDataFile);
                        }).ToList(), myMsDataFile, combinedParams,svec2 );
#endif
                std::vector<PeptideSpectralMatch*> tmppsms;
                for ( auto b : getParameters()->getAllPsms() ) {
                    if ( b->getFullFilePath() == origDataFile ) {
                        tmppsms.push_back(b);
                    }
                }
                LocalizationEngine tempVar(tmppsms, myMsDataFile, combinedParams, svec2);
                (&tempVar)->Run();
                getParameters()->getMyFileManager()->DoneWithFile(origDataFile);
                std::vector<std::string> svec3 = {getParameters()->getSearchTaskId(), "Individual Spectra Files",
                                                  origDataFile};
                ProgressEventArgs tempVar2(100, "Done with localization analysis!", svec3);
                ReportProgress(&tempVar2);
            }
        }
        
        // count different modifications observed
        auto tmppsms2 = getParameters()->getAllPsms();
        std::vector<std::string> svec4 = {getParameters()->getSearchTaskId()};
        auto tempVar3 = new ModificationAnalysisEngine(tmppsms2, getCommonParameters(), svec4);
        tempVar3->Run();                 
    }
    
    void PostSearchAnalysisTask::QuantificationAnalysis()
    {
        if (!getParameters()->getSearchParameters()->getDoQuantification())
        {
            return;
        }
        
        // pass quantification parameters to FlashLFQ
        Status("Quantifying...", getParameters()->getSearchTaskId());
        
        // construct file info for FlashLFQ
        auto spectraFileInfo = std::vector<SpectraFileInfo*>();
        
        // get experimental design info for normalization
        if (getParameters()->getSearchParameters()->getNormalize())
        {
            std::experimental::filesystem::path Path = getParameters()->getCurrentRawFileList().front();//.FullName;
            std::string assumedExperimentalDesignPath = Path.parent_path();
            assumedExperimentalDesignPath = assumedExperimentalDesignPath +"/" +
                GlobalVariables::getExperimentalDesignFileName();
            
            if (std::experimental::filesystem::exists(assumedExperimentalDesignPath))
            {
#ifdef ORIG
                auto experimentalDesign = File::ReadAllLines(assumedExperimentalDesignPath).ToDictionary([&] (std::any p){
                        p->Split('\t')[0];
                    }, [&] (std::any p)  {
                        return p;
                    });
#endif
                std::unordered_map<std::string, std::string> experimentalDesign;
                char delimiter = '\t';
                std::ifstream input(assumedExperimentalDesignPath);
                if ( input.is_open() ) {
                    std::string line;
                    while ( getline( input, line) ) {
                        auto svec = StringHelper::split(line, delimiter);
                        for ( auto s: svec ) {
                            experimentalDesign.emplace(s, line);
                        }
                    }
                }
                else {
                    std::cout << "Could not open file " << assumedExperimentalDesignPath << std::endl;
                }
                input.close();
                
                for (auto file : getParameters()->getCurrentRawFileList())
                {
                    //std::string filename = Path::GetFileNameWithoutExtension(file);
                    std::string filename = file.substr(0, file.find_last_of("."));
                    
                    auto expDesignForThisFile = experimentalDesign[filename];
                    char delimiter = '\t';
                    auto split = StringHelper::split(expDesignForThisFile, delimiter);
                    
                    std::string condition = split[1];
                    int biorep = std::stoi(split[2]);
                    int fraction = std::stoi(split[3]);
                    int techrep = std::stoi(split[4]);
                    
                    // experimental design info passed in here for each spectra file
                    auto tempVar = new SpectraFileInfo( file, condition, biorep - 1, techrep - 1, fraction - 1 );
                    spectraFileInfo.push_back(tempVar);
                    
                    getParameters()->getMyFileManager()->DoneWithFile(file);
                }
            }
            else
            {
                throw MetaMorpheusException("Could not find experimental design file at location:\n" +
                                            assumedExperimentalDesignPath);
            }
        }
        else
        {
            for (auto file : getParameters()->getCurrentRawFileList())
            {
                // experimental design info passed in here for each spectra file
                auto tempVar2 = new SpectraFileInfo(file, "", 0, 0, 0);
                spectraFileInfo.push_back(tempVar2);
                getParameters()->getMyFileManager()->DoneWithFile(file);
            }
        }
        
        // get PSMs to pass to FlashLFQ
#ifdef ORIG
        auto unambiguousPsmsBelowOnePercentFdr = getParameters()->getAllPsms().Where([&] (std::any p)   {
                return p::FdrInfo::QValue <= 0.01 && p::FdrInfo::QValueNotch <= 0.01 &&
                !p::IsDecoy && p::FullSequence != nullptr;
            }).ToList();
#endif
        std::vector<PeptideSpectralMatch*> unambiguousPsmsBelowOnePercentFdr;
        for ( auto p = getParameters()->getAllPsms().begin(); p != getParameters()->getAllPsms().end(); p++ ) {
            if ( (*p)->getFdrInfo()->getQValue() <= 0.01      &&
                 (*p)->getFdrInfo()->getQValueNotch() <= 0.01 &&
                 !(*p)->getIsDecoy()                          &&
                 (*p)->getFullSequence().length() != 0 ) {
                unambiguousPsmsBelowOnePercentFdr.push_back(*p);
            }
        }
        
        
#ifdef ORIG
        auto psmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy([&] (std::any p)  {
                p::FullFilePath;
            });
#endif
        std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f1 = [&]
            (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
            return l->getFullFilePath() < r->getFullFilePath(); } ;
        std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f2 = [&]
            (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
            return l->getFullFilePath() != r->getFullFilePath(); } ;
        std::vector<std::vector<PeptideSpectralMatch*>> psmsGroupedByFile = Group::GroupBy ( unambiguousPsmsBelowOnePercentFdr, f1, f2);
        
        // pass protein group info for each PSM
        std::unordered_map<PeptideSpectralMatch*, std::vector<FlashLFQ::ProteinGroup*>> psmToProteinGroups;
        if (getProteinGroups().size() > 0)
        {
            for (auto proteinGroup : getProteinGroups())
            {
#ifdef ORIG
                auto proteinsOrderedByAccession = proteinGroup->getProteins().OrderBy([&] (std::any p)   {
                        p::Accession;
                    });
#endif
                auto proteinsOrderedByAccession = proteinGroup->getProteins();
                std::sort(proteinsOrderedByAccession.begin(), proteinsOrderedByAccession.end(), [&]
                          (Protein *l, Protein *r) {
                              return l->getAccession() < r->getAccession();
                          });
                
#ifdef ORIG
                auto flashLfqProteinGroup = new FlashLFQ::ProteinGroup(proteinGroup->getProteinGroupName(),
                                                                       std::string::Join("|", proteinsOrderedByAccession->Select([&] (std::any p) {
                                                                                   p::GeneNames->Select([&] (std::any x) {
                                                                                           x::Item2;
                                                                                       }).FirstOrDefault();
                                                                               })), std::string::Join("|", proteinsOrderedByAccession->Select([&] (std::any p){
                                                                                           p::Organism;
                                                                                       }).Distinct()));
#endif
                std::string del = "|";
                std::vector<std::string> svec1, svec2;
                for ( auto p: proteinsOrderedByAccession ) {
                    svec1.push_back(std::get<1>(p->getGeneNames().front()) );
                    bool found  = false;
                    for ( auto q: svec2 ) {
                        if ( q == p->getOrganism() ) {
                            found = true;
                            break;
                        }
                    }
                    if (!found ) {
                        svec2.push_back(p->getOrganism() );
                    }                    
                }
                std::string s1 = StringHelper::join (svec1, del);
                std::string s2 = StringHelper::join (svec2, del);
                auto flashLfqProteinGroup = new FlashLFQ::ProteinGroup(proteinGroup->getProteinGroupName(),
                                                                       s1, s2 );
                
#ifdef ORIG
                //for (auto psm : proteinGroup->getAllPsmsBelowOnePercentFDR().Where([&] (std::any v)  {
                //            return v::FullSequence != nullptr;
                //        }))
#endif
                for (auto psm : proteinGroup->getAllPsmsBelowOnePercentFDR() )  
                {
                    if ( psm->getFullSequence().length() == 0 ) {
                        continue;
                    }
                    
                    std::vector<FlashLFQ::ProteinGroup*> flashLfqProteinGroups;
                    std::unordered_map<PeptideSpectralMatch*, std::vector<FlashLFQ::ProteinGroup*>>::const_iterator psmToProteinGroups_iterator = psmToProteinGroups.find(psm);
                    if (psmToProteinGroups_iterator != psmToProteinGroups.end())
                    {
                        flashLfqProteinGroups = psmToProteinGroups_iterator->second;
                        flashLfqProteinGroups.push_back(flashLfqProteinGroup);
                    }
                    else
                    {
                        //flashLfqProteinGroups = psmToProteinGroups_iterator->second;
                        std::vector<FlashLFQ::ProteinGroup*> svec1 = {flashLfqProteinGroup};
                        psmToProteinGroups.emplace(psm, svec1);
                    }
                }
                
                //C# TO C++ CONVERTER TODO TASK: A 'delete flashLfqProteinGroup' statement was not added since
                //flashLfqProteinGroup was passed to a method or constructor. Handle memory management manually.
            }
        }
        else
        {
            // if protein groups were not constructed, just use accession numbers
            auto accessionToPg = std::unordered_map<std::string, FlashLFQ::ProteinGroup*>();
            for (auto psm : unambiguousPsmsBelowOnePercentFdr)
            {
#ifdef ORIG
                auto proteins = psm.BestMatchingPeptides->Select([&] (std::any b){
                        b::Peptide::Protein;
                    }).Distinct();
#endif
                std::vector<Protein*> proteins;
                for ( auto b: psm->getBestMatchingPeptides() ) {
                    bool found = false;
                    for ( auto q: proteins ) {
                        if ( q->Equals(std::get<1>(b)->getProtein()) ) {
                            found = true;
                            break;
                        }
                    }
                    if ( !found ) {
                        proteins.push_back(std::get<1>(b)->getProtein() );
                    }
                }
                
                for (auto protein : proteins)
                {
                    if (accessionToPg.find(protein->getAccession()) == accessionToPg.end())
                    {
                        std::vector<std::string> svec;
                        for ( auto p:  protein->getGeneNames() ) {
                            bool found = false;
                            std::string tmps = std::get<1>(p);
                            for ( auto q: svec ) {
                                if ( tmps == q ) {
                                    found = true;
                                    break;
                                }
                            }
                            if ( !found ) {
                                svec.push_back(tmps);
                            }
                        }
                        std::string del = "|";
                        std::string s = StringHelper::join(svec, del );
                        FlashLFQ::ProteinGroup tempVar3(protein->getAccession(), s, protein->getOrganism());
                        accessionToPg.emplace(protein->getAccession(), &tempVar3);
                    }
                    
                    std::vector<FlashLFQ::ProteinGroup*> proteinGroups;
                    std::unordered_map<PeptideSpectralMatch*, std::vector<FlashLFQ::ProteinGroup*>>::const_iterator psmToProteinGroups_iterator = psmToProteinGroups.find(psm);
                    if (psmToProteinGroups_iterator != psmToProteinGroups.end())
                    {
                        proteinGroups = psmToProteinGroups_iterator->second;
                        proteinGroups.push_back(accessionToPg[protein->getAccession()]);
                    }
                    else
                    {
                        //proteinGroups = psmToProteinGroups_iterator->second;
                        //psmToProteinGroups.emplace(psm, std::vector<FlashLFQ::ProteinGroup*>(protein->getAccession()) );
                        std::vector<FlashLFQ::ProteinGroup*>vpg = {accessionToPg[protein->getAccession()]};
                        psmToProteinGroups.emplace(psm, vpg );
                    }
                }
            }
        }
        
        // some PSMs may not have protein groups (if 2 peptides are required to construct a protein group,
        // some PSMs will be left over)
        // the peptides should still be quantified but not considered for protein quantification
        auto undefinedPg = new FlashLFQ::ProteinGroup("UNDEFINED", "", "");
        //sort the unambiguous psms by protease to make MBR compatible with multiple proteases
        std::unordered_map<Protease*, std::vector<PeptideSpectralMatch*>> proteaseSortedPsms;
        std::unordered_map<Protease*, FlashLfqResults*> proteaseSortedFlashLFQResults;
        
        for (auto dp : getParameters()->getListOfDigestionParams())
        {
            if (proteaseSortedPsms.find(dp->getProtease()) == proteaseSortedPsms.end())
            {
                proteaseSortedPsms.emplace(dp->getProtease(), std::vector<PeptideSpectralMatch*>());
            }
        }
        for (auto psm : unambiguousPsmsBelowOnePercentFdr)
        {
            if (psmToProteinGroups.find(psm) == psmToProteinGroups.end())
            {
                psmToProteinGroups.emplace(psm, std::vector<FlashLFQ::ProteinGroup*> {undefinedPg});
            }
            
            proteaseSortedPsms[psm->digestionParams->getProtease()].push_back(psm);
        }
        
        // pass PSM info to FlashLFQ
        auto flashLFQIdentifications = std::vector<Identification*>();
        for (auto spectraFile : psmsGroupedByFile)
        {
#ifdef ORIG
            auto rawfileinfo = spectraFileInfo.Where([&] (std::any p) {
                    p::FullFilePathWithExtension->Equals(spectraFile->Key);
                }).front();
#endif
            SpectraFileInfo* rawfileinfo;
            for ( auto p: spectraFileInfo ) {
                if ( p->FullFilePathWithExtension == spectraFile[0]->getFullFilePath() ) {
                    rawfileinfo = p;
                    break;
                }
            }
            
            for (auto psm : spectraFile)
            {
                auto tempVar4 = new Identification (rawfileinfo, psm->getBaseSequence(),
                                                    psm->getFullSequence(),
                                                    psm->getPeptideMonisotopicMass().value(),
                                                    psm->getScanRetentionTime(),
                                                    psm->getScanPrecursorCharge(),
                                                    psmToProteinGroups[psm]);
                flashLFQIdentifications.push_back(tempVar4);
            }
        }
        
        // run FlashLFQ
        auto flashLfqEngine = new FlashLfqEngine(flashLFQIdentifications,
                                                 getParameters()->getSearchParameters()->getNormalize(),
                                                 false, 
                                                 getParameters()->getSearchParameters()->getMatchBetweenRuns(),
                                                 getParameters()->getSearchParameters()->getQuantifyPpmTol(),
                                                 5.0, 5.0, false, 2, false, true,
                                                 true,
                                                 GlobalVariables::getElementsLocation(),
                                                 getCommonParameters()->getMaxThreadsToUsePerFile());
        
        if (!flashLFQIdentifications.empty())
        {
            getParameters()->setFlashLfqResults(flashLfqEngine->Run());
        }
        
        //MultiProtease MBR capability code
        //Parameters.FlashLfqResults = null;
        
        // EDGAR: NOte: this block was already commented out in the C# version
        //foreach (var proteasePsms in proteaseSortedPsms)
        //{
        //    var flashLFQIdentifications = new List<Identification>();
        //    var proteasePsmsGroupedByFile = proteasePsms.Value.GroupBy(p => p.FullFilePath);
        //    foreach (var spectraFile in proteasePsmsGroupedByFile)
        //    {
        //        var rawfileinfo = spectraFileInfo.Where(p => p.FullFilePathWithExtension.Equals(spectraFile.Key)).First();
        
        //        foreach (var psm in spectraFile)
        //        {
        //            flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence,
        //                psm.PeptideMonisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroups[psm]));
        //        }
        //    }
        
        //    // run FlashLFQ
        //    var FlashLfqEngine = new FlashLFQEngine(
        //        allIdentifications: flashLFQIdentifications,
        //        normalize: Parameters.SearchParameters.Normalize,
        //        ppmTolerance: Parameters.SearchParameters.QuantifyPpmTol,
        //        matchBetweenRuns: Parameters.SearchParameters.MatchBetweenRuns,
        //        silent: true,
        //        optionalPeriodicTablePath: GlobalVariables.ElementsLocation);
        
        //    if (flashLFQIdentifications.Any())
        //    {
        //        //make specific to protease
        //        var results = FlashLfqEngine.Run();
        
        //        if (Parameters.FlashLfqResults == null)
        //        {
        //            Parameters.FlashLfqResults = results;
        //        }
        //        else
        //        {
        //            Parameters.FlashLfqResults.MergeResultsWith(results);
        //        }
        //    }
        //}
        
        // get protein intensity back from FlashLFQ
        if (getProteinGroups().size() > 0 && getParameters()->getFlashLfqResults() != nullptr)
        {
            for (auto proteinGroup : getProteinGroups())
            {
                proteinGroup->setFilesForQuantification(spectraFileInfo);
                proteinGroup->setIntensitiesByFile(std::unordered_map<SpectraFileInfo*, double>());
                
                for (auto spectraFile : proteinGroup->getFilesForQuantification())
                {
#ifdef ORIG
                    //std::any flashLfqProteinGroup;
                    //if (getParameters()->getFlashLfqResults()->ProteinGroups.TryGetValue(proteinGroup->getProteinGroupName(), flashLfqProteinGroup))
#endif
                    auto pg = getParameters()->getFlashLfqResults()->ProteinGroups;
                    
                    if ( pg.find(proteinGroup->getProteinGroupName()) != pg.end() )
                    {
                        auto flashLfqProteinGroup = pg[proteinGroup->getProteinGroupName()];
                        proteinGroup->getIntensitiesByFile().emplace(spectraFile, flashLfqProteinGroup->GetIntensity(spectraFile));
                    }
                    else
                    {
                        proteinGroup->getIntensitiesByFile().emplace(spectraFile, 0);
                    }
                }
            }
        }
        
        delete flashLfqEngine;
        //C# TO C++ CONVERTER TODO TASK: A 'delete undefinedPg' statement was not added since undefinedPg
        //was passed to a method or constructor. Handle memory management manually.
    }
    
    void PostSearchAnalysisTask::HistogramAnalysis()
    {
        if (getParameters()->getSearchParameters()->getDoHistogramAnalysis())
        {
#ifdef ORIG
            //auto limitedpsms_with_fdr = getParameters()->getAllPsms().Where([&] (std::any b) {
            //        (b::FdrInfo::QValue <= 0.01);
            //    }).ToList();
#endif
            std::vector<PeptideSpectralMatch*> limitedpsms_with_fdr;
            for ( auto p = getParameters()->getAllPsms().begin(); p != getParameters()->getAllPsms().end(); p++ ) {
                if ( (*p)->getFdrInfo()->getQValue() <= 0.01 ) {
                    limitedpsms_with_fdr.push_back(*p);
                }
            }
            
#ifdef ORIG
            //if (limitedpsms_with_fdr.Any([&] (std::any b )  {
            //            !b::IsDecoy;
            //        }))
#endif
            bool any_cond = false;
            for ( auto b: limitedpsms_with_fdr ) {
                if ( !b->getIsDecoy() ) {
                    any_cond=true;
                    break;
                }
            }                                  
            
            if ( any_cond) 
            {
                std::vector<std::string> svec = {getParameters()->getSearchTaskId()};
                Status("Running histogram analysis...", svec );
                auto myTreeStructure = new BinTreeStructure();
                myTreeStructure->GenerateBins(limitedpsms_with_fdr,
                                              getParameters()->getSearchParameters()->getHistogramBinTolInDaltons());
                auto writtenFile = getParameters()->getOutputFolder() + "/MassDifferenceHistogram.tsv";
                WriteTree(myTreeStructure, writtenFile);
                FinishedWritingFile(writtenFile, svec);
                
                //C# TO C++ CONVERTER TODO TASK: A 'delete myTreeStructure' statement was not added since
                //myTreeStructure was passed to a method or constructor. Handle memory management manually.
            }
        }
    }
    
    void PostSearchAnalysisTask::WritePsmResults()
    {
        Status("Writing results...", getParameters()->getSearchTaskId());
#ifdef ORIG
        std::vector<PeptideSpectralMatch*> filteredPsmListForOutput = getParameters()->getAllPsms().Where([&] (std::any p) {
                return p::FdrInfo::QValue <= getCommonParameters()->getQValueOutputFilter() && p::FdrInfo::QValueNotch <= getCommonParameters()->getQValueOutputFilter();
            }).ToList();
#endif
        std::vector<PeptideSpectralMatch*> filteredPsmListForOutput;
        for ( auto p = getParameters()->getAllPsms().begin(); p != getParameters()->getAllPsms().end(); p++ ) {
            if ( (*p)->getFdrInfo()->getQValue() <= getCommonParameters()->getQValueOutputFilter()      &&
                 (*p)->getFdrInfo()->getQValueNotch() <= getCommonParameters()->getQValueOutputFilter() ) {
                filteredPsmListForOutput.push_back(*p);
            }
        }
        
        
        if (!getParameters()->getSearchParameters()->getWriteDecoys())
        {
#ifdef ORIG
            filteredPsmListForOutput.RemoveAll([&] (std::any b) {
                    b::IsDecoy;                   
                });
#endif
            std::remove_if (filteredPsmListForOutput.begin(), filteredPsmListForOutput.end(), [&] (PeptideSpectralMatch* b)
                            {return b->getIsDecoy();} );
        }
        
        if (!getParameters()->getSearchParameters()->getWriteContaminants())
        {
#ifdef ORIG
            filteredPsmListForOutput.RemoveAll([&] (std::any b)   {
                    b::IsContaminant;
                });
#endif
            std::remove_if (filteredPsmListForOutput.begin(), filteredPsmListForOutput.end(), [&] (PeptideSpectralMatch* b)
                            {return b->getIsContaminant();} );
            
        }
        
        // write PSMs
        std::string writtenFile = getParameters()->getOutputFolder() + "/AllPSMs.psmtsv";
        auto p = getParameters()->getSearchParameters()->getModsToWriteSelection();
        WritePsmsToTsv(filteredPsmListForOutput, writtenFile, &p);
        std::vector<std::string> tmpvecx = {getParameters()->getSearchTaskId()};
        FinishedWritingFile(writtenFile, tmpvecx);
        
        // write PSMs for percolator
        writtenFile = getParameters()->getOutputFolder() + "/AllPSMs_FormattedForPercolator.tsv";
        WritePsmsForPercolator(filteredPsmListForOutput, writtenFile, getCommonParameters()->getQValueOutputFilter());
        FinishedWritingFile(writtenFile, tmpvecx);
        
        // write best (highest-scoring) PSM per peptide
        writtenFile = getParameters()->getOutputFolder() +  "/AllPeptides.psmtsv";
        std::vector<PeptideSpectralMatch*> peptides = getParameters()->getAllPsms().GroupBy([&] (std::any b) {
                b::FullSequence;
            })->Select([&] (std::any b)  {
                    b::FirstOrDefault();
		}).ToList();
        WritePsmsToTsv(filteredPsmListForOutput.GroupBy([&] (std::any b)  {
                    b::FullSequence;
		})->Select([&] (std::any b) {
                        b::FirstOrDefault();
                    }).ToList(), writtenFile, getParameters()->getSearchParameters()->getModsToWriteSelection());
        std::vector<std::string> svec2 = {getParameters()->getSearchTaskId()};
        FinishedWritingFile(writtenFile, svec2);
        
        // write summary text
        getParameters()->getSearchTaskResults()->AddNiceText("All target PSMS within 1% FDR: " + getParameters()->getAllPsms().size()([&] (std::any a)    {
                    return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
                }));
        getParameters()->getSearchTaskResults()->AddNiceText("All target peptides within 1% FDR: " + peptides.size()([&] (std::any a){
                    return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
		}));
        if (getParameters()->getSearchParameters()->getDoParsimony())
        {
            getParameters()->getSearchTaskResults()->AddNiceText("All target protein groups within 1% FDR: " + std::to_string(getProteinGroups().size()([&] (std::any b) {
                            return b::QValue <= 0.01 && !b::IsDecoy;
			})) + "\r\n");
        }
        
        setPsmsGroupedByFile(filteredPsmListForOutput.GroupBy([&] (std::any p) {
                    p::FullFilePath;
		}));
        
        for (auto file : getPsmsGroupedByFile())
        {
            // write summary text
            auto psmsForThisFile = file; //->ToList();
            std::string fname = (*file->begin()).FullFilePath;
            std::string strippedFileName = fname.substr(0, fname.find_last_of(".")); 
            auto peptidesForFile = psmsForThisFile.GroupBy([&] (std::any b)  {
                    b::FullSequence;
                })->Select([&] (std::any b)   {
                        b::FirstOrDefault();
                    }).ToList();
            
            getParameters()->getSearchTaskResults()->AddNiceText("MS2 spectra in " + strippedFileName + ": " + std::to_string(getParameters()->getNumMs2SpectraPerFile()[strippedFileName][0]));
            getParameters()->getSearchTaskResults()->AddNiceText("Precursors fragmented in " + strippedFileName + ": " + std::to_string(getParameters()->getNumMs2SpectraPerFile()[strippedFileName][1]));
            getParameters()->getSearchTaskResults()->AddNiceText("Target PSMs within 1% FDR in " + strippedFileName + ": " + psmsForThisFile.size()([&] (std::any a) {
                        return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
                    }));
            getParameters()->getSearchTaskResults()->AddNiceText("Target peptides within 1% FDR in " +
                                                                 strippedFileName + ": " +
                                                                 std::to_string(peptidesForFile.size()([&] (std::any a)  {
                                                                             return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
                                                                         })) + "\r\n");
            
            // writes all individual spectra file search results to subdirectory
            if (getParameters()->getCurrentRawFileList().size() > 1)
            {
                // create individual files subdirectory
                FileSystem::createDirectory(getParameters()->getIndividualResultsOutputFolder());
                
                // write PSMs
                writtenFile = getParameters()->getIndividualResultsOutputFolder() + "/" + strippedFileName + "_PSMs.psmtsv";
                WritePsmsToTsv(psmsForThisFile, writtenFile, getParameters()->getSearchParameters()->getModsToWriteSelection());
                std::vector<std::string> tmpvec = {getParameters()->getSearchTaskId(),
                                                   "Individual Spectra Files", file->front().getFullFilePath()};
                FinishedWritingFile(writtenFile, tmpvec );
                
                // write PSMs for percolator
                writtenFile = getParameters()->getIndividualResultsOutputFolder() + "/" + strippedFileName +
                    "_PSMsFormattedForPercolator.tsv";
                WritePsmsForPercolator(psmsForThisFile, writtenFile, getCommonParameters()->getQValueOutputFilter());
                std::vector<std::string> tmpvec2 = {getParameters()->getSearchTaskId(),
                                                    "Individual Spectra Files", file->front().getFullFilePath()};
                FinishedWritingFile(writtenFile, tmpvec2 );
                
                // write best (highest-scoring) PSM per peptide
                writtenFile = getParameters()->getIndividualResultsOutputFolder() + "/" + strippedFileName + "_Peptides.psmtsv";
                WritePsmsToTsv(peptidesForFile, writtenFile, getParameters()->getSearchParameters()->getModsToWriteSelection());
                std::vector<std::string> tmpvec3 = {getParameters()->getSearchTaskId(), "Individual Spectra Files",
                                                    file->front().getFullFilePath()};
                FinishedWritingFile(writtenFile, tmpvec3);
            }
        }
    }
    
    void PostSearchAnalysisTask::WriteProteinResults()
    {
        if (getParameters()->getSearchParameters()->getDoParsimony())
        {
            // write protein groups to tsv
            std::string writtenFile = getParameters()->getOutputFolder() + "/AllProteinGroups.tsv";
            auto tvar = getProteinGroups();
            auto tvar2 = getCommonParameters()->getQValueOutputFilter();
            std::vector<std::string> tvec = {getParameters()->getSearchTaskId()}; 
            WriteProteinGroupsToTsv( tvar, writtenFile, tvec, tvar2 );
            
            // write all individual file results to subdirectory
            // local protein fdr, global parsimony, global psm fdr
            if (getParameters()->getCurrentRawFileList().size() > 1 ||
                getParameters()->getSearchParameters()->getWriteMzId() ||
                getParameters()->getSearchParameters()->getWritePepXml())
            {
                FileSystem::createDirectory(getParameters()->getIndividualResultsOutputFolder());
                
                for (auto fullFilePath : getPsmsGroupedByFile().Select([&] (std::any v)    {
                            v::Key;
                        }))
                {
                    std::string strippedFileName = Path::GetFileNameWithoutExtension(fullFilePath);
                    
                    std::vector<PeptideSpectralMatch*> psmsForThisFile = getPsmsGroupedByFile().Where([&] (std::any p)  {
                            return p->Key == fullFilePath;
                        }).SelectMany([&] (std::any g)   {
                                return g;
                            }).ToList();
                    auto subsetProteinGroupsForThisFile = getProteinGroups().Select([&] (std::any p)  {
                            p::ConstructSubsetProteinGroup(fullFilePath);
                        }).ToList();
                    
                    std::vector<std::string> tmpvec = {getParameters()->getSearchTaskId(),
                                                       "Individual Spectra Files", fullFilePath};
                    ProteinScoringAndFdrEngine tempVar(subsetProteinGroupsForThisFile, psmsForThisFile,
                                                       getParameters()->getSearchParameters()->getNoOneHitWonders(),
                                                       getParameters()->getSearchParameters()->getModPeptidesAreDifferent(),
                                                       false, getCommonParameters(),
                                                       tmpvec );
                    ProteinScoringAndFdrResults *subsetProteinScoringAndFdrResults = static_cast<ProteinScoringAndFdrResults*>((&tempVar)->Run());
                    
                    subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults->SortedAndScoredProteinGroups;
                    
                    getParameters()->getSearchTaskResults()->AddNiceText("Target protein groups within 1 % FDR in " + strippedFileName + ": " + subsetProteinGroupsForThisFile.size()([&] (std::any b)  {
                                return b::QValue <= 0.01 && !b::IsDecoy;
                            }));
                    
                    // write individual spectra file protein groups results to tsv
                    if (getParameters()->getCurrentRawFileList().size() > 1)
                    {
                        writtenFile =  getParameters()->getIndividualResultsOutputFolder() +"/" +
                            strippedFileName + "_ProteinGroups.tsv";
                        std::vector<std::string> tmpvec2;
                        tmpvec2.push_back(getParameters()->getSearchTaskId());
                        tmpvec2.push_back("Individual Spectra Files");
                        tmpvec2.push_back(fullFilePath);
                        WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, writtenFile, tmpvec2,
                                                getCommonParameters()->getQValueOutputFilter());
                    }
                    
                    // write mzID
                    if (getParameters()->getSearchParameters()->getWriteMzId())
                    {
                        std::vector<std::string> tmpvec2a;
                        tmpvec2a.push_back(getParameters()->getSearchTaskId());
                        tmpvec2a.push_back("Individual Spectra Files");
                        tmpvec2a.push_back(fullFilePath);
                        Status("Writing mzID...", tmpvec2a);
                        
                        auto mzidFilePath = getParameters()->getIndividualResultsOutputFolder() + "/" +
                            strippedFileName + ".mzID";
                        MzIdentMLWriter::WriteMzIdentMl(psmsForThisFile, subsetProteinGroupsForThisFile,
                                                        getParameters()->getVariableModifications(),
                                                        getParameters()->getFixedModifications(),
                                                        {getCommonParameters()->getDigestionParams()->getProtease()},
                                                        getCommonParameters()->getQValueOutputFilter(),
                                                        getCommonParameters()->getProductMassTolerance(),
                                                        getCommonParameters()->getPrecursorMassTolerance(),
                                                        getCommonParameters()->getDigestionParams()->getMaxMissedCleavages(),
                                                        mzidFilePath);
                        
                        FinishedWritingFile(mzidFilePath, tmpvec2a);
                    }
                    
                    // write pepXML
                    if (getParameters()->getSearchParameters()->getWritePepXml())
                    {
                        std::vector<std::string> svec;
                        svec.push_back(getParameters()->getSearchTaskId());
                        svec.push_back("Individual Spectra Files");
                        svec.push_back(fullFilePath);

                        Status("Writing pepXML...", svec);
                        
                        auto pepXMLFilePath = getParameters()->getIndividualResultsOutputFolder() + "/" +
                            strippedFileName + ".pep.XM";
                        auto tempvar = getParameters()->getDatabaseFilenameList();
                        auto tvar = getParameters()->getVariableModifications();
                        auto tvar2 = getParameters()->getFixedModifications();
                        PepXMLWriter::WritePepXml(psmsForThisFile, tempvar, tvar, tvar2,
                                                  getCommonParameters(), pepXMLFilePath,
                                                  getCommonParameters()->getQValueOutputFilter());
                        
                        FinishedWritingFile(pepXMLFilePath, svec );
                    }
                    
                    std::vector<std::string> tmpvec5;
                    tmpvec5.push_back(getParameters()->getSearchTaskId());
                    tmpvec5.push_back("Individual Spectra Files");
                    tmpvec5.push_back(fullFilePath);

                    ProgressEventArgs tempVar2(100, "Done!", tmpvec5 );
                    ReportProgress(&tempVar2);
                }
            }
        }
    }
    
    void PostSearchAnalysisTask::WriteQuantificationResults()
    {
        if (getParameters()->getSearchParameters()->getDoQuantification() &&
            getParameters()->getFlashLfqResults() != nullptr)
        {
            // write peaks
            std::vector<std::string> tmpvec = {getParameters()->getSearchTaskId()};
            WritePeakQuantificationResultsToTsv(getParameters()->getFlashLfqResults(),
                                                getParameters()->getOutputFolder(),
                                                "AllQuantifiedPeaks",
                                                tmpvec);
            
            // write peptide quant results
            WritePeptideQuantificationResultsToTsv(getParameters()->getFlashLfqResults(),
                                                   getParameters()->getOutputFolder(),
                                                   "AllQuantifiedPeptides",
                                                   tmpvec );
            
            // write individual results
            if (getParameters()->getCurrentRawFileList().size() > 1)
            {
                for (auto file : getParameters()->getFlashLfqResults()->Peaks)
                {
                    std::vector<std::string> vec3 = {getParameters()->getSearchTaskId(),
                                                     "Individual Spectra Files",
                                                     file->Key->FullFilePathWithExtension}
                        WritePeakQuantificationResultsToTsv(getParameters()->getFlashLfqResults(),
                                                            getParameters()->getIndividualResultsOutputFolder(),
                                                            file->Key->FilenameWithoutExtension + "_QuantifiedPeaks",
                                                            vec3);
                }
            }
        }
    }
    
    void PostSearchAnalysisTask::WritePrunedDatabase()
    {
        if (getParameters()->getSearchParameters()->getWritePrunedDatabase())
        {
            std::vector<std::string> tempvec = {getParameters()->getSearchTaskId()};
            Status("Writing Pruned Database...", tempvec);
            std::unordered_set<Modification*> modificationsToWriteIfBoth;
            std::unordered_set<Modification*> modificationsToWriteIfInDatabase;
            std::unordered_set<Modification*> modificationsToWriteIfObserved;
            
#ifdef ORIG
            auto confidentPsms = getParameters()->getAllPsms().Where([&] (std::any b)   {
                    return b::FdrInfo::QValueNotch <= 0.01 && b::FdrInfo::QValue <= 0.01 &&
                    !b::IsDecoy && b::BaseSequence != nullptr;
                }).ToList();
#endif
            std::vector<PeptideSpectralMatch* > confidentPsms;
            for ( auto b:  getParameters()->getAllPsms() ) {
                if ( b->getFdrInfo()->getQValueNotch() <= 0.01 &&
                     b->getFdrInfo()->getQValue()      <= 0.01 &&
                     !b->getIsDecoy                            &&
                     b->getBaseSequence().length() != 0 ){
                    confidentPsms.push_back(b);
                }
            }
            
            auto proteinToConfidentBaseSequences = std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>();
            
            // associate all confident PSMs with all possible proteins they could be digest products of (before or after parsimony)
            for (auto psm : confidentPsms)
            {
#ifdef ORIG
                auto myPepsWithSetMods = psm->BestMatchingPeptides->Select([&] (std::any p)  {
                        p::Peptide;
                    });
#endif
                std::vector<PeptideWithSetModifications *> myPepsWithSetMods;
                for ( auto p: psm->getBestMatchingPeptides() ) {
                    myPepsWithSetMods.push_back(std::get<1>(p) );
                }
                
                
                for (auto peptide : myPepsWithSetMods)
                {
                    std::vector<PeptideWithSetModifications*> myPepList;
                    std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>::const_iterator proteinToConfidentBaseSequences_iterator = proteinToConfidentBaseSequences.find(peptide.Protein.NonVariantProtein);
                    if (proteinToConfidentBaseSequences_iterator != proteinToConfidentBaseSequences.end())
                    {
                        myPepList = proteinToConfidentBaseSequences_iterator->second;
                        myPepList.push_back(peptide);
                    }
                    else
                    {
                        //myPepList = proteinToConfidentBaseSequences_iterator->second;
                        std::vector<PeptideWithSetModifications*> tvec = {peptide};
                        proteinToConfidentBaseSequences.emplace(peptide->Protein.NonVariantProtein, tvec);
                    }
                }
            }
            
            // Add user mod selection behavours to Pruned DB
            for (auto modType : getParameters()->getSearchParameters()->getModsToWriteSelection())
            {
                for (Modification *mod : GlobalVariables::getAllModsKnown().Where([&] (std::any b)    {
                            b::ModificationType->Equals(modType.Key);
                        }))
                {
                    if (modType.Value == 1) // Write if observed and in database
                    {
                        modificationsToWriteIfBoth.insert(mod);
                    }
                    if (modType.Value == 2) // Write if in database
                    {
                        modificationsToWriteIfInDatabase.insert(mod);
                    }
                    if (modType.Value == 3) // Write if observed
                    {
                        modificationsToWriteIfObserved.insert(mod);
                    }
                }
            }
            
            //generates dictionary of proteins with only localized modifications
#ifdef ORIG
            auto ModPsms = getParameters()->getAllPsms().Where([&] (std::any b) {
                    return b::FdrInfo::QValueNotch <= 0.01 && b::FdrInfo::QValue <= 0.01 &&
                    !b::IsDecoy && b::FullSequence != nullptr;
                }).ToList();
#endif
            std::vector<PeptideSpectralMatch* > ModPsms;
            for ( auto b:  getParameters()->getAllPsms() ) {
                if ( b->getFdrInfo()->getQValueNotch() <= 0.01 &&
                     b->getFdrInfo()->getQValue()      <= 0.01 &&
                     !b->getIsDecoy                            &&
                     b->getFullSequence().length() != 0 ){
                    ModPsms.push_back(b);
                }
            }
            
            std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>> proteinToConfidentModifiedSequences; 
            
            for (auto psm : ModPsms)
            {
#ifdef ORIG
                auto myPepsWithSetMods = psm->BestMatchingPeptides->Select([&] (std::any p)     {
                        p::Peptide;
                    });
#endif
                std::vector<PeptideWithSetModifications *> myPepsWithSetMods;
                for ( auto p: psm->getBestMatchingPeptides() ) {
                    myPepsWithSetMods.push_back( std::get<1>(p) );
                }
                
                for (auto peptide : myPepsWithSetMods)
                {
                    std::vector<PeptideWithSetModifications*> myPepList;
                    std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>::const_iterator proteinToConfidentModifiedSequences_iterator = proteinToConfidentModifiedSequences.find(peptide->getProtein()->getNonVariantProtein());
                    if (proteinToConfidentModifiedSequences_iterator != proteinToConfidentModifiedSequences.end())
                    {
                        myPepList = proteinToConfidentModifiedSequences_iterator->second;
                        myPepList.push_back(peptide);
                    }
                    else
                    {
                        //myPepList = proteinToConfidentModifiedSequences_iterator->second;
                        std::vector<PeptideWithSetModifications*> tvar = {peptide};
                        proteinToConfidentModifiedSequences.emplace(peptide->getProtein()->getNonVariantProtein(), tvar);
                    }
                }
            }
            
            // mods included in pruned database will only be confidently localized mods (peptide's FullSequence != null)
            for (auto nonVariantProtein : getParameters()->getProteinList().Select([&] (std::any p)  {
                        p::NonVariantProtein;
                    }).Distinct())
            {
                if (!nonVariantProtein::IsDecoy)
                {
                    std::vector<PeptideWithSetModifications*> psms;
                    std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>::const_iterator proteinToConfidentModifiedSequences_iterator = proteinToConfidentModifiedSequences.find(nonVariantProtein);
                    psms = proteinToConfidentModifiedSequences_iterator->second;
                    
                    // sequence variant is null if mod is not on a variant
                    std::unordered_set<(int, Modification, SequenceVariation)*> modsObservedOnThisProtein;
                    
                    for (auto psm : (psms.size() != 0 ) ? psms : std::vector<PeptideWithSetModifications*>())
                    {
                        for (auto idxModKV : psm->AllModsOneIsNterminus)
                        {
                            int proteinIdx = GetOneBasedIndexInProtein(idxModKV->Key, psm);
                            SequenceVariation *relevantVariant = psm->Protein.AppliedSequenceVariations.FirstOrDefault([&] (std::any sv) {
                                    VariantApplication::IsSequenceVariantModification(sv, proteinIdx);
                                });
                            SequenceVariation *unappliedVariant = relevantVariant == nullptr ? nullptr : psm->Protein.SequenceVariations.FirstOrDefault([&] (std::any sv)  {
                                    return sv::Description != nullptr &&
                                    sv::Description->Equals(relevantVariant->Description);
                                });
                            modsObservedOnThisProtein.insert((VariantApplication::RestoreModificationIndex(psm->Protein, proteinIdx), idxModKV->Value, unappliedVariant));
                        }
                    }
                    
                    std::unordered_map<(SequenceVariation, int)*, std::vector<Modification*>> modsToWrite;
                    
                    //Add if observed (regardless if in database)
                    for (auto observedMod : modsObservedOnThisProtein)
                    {
                        auto tempMod = observedMod->Item2;
                        
                        if (std::find(modificationsToWriteIfObserved.begin(), modificationsToWriteIfObserved.end(), tempMod) != modificationsToWriteIfObserved.end())
                        {
                            auto svIdxKey = (observedMod->Item3, observedMod->Item1);
                            if (modsToWrite.find(svIdxKey) == modsToWrite.end())
                            {
                                std::vector<Modification*> tvec = {observedMod->Item2};
                                modsToWrite.emplace(svIdxKey, tvec);
                            }
                            else
                            {
                                modsToWrite[svIdxKey].push_back(observedMod->Item2);
                            }
                        }
                    }
                    
                    // Add modification if in database (two cases: always or if observed)
                    for (auto modkv : nonVariantProtein::OneBasedPossibleLocalizedModifications)
                    {
                        for (auto mod : modkv->Value)
                        {
                            //Add if always In Database or if was observed and in database and not set to not include
                            if (std::find(modificationsToWriteIfInDatabase.begin(), modificationsToWriteIfInDatabase.end(), mod) != modificationsToWriteIfInDatabase.end() ||
                                (std::find(modificationsToWriteIfBoth.begin(), modificationsToWriteIfBoth.end(), mod) != modificationsToWriteIfBoth.end() &&
                                 std::find(modsObservedOnThisProtein.begin(), modsObservedOnThisProtein.end(), (modkv->Key, mod, nullptr)) != modsObservedOnThisProtein.end())))	{
                            if (modsToWrite.find((nullptr, modkv->Key)) == modsToWrite.end())
                            {
                                modsToWrite.emplace((nullptr, modkv->Key), std::vector<Modification*> {mod});
                            }
                            else
                            {
                                modsToWrite[(nullptr, modkv->Key)].push_back(mod);
                            }
                        }
                    }
                }
                
                // Add variant modification if in database (two cases: always or if observed)
                for (SequenceVariation *sv : nonVariantProtein::SequenceVariations)
                {
                    for (auto modkv : sv->OneBasedModifications)
                    {
                        for (auto mod : modkv->Value)
                        {
                            //Add if always In Database or if was observed and in database and not set to not include
                            if (std::find(modificationsToWriteIfInDatabase.begin(), modificationsToWriteIfInDatabase.end(), mod) != modificationsToWriteIfInDatabase.end() ||
                                (std::find(modificationsToWriteIfBoth.begin(), modificationsToWriteIfBoth.end(), mod) != modificationsToWriteIfBoth.end() &&
                                 std::find(modsObservedOnThisProtein.begin(), modsObservedOnThisProtein.end(), (modkv->Key, mod, sv)) != modsObservedOnThisProtein.end())))
                        {
                            if (modsToWrite.find((sv, modkv->Key)) == modsToWrite.end())
                            {
                                modsToWrite.emplace((sv, modkv->Key), std::vector<Modification*> {mod});
                            }
                            else
                            {
                                modsToWrite[(sv, modkv->Key)].push_back(mod);
                            }
                        }
                    }
                }
            }
            
            if (proteinToConfidentBaseSequences.find(nonVariantProtein::NonVariantProtein) != proteinToConfidentBaseSequences.end())
            {
                // adds confidently localized and identified mods
                nonVariantProtein::OneBasedPossibleLocalizedModifications->Clear();
                for (auto kvp : modsToWrite.Where([&] (std::any kv)  {
                            return kv::Key->Item1 == nullptr;
                        }))
                {
                    nonVariantProtein::OneBasedPossibleLocalizedModifications->Add(kvp::Key->Item2, kvp->Value);
                }
                for (auto sv : nonVariantProtein::SequenceVariations)
                {
                    sv->OneBasedModifications->Clear();
                    for (auto kvp : modsToWrite.Where([&] (std::any kv)  {
                                return kv::Key->Item1 != nullptr && kv::Key->Item1->Equals(sv);
                            }))
                    {
                        sv->OneBasedModifications->Add(kvp::Key->Item2, kvp->Value);
                    }
                }
            }
        }
        
        
        //writes all proteins
        if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)  {
                    !b::IsContaminant;
                }))
        {
            std::string outputXMLdbFullName = FileSystem::combine(getParameters()->getOutputFolder(), std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)  {
                            !b::IsContaminant;
                        })->Select([&] (std::any b)  {
                                Path::GetFileNameWithoutExtension(b::FilePath);
                            })) + "pruned.xml");
            ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(), getParameters()->getProteinList().Select([&] (std::any p)  {
                        p::NonVariantProtein;
                    }).Where([&] (std::any b) {
                            return !b::IsDecoy && !b::IsContaminant;
                        }).ToList(), outputXMLdbFullName);
            FinishedWritingFile(outputXMLdbFullName, std::vector<std::string> {getParameters()->getSearchTaskId()});
        }
        if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b) {
                    b::IsContaminant;
                }))
        {
            std::string outputXMLdbFullNameContaminants = FileSystem::combine(getParameters()->getOutputFolder(), std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)  {
                            b::IsContaminant;
                        })->Select([&] (std::any b) {
                                Path::GetFileNameWithoutExtension(b::FilePath);
                            })) + "pruned.xml");
            ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(), getParameters()->getProteinList().Select([&] (std::any p)  {
                        p::NonVariantProtein;
                    }).Where([&] (std::any b) {
                            return !b::IsDecoy && b::IsContaminant;
                        }).ToList(), outputXMLdbFullNameContaminants);
            std::vector<std::string> tvec = {getParameters()->getSearchTaskId()};
            FinishedWritingFile(outputXMLdbFullNameContaminants, tvec);
        }
        
        //writes only detected proteins
        if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)  {
                    !b::IsContaminant;
                }))
        {
            std::string outputXMLdbFullName = FileSystem::combine(getParameters()->getOutputFolder(),
                                                                  std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)   {
                                                                              !b::IsContaminant;
                                                                          })->Select([&] (std::any b) {
                                                                                  Path::GetFileNameWithoutExtension(b::FilePath);
                                                                              })) + "proteinPruned.xml");
            ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(), proteinToConfidentBaseSequences.Keys->Where([&] (std::any b) {
                        return !b::IsDecoy && !b::IsContaminant;
                    }).ToList(), outputXMLdbFullName);
            FinishedWritingFile(outputXMLdbFullName, std::vector<std::string> {getParameters()->getSearchTaskId()});
        }
        if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)  {
                    b::IsContaminant;
                }))
        {
            std::string outputXMLdbFullNameContaminants = FileSystem::combine(getParameters()->getOutputFolder(), std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b) {
                            b::IsContaminant;
                        })->Select([&] (std::any b)  {
                                Path::GetFileNameWithoutExtension(b::FilePath);
                            })) + "proteinPruned.xml");
            ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(), proteinToConfidentBaseSequences.Keys->Where([&] (std::any b)  {
                        return !b::IsDecoy && b::IsContaminant;
                    }).ToList(), outputXMLdbFullNameContaminants);
            FinishedWritingFile(outputXMLdbFullNameContaminants, std::vector<std::string> {getParameters()->getSearchTaskId()});
        }
        
    }
    
    int PostSearchAnalysisTask::GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications *peptideWithSetModifications)
    {
        if (oneIsNterminus == 1)
        {
            return peptideWithSetModifications->OneBasedStartResidueInProtein;
        }
        if (oneIsNterminus == peptideWithSetModifications->Length + 2)
        {
            return peptideWithSetModifications->OneBasedEndResidueInProtein;
        }
        return peptideWithSetModifications->OneBasedStartResidueInProtein + oneIsNterminus - 2;
    }
    
    void PostSearchAnalysisTask::WriteTree(BinTreeStructure *myTreeStructure, const std::string &writtenFile)
    {
        StreamWriter output = StreamWriter(writtenFile);
        output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tUnimodDiffs\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tprotNtermLocFrac\tpepNtermLocFrac\tpepCtermLocFrac\tprotCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot");
        for (Bin *bin : myTreeStructure->getFinalBins().OrderByDescending([&] (std::any b) {
                    b->Count;
                }))
        {
            output.WriteLine(bin->MassShift.ToString("F4") + "\t" +
                             bin->Count.ToString() + "\t" +
                             bin->CountDecoy.ToString() + "\t" +
                             bin->CountTarget.ToString() + "\t" +
                             bin->LocalizeableTarget.ToString() + "\t" +
                             (bin->CountTarget - bin->LocalizeableTarget).ToString() + "\t" +
                             (bin->Count == 0 ? NAN : static_cast<double>(bin->CountDecoy) / bin->Count).ToString("F3") + "\t" +
                             (Normal::CDF(0, 1, bin->ComputeZ(0.01))).ToString("F3") + "\t" +
                             (Normal::CDF(0, 1, bin->ComputeZ(0.255))).ToString("F3") + "\t" +
                             (bin->CountTarget == 0 ? NAN : static_cast<double>(bin->LocalizeableTarget) / bin->CountTarget).ToString("F3") + "\t" +
                             bin->Mine + "\t" + bin->UnimodId + "\t" + bin->UnimodFormulas + "\t" +
                             bin->UnimodDiffs + "\t" + bin->AA + "\t" + bin->Combos + "\t" +
                             std::string::Join(",", bin->ModsInCommon::OrderByDescending([&] (std::any b){
                                         b->Value;
                                     }).Where([&] (std::any b)	{
                                             return b->Value > bin->CountTarget / 10.0;
                                         })->Select([&] (std::any b) {
                                                 return b::Key + ":" + (static_cast<double>(b->Value) / bin->CountTarget).ToString("F3");
                                             })) + "\t" +
                             std::string::Join(",", bin->AAsInCommon::OrderByDescending([&] (std::any b){
                                         b->Value;
                                     }).Where([&] (std::any b) {
                                             return b->Value > bin->CountTarget / 10.0;
                                         })->Select([&] (std::any b)	{
                                                 return b::Key + ":" + (static_cast<double>(b->Value) / bin->CountTarget).ToString("F3");
                                             })) + "\t" +
                             std::string::Join(",", bin->ResidueCount::OrderByDescending([&] (std::any b) {
                                         b->Value;
                                     })->Select([&] (std::any b){
                                             return b::Key + ":" + b->Value;
                                         })) + "\t" +
                             (bin->LocalizeableTarget == 0 ? NAN : static_cast<double>(bin->ProtNlocCount) / bin->LocalizeableTarget).ToString("F3") + "\t" +
                             (bin->LocalizeableTarget == 0 ? NAN : static_cast<double>(bin->PepNlocCount) / bin->LocalizeableTarget).ToString("F3") + "\t" +
                             (bin->LocalizeableTarget == 0 ? NAN : static_cast<double>(bin->PepClocCount) / bin->LocalizeableTarget).ToString("F3") + "\t" + (bin->LocalizeableTarget == 0 ? NAN : static_cast<double>(bin->ProtClocCount) / bin->LocalizeableTarget).ToString("F3") + "\t" +
                             (bin->FracWithSingle).ToString("F3") + "\t" +
                             (static_cast<double>(bin->Overlapping) / bin->CountTarget).ToString("F3") + "\t" +
                             (bin->MedianLength).ToString("F3") + "\t" + bin->UniprotID);
        }
        
    }
    
    void PostSearchAnalysisTask::WritePsmsForPercolator(std::vector<PeptideSpectralMatch*> &psmList,
                                                        const std::string &writtenFileForPercolator,
                                                        double qValueCutoff)
    {
        StreamWriter output = StreamWriter(writtenFileForPercolator);
        output.WriteLine("SpecId\tLabel\tScanNr\tF1\tF2\tPeptide\tProteins");
        output.WriteLine("DefaultDirection\t-\t-\t1\t1\t\t");
        for (int i = 0; i < psmList.size(); i++)
        {
            auto psm = psmList[i];
            
            if (psm->getFdrInfo()->getQValue() > qValueCutoff || psm->getFdrInfo()->getQValueNotch() > qValueCutoff)
            {
                continue;
            }
            
            output.Write(std::to_string(i));
            output.Write('\t' + std::to_string(psm->getIsDecoy() ? -1 : 1));
            output.Write('\t' + std::to_string(psm->getScanNumber()));
            
            // Features
            output.Write(StringHelper::toString('\t') + std::string::Join("\t", psm->getFeatures()));
            
            // HACKY: Ignores all ambiguity
            auto pwsm = psm->BestMatchingPeptides.front().Peptide;
            
            output.Write('\t' + (pwsm->PreviousAminoAcid + "." + pwsm->FullSequence + "." + pwsm->NextAminoAcid).ToString());
            output.Write('\t' + (pwsm->Protein.Accession).ToString());
            output.WriteLine();
        }
    }
    
    
    void PostSearchAnalysisTask::WriteProteinGroupsToTsv(std::vector<EngineLayer::ProteinGroup*> &proteinGroups,
                                                         const std::string &filePath,
                                                         std::vector<std::string> &nestedIds,
                                                         double qValueCutoff)
    {
        if (proteinGroups.size() > 0 && proteinGroups.Any())
        {
            StreamWriter output = StreamWriter(filePath);
            output.WriteLine(proteinGroups.front().GetTabSeparatedHeader());
            for (int i = 0; i < proteinGroups.size(); i++)
            {
                if ((!getParameters()->getSearchParameters()->getWriteDecoys() &&
                     proteinGroups[i]->getIsDecoy()) ||
                    (!getParameters()->getSearchParameters()->getWriteContaminants() &&
                     proteinGroups[i]->getIsContaminant()))
                {
                    continue;
                }
                else if (proteinGroups[i]->getQValue() <= qValueCutoff)
                {
                    output.WriteLine(proteinGroups[i]);
                }
            }
        }
        
        FinishedWritingFile(filePath, nestedIds);
    }
    
    void PostSearchAnalysisTask::WritePeptideQuantificationResultsToTsv(FlashLfqResults *flashLFQResults,
                                                                        const std::string &outputFolder,
                                                                        const std::string &fileName,
                                                                        std::vector<std::string> &nestedIds)
    {
        auto fullSeqPath = outputFolder +"/" + fileName + ".tsv";
        
        flashLFQResults->WriteResults(nullptr, fullSeqPath, nullptr);
        
        FinishedWritingFile(fullSeqPath, nestedIds);
    }
    
    void PostSearchAnalysisTask::WritePeakQuantificationResultsToTsv(FlashLfqResults *flashLFQResults,
                                                                     const std::string &outputFolder,
                                                                     const std::string &fileName,
                                                                     std::vector<std::string> &nestedIds)
    {
        auto peaksPath = outputFolder + "/" +  fileName + ".tsv";
        
        flashLFQResults->WriteResults(peaksPath, nullptr, nullptr);
        
        FinishedWritingFile(peaksPath, nestedIds);
    }
}
