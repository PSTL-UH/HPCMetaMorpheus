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
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include <experimental/filesystem>
#include "Group.h"
#include "Math.h"

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
    
    std::vector<std::pair<std::string, std::vector<PeptideSpectralMatch*>>>& PostSearchAnalysisTask::getPsmsGroupedByFile()
    {
        return privatePsmsGroupedByFile;
    }
    
    
    void PostSearchAnalysisTask::setPsmsGroupedByFile(const std::vector<std::pair<std::string, std::vector<PeptideSpectralMatch*>>> &value)    {
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
                bool found = false;
                for ( auto t : tvec ) {
                    if ( t[0]->getFullFilePath() == psm->getFullFilePath()                      &&
                         t[0]->getScanNumber()   == psm->getScanNumber()                        &&
                         t[0]->getPeptideMonisotopicMass() == psm->getPeptideMonisotopicMass() ) {
                        t.push_back(psm);
                        found = true;
                        break;
                    }
                }
                if ( !found )  {
                    std::vector<PeptideSpectralMatch *> *t = new std::vector<PeptideSpectralMatch *>;
                    t->push_back(psm);
                    tvec.push_back(*t);
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
        auto tempVar = new FdrAnalysisEngine ( tmppsms, massDiffAcceptorNumNotches, getCommonParameters(), svec);
        tempVar->Run();
        
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
        auto tempVar = new ProteinParsimonyEngine (tmppsms, getParameters()->getSearchParameters()->getModPeptidesAreDifferent(),
                                       getCommonParameters(), svec );
        ProteinParsimonyResults *proteinAnalysisResults = static_cast<ProteinParsimonyResults*>(tempVar->Run());
        
        // score protein groups and calculate FDR
        auto tmp = proteinAnalysisResults->getProteinGroups();
        tmppsms = getParameters()->getAllPsms();

        std::cout << "ProteinAnalysis: after ParsimonyEngine: tmp size " << tmp.size() << " tmppsms size " << tmppsms.size() << std::endl;
        for ( auto p : tmp ) {
            std::cout << p->ToString() << std::endl;
        }
        
        auto tempVar2 = new ProteinScoringAndFdrEngine (tmp, tmppsms,
                                            getParameters()->getSearchParameters()->getNoOneHitWonders(),
                                            getParameters()->getSearchParameters()->getModPeptidesAreDifferent(),
                                            true, getCommonParameters(), svec);
        ProteinScoringAndFdrResults *proteinScoringAndFdrResults = static_cast<ProteinScoringAndFdrResults*>(tempVar2->Run());
        
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
                auto tempVar = new LocalizationEngine (tmppsms, myMsDataFile, combinedParams, svec2);
                tempVar->Run();
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
        for ( auto p : getParameters()->getAllPsms() ) {
            if ( p->getFdrInfo()->getQValue() <= 0.01      &&
                 p->getFdrInfo()->getQValueNotch() <= 0.01 &&
                 !p->getIsDecoy()                          &&
                 p->getFullSequence().length() != 0 ) {
                unambiguousPsmsBelowOnePercentFdr.push_back(p);
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
                std::vector<Protein*> proteinsOrderedByAccession(proteinGroup->getProteins().begin(), proteinGroup->getProteins().end());
                std::sort( proteinsOrderedByAccession.begin(), proteinsOrderedByAccession.end(), [&] (Protein *l, Protein *r) {
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
                    if ( p->getGeneNames().size() > 0 ) {
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
                        psmToProteinGroups.emplace(psm, flashLfqProteinGroups);
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
                        psmToProteinGroups.emplace(psm, proteinGroups);
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
            for ( auto p : getParameters()->getAllPsms() ) {
                if ( p->getFdrInfo()->getQValue() <= 0.01 ) {
                    limitedpsms_with_fdr.push_back(p);
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
                
                delete myTreeStructure;
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
        for ( auto p : getParameters()->getAllPsms() ) {
            if ( p->getFdrInfo()->getQValue() <= getCommonParameters()->getQValueOutputFilter()      &&
                 p->getFdrInfo()->getQValueNotch() <= getCommonParameters()->getQValueOutputFilter() ) {
                filteredPsmListForOutput.push_back(p);
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
#ifdef ORIG
        // EDGAR: I don't think this needs a groupby, only first element is used. more like a distinct.
        std::vector<PeptideSpectralMatch*> peptides = getParameters()->getAllPsms().GroupBy([&] (std::any b) {
                b::FullSequence;
            })->Select([&] (std::any b)  {
                    b::FirstOrDefault();
		}).ToList();
#endif
        std::vector<PeptideSpectralMatch*> peptides;
        for ( auto b:  getParameters()->getAllPsms() ) {
            bool found=  false;
            for ( auto p: peptides ) {
                if (b->getFullSequence() == p->getFullSequence() ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                peptides.push_back(b);
            }
        }
        
#ifdef ORIG
        //EDGAR: the same here.
        WritePsmsToTsv(filteredPsmListForOutput.GroupBy([&] (std::any b)  {
                    b::FullSequence;
		})->Select([&] (std::any b) {
                        b::FirstOrDefault();
                    }).ToList(), writtenFile, getParameters()->getSearchParameters()->getModsToWriteSelection());
#endif
        std::vector<PeptideSpectralMatch*> tmppeptides;
        for ( auto b:  filteredPsmListForOutput ) {
            bool found=  false;
            for ( auto p: tmppeptides ) {
                if (b->getFullSequence() == p->getFullSequence() ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                tmppeptides.push_back(b);
            }
        }
        auto tempvar = getParameters()->getSearchParameters()->getModsToWriteSelection();
        WritePsmsToTsv(tmppeptides, writtenFile, &tempvar);        
        std::vector<std::string> svec2 = {getParameters()->getSearchTaskId()};
        FinishedWritingFile(writtenFile, svec2);
        
        // write summary text
#ifdef ORIG
        getParameters()->getSearchTaskResults()->AddNiceText("All target PSMS within 1% FDR: " + getParameters()->getAllPsms().size()([&] (std::any a)    {
                    return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
                }));
#endif
        int count=0;
        for ( auto a :  getParameters()->getAllPsms() ) {
            if ( a->getFdrInfo()->getQValue() <= 0.01 && !a->getIsDecoy() ) {
                count++;
            }
        }
        getParameters()->getSearchTaskResults()->AddNiceText("All target PSMS within 1% FDR: " + std::to_string(count));
                                                             
#ifdef ORIG
        getParameters()->getSearchTaskResults()->AddNiceText("All target peptides within 1% FDR: " + peptides.size()([&] (std::any a){
                    return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
		}));
#endif
        count=0;
        for ( auto a :  peptides ) {
            if ( a->getFdrInfo()->getQValue() <= 0.01 && !a->getIsDecoy() ) {
                count++;
            }
        }
        getParameters()->getSearchTaskResults()->AddNiceText("All target peptides within 1% FDR: " + std::to_string(count));

        if (getParameters()->getSearchParameters()->getDoParsimony())
        {
#ifdef ORIG
            getParameters()->getSearchTaskResults()->AddNiceText("All target protein groups within 1% FDR: " + std::to_string(getProteinGroups().size()([&] (std::any b) {
                            return b::QValue <= 0.01 && !b::IsDecoy;
			})) + "\r\n");
#endif
            count=0;
            for ( auto b :  getProteinGroups() ) {
                if ( b->getQValue() <= 0.01 && !b->getIsDecoy() ) {
                    count++;
                }
            }
            getParameters()->getSearchTaskResults()->AddNiceText("All target protein groups within 1% FDR: " +
                                                                 std::to_string(count) + "\r\n");                       
        }
        
#ifdef ORIG
        setPsmsGroupedByFile(filteredPsmListForOutput.GroupBy([&] (std::any p) {
                    p::FullFilePath;
		}));
#endif
        auto tmparg = new std::vector<std::pair<std::string, std::vector<PeptideSpectralMatch*>>>;
        std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f1 = [&]
            (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
            return l->getFullFilePath() < r->getFullFilePath(); } ;
        std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f2 = [&]
            (PeptideSpectralMatch *l, PeptideSpectralMatch *r) {
            return l->getFullFilePath() != r->getFullFilePath(); } ;
        std::vector<std::vector<PeptideSpectralMatch*>> tmppsms= Group::GroupBy ( filteredPsmListForOutput, f1, f2);
        for ( auto p: tmppsms ) {
            tmparg->push_back(std::make_pair(p[0]->getFullFilePath(), p));
        }
        setPsmsGroupedByFile(*tmparg);

        for (auto file : getPsmsGroupedByFile() )
        {
            // write summary text
            auto psmsForThisFile = file.second; //->ToList();
            std::string fname = file.first;
            std::string strippedFileName = fname.substr(0, fname.find_last_of(".")); 
#ifdef ORIG
            auto peptidesForFile = psmsForThisFile.GroupBy([&] (std::any b)  {
                    b::FullSequence;
                })->Select([&] (std::any b)   {
                        b::FirstOrDefault();
                    }).ToList();
#endif
            std::vector<PeptideSpectralMatch*> peptidesForFile;
            for ( auto b : psmsForThisFile ) {
                bool found = false;
                for ( auto q: peptidesForFile ) {
                    if ( b->getFullSequence() == q->getFullSequence() ) {
                        found = true;
                        break;
                    }
                }
                if ( !found ) {
                    peptidesForFile.push_back(b);
                }
            }
            
            getParameters()->getSearchTaskResults()->AddNiceText("MS2 spectra in " + strippedFileName + ": " +
                                       std::to_string(getParameters()->getNumMs2SpectraPerFile()[strippedFileName][0]));
            getParameters()->getSearchTaskResults()->AddNiceText("Precursors fragmented in " + strippedFileName +
                                  ": " + std::to_string(getParameters()->getNumMs2SpectraPerFile()[strippedFileName][1]));
#ifdef ORIG
            getParameters()->getSearchTaskResults()->AddNiceText("Target PSMs within 1% FDR in " + strippedFileName + ": " +
                                       psmsForThisFile.size()([&] (std::any a) {
                        return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
                    }));
#endif
            int count =0;
            for ( auto a: psmsForThisFile ) {
                if ( a->getFdrInfo()->getQValue() <= 0.01 && !a->getIsDecoy() ) {
                    count++;
                }
            }
            getParameters()->getSearchTaskResults()->AddNiceText("Target PSMs within 1% FDR in " + strippedFileName + ": "
                                                                 + std::to_string(count));
#ifdef ORIG
            getParameters()->getSearchTaskResults()->AddNiceText("Target peptides within 1% FDR in " +
                                                                 strippedFileName + ": " +
                                                                 std::to_string(peptidesForFile.size()([&] (std::any a)  {
                                                                             return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
                                                                         })) + "\r\n");
#endif
            count = 0;
            for ( auto a: peptidesForFile ) {
                if ( a->getFdrInfo()->getQValue() <= 0.01 && !a->getIsDecoy() ) {
                    count++;
                }
            }
            getParameters()->getSearchTaskResults()->AddNiceText("Target peptides within 1% FDR in " +
                                                                 strippedFileName + ": " + std::to_string(count));
            
            // writes all individual spectra file search results to subdirectory
            if (getParameters()->getCurrentRawFileList().size() > 1)
            {
                // create individual files subdirectory
                FileSystem::createDirectory(getParameters()->getIndividualResultsOutputFolder());
                
                // write PSMs
                writtenFile = getParameters()->getIndividualResultsOutputFolder() + "/" + strippedFileName + "_PSMs.psmtsv";
                auto targ = getParameters()->getSearchParameters()->getModsToWriteSelection();
                WritePsmsToTsv(psmsForThisFile, writtenFile, &targ);
                std::vector<std::string> tmpvec = {getParameters()->getSearchTaskId(),
                                                   "Individual Spectra Files", file.second.front()->getFullFilePath()};
                FinishedWritingFile(writtenFile, tmpvec );
                
                // write PSMs for percolator
                writtenFile = getParameters()->getIndividualResultsOutputFolder() + "/" + strippedFileName +
                    "_PSMsFormattedForPercolator.tsv";
                WritePsmsForPercolator(psmsForThisFile, writtenFile, getCommonParameters()->getQValueOutputFilter());
                FinishedWritingFile(writtenFile, tmpvec );
                
                // write best (highest-scoring) PSM per peptide
                writtenFile = getParameters()->getIndividualResultsOutputFolder() + "/" + strippedFileName + "_Peptides.psmtsv";
                auto targ2 = getParameters()->getSearchParameters()->getModsToWriteSelection();
                WritePsmsToTsv(peptidesForFile, writtenFile, &targ2);
                FinishedWritingFile(writtenFile, tmpvec);
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
                
#ifdef ORIG
                //for (auto fullFilePath : getPsmsGroupedByFile().Select([&] (std::any v)    {
                //            v::Key;
                //        }))
#endif
                std::vector<std::string> tmpFilePath;
                for ( auto v: getPsmsGroupedByFile() ) {
                    tmpFilePath.push_back(v.first);
                }

                for (auto fullFilePath : tmpFilePath )
                {
                    std::string strippedFileName = fullFilePath.substr(0,fullFilePath.find_last_of("."));
#ifdef ORIG
                    std::vector<PeptideSpectralMatch*> psmsForThisFile = getPsmsGroupedByFile().Where([&] (std::any p)  {
                            return p->Key == fullFilePath;
                        }).SelectMany([&] (std::any g)   {
                                return g;
                            }).ToList();
#endif
                    std::vector<PeptideSpectralMatch*> psmsForThisFile;
                    for ( auto p: getPsmsGroupedByFile() ) {
                        if ( p.first == fullFilePath ) {
                            for ( auto g : p.second ) {
                                psmsForThisFile.push_back(g);
                            }
                        }   
                    }
                    
#ifdef ORIG
                    auto subsetProteinGroupsForThisFile = getProteinGroups().Select([&] (std::any p)  {
                            p::ConstructSubsetProteinGroup(fullFilePath);
                        }).ToList();
#endif
                    std::vector<EngineLayer::ProteinGroup*> subsetProteinGroupsForThisFile;
                    for ( auto p: getProteinGroups() ) {
                        subsetProteinGroupsForThisFile.push_back(p->ConstructSubsetProteinGroup(fullFilePath));
                    }
                    
                    std::vector<std::string> tmpvec = {getParameters()->getSearchTaskId(),
                                                       "Individual Spectra Files", fullFilePath};
                    auto tempVar = new ProteinScoringAndFdrEngine (subsetProteinGroupsForThisFile, psmsForThisFile,
                                                       getParameters()->getSearchParameters()->getNoOneHitWonders(),
                                                       getParameters()->getSearchParameters()->getModPeptidesAreDifferent(),
                                                       false, getCommonParameters(),
                                                       tmpvec );
                    ProteinScoringAndFdrResults *subsetProteinScoringAndFdrResults = static_cast<ProteinScoringAndFdrResults*>(tempVar->Run());
                    
                    subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults->SortedAndScoredProteinGroups;
                    
#ifdef ORIG
                    getParameters()->getSearchTaskResults()->AddNiceText("Target protein groups within 1 % FDR in " +
                                                                         strippedFileName + ": " +
                                                                         subsetProteinGroupsForThisFile.size()([&] (std::any b)  {
                                                                                 return b::QValue <= 0.01 && !b::IsDecoy;
                                                                             }));
#endif
                    int count=0;
                    for ( auto b: subsetProteinGroupsForThisFile) {
                        if ( b->getQValue() <= 0.01 && !b->getIsDecoy() ) {
                            count++;
                        }
                    }
                    getParameters()->getSearchTaskResults()->AddNiceText("Target protein groups within 1 % FDR in " +
                                                                         strippedFileName + ": " + std::to_string(count));
                    
                    // write individual spectra file protein groups results to tsv
                    if (getParameters()->getCurrentRawFileList().size() > 1)
                    {
                        writtenFile =  getParameters()->getIndividualResultsOutputFolder() +"/" +
                            strippedFileName.substr(strippedFileName.find_last_of("/")) + "_ProteinGroups.tsv";
                        std::vector<std::string> tmpvec2;
                        tmpvec2.push_back(getParameters()->getSearchTaskId());
                        tmpvec2.push_back("Individual Spectra Files");
                        tmpvec2.push_back(fullFilePath);

                        auto targ = getCommonParameters()->getQValueOutputFilter();
                        WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, writtenFile,
                                                tmpvec2, targ);
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
                            strippedFileName.substr(strippedFileName.find_last_of("/")) + ".mzID";
                        auto targ1 = getParameters()->getVariableModifications();
                        auto targ2 = getParameters()->getFixedModifications();
                        std::vector<Protease*> vecarg1 = {getCommonParameters()->getDigestionParams()->getProtease()};
                        MzIdentMLWriter::WriteMzIdentMl(psmsForThisFile, subsetProteinGroupsForThisFile,
                                                        targ1, targ2, vecarg1,
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
                            strippedFileName.substr(strippedFileName.find_last_of("/")) + ".pep.XM";
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
                                                     file.first->FullFilePathWithExtension};
                    WritePeakQuantificationResultsToTsv(getParameters()->getFlashLfqResults(),
                                                        getParameters()->getIndividualResultsOutputFolder(),
                                                        file.first->FilenameWithoutExtension + "_QuantifiedPeaks",
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
                     !b->getIsDecoy()                          &&
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
                    std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>::const_iterator proteinToConfidentBaseSequences_iterator = proteinToConfidentBaseSequences.find(peptide->getProtein()->getNonVariantProtein());
                    if (proteinToConfidentBaseSequences_iterator != proteinToConfidentBaseSequences.end())
                    {
                        myPepList = proteinToConfidentBaseSequences_iterator->second;
                        myPepList.push_back(peptide);
                        proteinToConfidentBaseSequences.emplace(peptide->getProtein()->getNonVariantProtein(),myPepList );
                    }
                    else
                    {
                        //myPepList = proteinToConfidentBaseSequences_iterator->second;
                        std::vector<PeptideWithSetModifications*> tvec = {peptide};
                        proteinToConfidentBaseSequences.emplace(peptide->getProtein()->getNonVariantProtein(), tvec);
                    }
                }
            }
            
            // Add user mod selection behavours to Pruned DB
            for (auto modType : getParameters()->getSearchParameters()->getModsToWriteSelection())
            {
#ifdef ORIG
                //for (Modification *mod : GlobalVariables::getAllModsKnown().Where([&] (std::any b)    {
                //            b::ModificationType->Equals(modType.Key);
                //        }))
#endif
                for (Modification *mod : GlobalVariables::getAllModsKnown() )
                {
                    if ( mod->getModificationType() != modType.first ) {
                        continue;
                    }
                    if (modType.second == 1) // Write if observed and in database
                    {
                        modificationsToWriteIfBoth.insert(mod);
                    }
                    if (modType.second == 2) // Write if in database
                    {
                        modificationsToWriteIfInDatabase.insert(mod);
                    }
                    if (modType.second == 3) // Write if observed
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
                     !b->getIsDecoy()                          &&
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
                    std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>::const_iterator proteinToConfidentModifiedSequences_iterator =
                        proteinToConfidentModifiedSequences.find(peptide->getProtein()->getNonVariantProtein());
                    if (proteinToConfidentModifiedSequences_iterator != proteinToConfidentModifiedSequences.end())
                    {
                        myPepList = proteinToConfidentModifiedSequences_iterator->second;
                        myPepList.push_back(peptide);
                        proteinToConfidentModifiedSequences.emplace(peptide->getProtein()->getNonVariantProtein(), myPepList );
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
#ifdef ORIG
            //for (auto nonVariantProtein : getParameters()->getProteinList().Select([&] (std::any p)  {
            //            p::NonVariantProtein;
            //        }).Distinct())
#endif
            std::vector<Protein*> tmpVarProtein;
            for (auto p: getParameters()->getProteinList() ) {
                bool found = false;
                for ( auto q : tmpVarProtein ) {
                    if ( p->getNonVariantProtein()->Equals(q) ) {
                        found = true;
                        break;
                    }
                }
                if ( !found ) {
                    tmpVarProtein.push_back( p->getNonVariantProtein());
                }
            }

            for ( auto nonVariantProtein : tmpVarProtein )
            {
                if (!nonVariantProtein->getIsDecoy())
                {
                    std::vector<PeptideWithSetModifications*> psms;
                    std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>::const_iterator proteinToConfidentModifiedSequences_iterator =
                        proteinToConfidentModifiedSequences.find(nonVariantProtein);
                    psms = proteinToConfidentModifiedSequences_iterator->second;
                    
                    // sequence variant is null if mod is not on a variant
                    //std::unordered_set<std::tuple<int, Modification*, SequenceVariation*>> modsObservedOnThisProtein;
                    PSATTuple1_set modsObservedOnThisProtein;
                    
                    for (auto psm : (psms.size() != 0 ) ? psms : std::vector<PeptideWithSetModifications*>())
                    {
                        for (auto idxModKV : psm->getAllModsOneIsNterminus())
                        {
                            int proteinIdx = GetOneBasedIndexInProtein(idxModKV.first, psm);
#ifdef ORIG
                            SequenceVariation *relevantVariant = psm->Protein.AppliedSequenceVariations.FirstOrDefault([&]
                                                                                                                       (std::any sv) {
                                    VariantApplication::IsSequenceVariantModification(sv, proteinIdx);
                                });
#endif
                            SequenceVariation *relevantVariant=nullptr;
                            for ( auto sv : psm->getProtein()->getAppliedSequenceVariations() ) {
                                if ( VariantApplication::IsSequenceVariantModification(sv, proteinIdx) ) {
                                    relevantVariant = sv;
                                    break;
                                }
                            }
#ifdef ORIG
                            SequenceVariation *unappliedVariant = relevantVariant == nullptr ? nullptr : psm->Protein.SequenceVariations.FirstOrDefault([&] (std::any sv)  {
                                    return sv::Description != nullptr &&
                                    sv::Description->Equals(relevantVariant->Description);
                                });
#endif
                            SequenceVariation *unappliedVariant = nullptr;
                            if ( relevantVariant != nullptr ) {
                                for ( auto sv : psm->getProtein()->getSequenceVariations() ) {
                                    if ( sv->getDescription() != nullptr &&
                                         sv->getDescription()->Equals(relevantVariant->getDescription()) ) {
                                        unappliedVariant = sv;
                                        break;
                                    }
                                }
                            }
                            
                            modsObservedOnThisProtein.insert(std::make_tuple(
                                                                 VariantApplication::RestoreModificationIndex(psm->getProtein(),
                                                                                                              proteinIdx),
                                                                 idxModKV.second,
                                                                 unappliedVariant));
                        }
                    }
                    
                    //std::unordered_map<std::tuple<SequenceVariation*, int>, std::vector<Modification*>> modsToWrite;
                    PSATTuple2_map modsToWrite;
                    
                    //Add if observed (regardless if in database)
                    for (auto observedMod : modsObservedOnThisProtein)
                    {
                        auto tempMod = std::get<1>(observedMod);
                        
                        if (std::find(modificationsToWriteIfObserved.begin(), modificationsToWriteIfObserved.end(), tempMod) !=
                            modificationsToWriteIfObserved.end())
                        {
                            auto svIdxKey = std::make_tuple(std::get<2>(observedMod), std::get<0>(observedMod));
                            if (modsToWrite.find(svIdxKey) == modsToWrite.end())
                            {
                                std::vector<Modification*> tvec = {std::get<1>(observedMod)};
                                modsToWrite.emplace(svIdxKey, tvec);
                            }
                            else
                            {
                                modsToWrite[svIdxKey].push_back(std::get<1>(observedMod));
                            }
                        }
                    }
                    
                    // Add modification if in database (two cases: always or if observed)
                    for (auto modkv : nonVariantProtein->getOneBasedPossibleLocalizedModifications() )
                    {
                        for (auto mod : modkv.second)
                        {
                            //Add if always In Database or if was observed and in database and not set to not include
                            if (std::find(modificationsToWriteIfInDatabase.begin(), modificationsToWriteIfInDatabase.end(), mod) !=
                                modificationsToWriteIfInDatabase.end()                                                          ||
                                (std::find(modificationsToWriteIfBoth.begin(), modificationsToWriteIfBoth.end(), mod) !=
                                 modificationsToWriteIfBoth.end()                                                               &&
                                 std::find(modsObservedOnThisProtein.begin(), modsObservedOnThisProtein.end(),
                                           std::make_tuple(modkv.first, mod, nullptr) ) != modsObservedOnThisProtein.end()))  {
                                if (modsToWrite.find(std::make_tuple(nullptr, modkv.first)) == modsToWrite.end())
                                {
                                    modsToWrite.emplace(std::make_tuple(nullptr, modkv.first), std::vector<Modification*> {mod});
                                }
                                else
                                {
                                    modsToWrite[std::make_tuple(nullptr, modkv.first)].push_back(mod);
                                }
                            }
                        }
                    }
                    
                    // Add variant modification if in database (two cases: always or if observed)
                    for (SequenceVariation *sv : nonVariantProtein->getSequenceVariations())
                    {
                        for (auto modkv : sv->getOneBasedModifications())
                        {
                            for (auto mod : modkv.second)
                            {
                                //Add if always In Database or if was observed and in database and not set to not include
                                if (std::find(modificationsToWriteIfInDatabase.begin(), modificationsToWriteIfInDatabase.end(), mod) !=
                                    modificationsToWriteIfInDatabase.end()                                                   ||
                                    (std::find(modificationsToWriteIfBoth.begin(), modificationsToWriteIfBoth.end(), mod) !=
                                     modificationsToWriteIfBoth.end()                                                        &&
                                     std::find(modsObservedOnThisProtein.begin(), modsObservedOnThisProtein.end(),
                                               std::make_tuple(modkv.first, mod, sv)) != modsObservedOnThisProtein.end())     )
                                {
                                    if (modsToWrite.find(std::make_tuple(sv, modkv.first)) == modsToWrite.end())
                                    {
                                        modsToWrite.emplace(std::make_tuple(sv, modkv.first), std::vector<Modification*> {mod});
                                    }
                                    else
                                    {
                                        modsToWrite[std::make_tuple(sv, modkv.first)].push_back(mod);
                                    }
                                }
                            }
                        }
                    }
                    
                    if (proteinToConfidentBaseSequences.find(nonVariantProtein->getNonVariantProtein()) !=
                        proteinToConfidentBaseSequences.end())
                    {
                        // adds confidently localized and identified mods
                        nonVariantProtein->getOneBasedPossibleLocalizedModifications().clear();
#ifdef ORIG
                        //for (auto kvp : modsToWrite.Where([&] (std::any kv)  {
                        //            return kv::Key->Item1 == nullptr;
                        //        }))
#endif
                        for (auto kvp : modsToWrite)
                        {
                            if ( std::get<0>(kvp.first) != nullptr )
                            {
                                continue;
                            }
                            
                            nonVariantProtein->getOneBasedPossibleLocalizedModifications()[std::get<1>(kvp.first)] = kvp.second;
                        }
                        for (auto sv : nonVariantProtein->getSequenceVariations())
                        {
                            sv->getOneBasedModifications().clear();
#ifdef ORIG
                            //for (auto kvp : modsToWrite.Where([&] (std::any kv)  {
                            //            return kv::Key->Item1 != nullptr && kv::Key->Item1->Equals(sv);
                            //        }))
#endif
                            for (auto kvp : modsToWrite )
                            {
                                if ( std::get<0>(kvp.first) == nullptr ||
                                     !(std::get<0>(kvp.first))->Equals(sv) ) {
                                    continue;
                                }
                                sv->getOneBasedModifications()[std::get<1>(kvp.first)] =  kvp.second;
                            }
                        }
                    }
                }
                
        
                //writes all proteins
#ifdef ORIG
                //if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)  {
                //    !b::IsContaminant;
                //}))
#endif
                bool cond = false;
                for ( auto b: getParameters()->getDatabaseFilenameList() ) {
                    if ( !b->getIsContaminant() ) {
                        cond = true;
                        break;
                    }
                }
                bool cond2 = false;
                for ( auto b: getParameters()->getDatabaseFilenameList() ) {
                    if ( b->getIsContaminant() ) {
                        cond2 = true;
                        break;
                    }
                }

                //Edgar: calculcate the basename only once to avoid repetitive code.
                std::string outputXMLdbBaseNameWO, outputXMLdbBaseName;
                if ( cond) {
                    std::vector<std::string> fnameVec;
                    for ( auto b: getParameters()->getDatabaseFilenameList() ) {
                        if ( !b->getIsContaminant() ) {
                            std::string fname = b->getFilePath();
                            std::string fname_short = fname.substr(0, fname.find_last_of("."));
                            fnameVec.push_back(fname_short);
                        }
                    }
                    std::string del="-";
                    outputXMLdbBaseNameWO = getParameters()->getOutputFolder() + "/" + StringHelper::join(fnameVec, del);
                }
                if ( cond2 )  {
                    std::vector<std::string> fnameVec;
                    for ( auto b: getParameters()->getDatabaseFilenameList() ) {
                        if ( b->getIsContaminant() ) {
                            std::string fname = b->getFilePath();
                            std::string fname_short = fname.substr(0, fname.find_last_of("."));
                            fnameVec.push_back(fname_short);
                        }
                    }
                    std::string del="-";
                    outputXMLdbBaseName = getParameters()->getOutputFolder() + "/" + StringHelper::join(fnameVec, del);
                }
                
                if ( cond) 
                {
#ifdef ORIG
                    std::string outputXMLdbFullName = FileSystem::combine(getParameters()->getOutputFolder(),
                                                         std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)  {
                                                                     !b::IsContaminant;
                                                                 })->Select([&] (std::any b)  {
                                                                         Path::GetFileNameWithoutExtension(b::FilePath);
                                                                     })) + "pruned.xml");
#endif
                    std::string outputXMLdbFullName =  outputXMLdbBaseNameWO + "pruned.xml";
                    
#ifdef ORG
                    ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(), getParameters()->getProteinList().Select([&] (std::any p)  {
                        p::NonVariantProtein;
                    }).Where([&] (std::any b) {
                            return !b::IsDecoy && !b::IsContaminant;
                        }).ToList(), outputXMLdbFullName);
#endif
                    std::vector<Protein*> proteinList;
                    for ( auto p: getParameters()->getProteinList() ) {
                        auto b = p->getNonVariantProtein();
                        if ( !b->getIsDecoy() && !b->getIsContaminant() ) {
                            proteinList.push_back(p);
                        }
                    }
                    std::unordered_map<std::string,ModDbTuple_set> tmpmap;
                    ProteinDbWriter::WriteXmlDatabase(tmpmap, proteinList, outputXMLdbFullName);

                    std::vector<std::string> tmpvec = {getParameters()->getSearchTaskId()};
                    FinishedWritingFile(outputXMLdbFullName, tmpvec);
                }
            
#ifdef ORIG
                //if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b) {
                //        b::IsContaminant;
                //    }))
#endif
                if ( cond2 )
                {
#ifdef ORIG
                    std::string outputXMLdbFullNameContaminants = FileSystem::combine(getParameters()->getOutputFolder(),
                                                  std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)  {
                                                              b::IsContaminant;
                                                          })->Select([&] (std::any b) {
                                                                  Path::GetFileNameWithoutExtension(b::FilePath);
                                                              })) + "pruned.xml");
#endif
                    std::string outputXMLdbFullNameContaminants = outputXMLdbBaseName +  "pruned.xml";
#ifdef ORIG                    
                    ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(),
                                                      getParameters()->getProteinList().Select([&] (std::any p)  {
                            p::NonVariantProtein;
                        }).Where([&] (std::any b) {
                                return !b::IsDecoy && b::IsContaminant;
                            }).ToList(), outputXMLdbFullNameContaminants);
#endif
                    std::vector<Protein*> proteinList;
                    for ( auto p: getParameters()->getProteinList() ) {
                        auto b = p->getNonVariantProtein();
                        if ( !b->getIsDecoy() && b->getIsContaminant() ) {
                            proteinList.push_back(p);
                        }
                    }
                    std::unordered_map<std::string, ModDbTuple_set> tmpmap;
                    ProteinDbWriter::WriteXmlDatabase(tmpmap, proteinList, outputXMLdbFullNameContaminants);
                    std::vector<std::string> tmpvec = {getParameters()->getSearchTaskId()};
                    FinishedWritingFile(outputXMLdbFullNameContaminants, tmpvec);
                }
            
                //writes only detected proteins
#ifdef ORIG
                //if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)  {
                //        !b::IsContaminant;
                //    }))
#endif
                if ( cond) 
                {
#ifdef ORIG
                    std::string outputXMLdbFullName = FileSystem::combine(getParameters()->getOutputFolder(),
                                           std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)   {
                                         !b::IsContaminant;
                                     })->Select([&] (std::any b) {
                                             Path::GetFileNameWithoutExtension(b::FilePath);
                                         })) + "proteinPruned.xml");
#endif
                    std::string outputXMLdbFullName = outputXMLdbBaseNameWO + "proteinPruned.xml";

#ifdef ORIG
                    ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(),
                                                      proteinToConfidentBaseSequences.Keys->Where([&] (std::any b) {
                            return !b::IsDecoy && !b::IsContaminant;
                        }).ToList(), outputXMLdbFullName);
#endif
                    std::vector<Protein*> proteinList;
                    for ( auto p:  proteinToConfidentBaseSequences ) {
                        auto b = std::get<0>(p);
                        if ( !b->getIsDecoy() && !b->getIsContaminant() ) {
                            proteinList.push_back(b);
                        }
                    }
                    std::unordered_map<std::string, ModDbTuple_set> tmpmap;
                    ProteinDbWriter::WriteXmlDatabase(tmpmap, proteinList, outputXMLdbFullName);
                
                    std::vector<std::string> tmpvec = {getParameters()->getSearchTaskId()};
                    FinishedWritingFile(outputXMLdbFullName, tmpvec);
                }
            
#ifdef ORIG
                //if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)  {
                //        b::IsContaminant;
                //    }))
#endif
                if ( cond2 ) 
                {
#ifdef ORIG
                    std::string outputXMLdbFullNameContaminants = FileSystem::combine(getParameters()->getOutputFolder(),
                                       std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b) {
                                                   b::IsContaminant;
                                               })->Select([&] (std::any b)  {
                                                       Path::GetFileNameWithoutExtension(b::FilePath);
                                                   })) + "proteinPruned.xml");
#endif
                    std::string outputXMLdbFullNameContaminants = outputXMLdbBaseName + "proteinPruned.xml";
#ifdef ORIG            
                    ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(),
                                                      proteinToConfidentBaseSequences.Keys->Where([&] (std::any b)  {
                                                              return !b::IsDecoy && b::IsContaminant;
                                                          }).ToList(), outputXMLdbFullNameContaminants);
#endif
                    std::vector<Protein*> proteinList;
                    for ( auto p:  proteinToConfidentBaseSequences ) {
                        auto b = std::get<0>(p);
                        if ( !b->getIsDecoy() && b->getIsContaminant() ) {
                            proteinList.push_back(b);
                        }
                    }
                    std::unordered_map<std::string, ModDbTuple_set> tmpmap;
                    ProteinDbWriter::WriteXmlDatabase(tmpmap, proteinList, outputXMLdbFullNameContaminants);
                    
                    std::vector<std::string> tmpvec = {getParameters()->getSearchTaskId()};
                    FinishedWritingFile(outputXMLdbFullNameContaminants, tmpvec);
                }
            
            }
        }
    }
        
    int PostSearchAnalysisTask::GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications *peptideWithSetModifications)
    {
        if (oneIsNterminus == 1)
        {
            return peptideWithSetModifications->getOneBasedStartResidueInProtein();
        }
        if (oneIsNterminus == peptideWithSetModifications->getLength() + 2)
        {
            return peptideWithSetModifications->getOneBasedEndResidueInProtein();
        }
        return peptideWithSetModifications->getOneBasedStartResidueInProtein() + oneIsNterminus - 2;
    }
    
    void PostSearchAnalysisTask::WriteTree(BinTreeStructure *myTreeStructure, const std::string &writtenFile)
    {
        std::ofstream output(writtenFile);
        output << "MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tUnimodDiffs\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tprotNtermLocFrac\tpepNtermLocFrac\tpepCtermLocFrac\tprotCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot" << std::endl;

#ifdef ORIG
        //  for (Bin *bin : myTreeStructure->getFinalBins().OrderByDescending([&] (std::any b) {
        //            b->Count;
        //        }))
#endif
        auto tmpBins =  myTreeStructure->getFinalBins();
        std::sort(tmpBins.begin(), tmpBins.end(), [&] (Bin *l, Bin *r) {
                return l->getCount() > r->getCount(); });
        for ( Bin *bin: tmpBins ) 
        {         
#ifdef ORIG
            std::string tmpstring1 = std::string::Join(",", bin->ModsInCommon::OrderByDescending([&] (std::any b){
                        b->Value;
                    }).Where([&] (std::any b)	{
                            return b->Value > bin->CountTarget / 10.0;
                        })->Select([&] (std::any b) {
                                return b::Key + ":" + (static_cast<double>(b->Value) / bin->CountTarget).ToString("F3");
                            }));
#endif
            std::vector<std::tuple<std::string, int>> tmpvec;
            for ( auto b: bin->ModsInCommon ) {
                if ( std::get<1>(b) > (bin->getCountTarget()/10.0) ) {
                    tmpvec.push_back(b);
                }
            }
            std::sort(tmpvec.begin(), tmpvec.end(), [&] (std::tuple<std::string, int> l, std::tuple<std::string, int> r) {
                    return std::get<1>(l) > std::get<1>(r);
                });
            
            std::vector<std::string> tmpstringvec;
            for ( auto b: tmpvec ) {
                std::stringstream ss;
                ss <<  std::get<0>(b) << ":" << std::to_string(static_cast<double>(std::get<1>(b)) / bin->getCountTarget());
                std::string s = ss.str();
                tmpstringvec.push_back(s);
            }
            std::string del = ",";
            std::string tmpstring1 = StringHelper::join(tmpstringvec, del);       
#ifdef ORIG
            std::string tmpstring2 = std::string::Join(",", bin->AAsInCommon::OrderByDescending([&] (std::any b){
                        b->Value;
                    }).Where([&] (std::any b) {
                            return b->Value > bin->CountTarget / 10.0;
                        })->Select([&] (std::any b)	{
                                return b::Key + ":" + (static_cast<double>(b->Value) / bin->CountTarget).ToString("F3");
                            }));
#endif
            std::vector<std::tuple<char, int>> tmpvec2;
            for ( auto b: bin->getAAsInCommon() ) {
                if ( std::get<1>(b) > (bin->getCountTarget()/10.0) ) {
                    tmpvec2.push_back(b);
                }
            }
            std::sort(tmpvec2.begin(), tmpvec2.end(), [&] (std::tuple<char, int> l, std::tuple<char, int> r) {
                    return std::get<1>(l) > std::get<1>(r);
                });
            tmpstringvec.clear();
            for ( auto b: tmpvec2 ) {
                std::stringstream ss;
                ss <<  std::get<0>(b) << ":" << std::to_string(static_cast<double>(std::get<1>(b)) / bin->getCountTarget());
                std::string s = ss.str();
                tmpstringvec.push_back(s);
            }
            std::string tmpstring2 = StringHelper::join(tmpstringvec, del);
#ifdef ORIG
            std::string tmpstring3 = std::string::Join(",", bin->ResidueCount::OrderByDescending([&] (std::any b) {
                        b->Value;
                    })->Select([&] (std::any b){
                            return b::Key + ":" + b->Value;
                        }));
#endif
            tmpvec2.clear();
            for ( auto b: bin->ResidueCount ) {
                tmpvec2.push_back(b);
            }
            std::sort(tmpvec2.begin(), tmpvec2.end(), [&] (std::tuple<char, int> l, std::tuple<char, int> r) {
                    return std::get<1>(l) > std::get<1>(r);
                });
            tmpstringvec.clear();
            for ( auto b: tmpvec2 ) {
                std::stringstream ss;
                ss <<  std::get<0>(b) << ":" << std::get<1>(b);
                std::string s = ss.str();
                tmpstringvec.push_back(s);
            }
            std::string tmpstring3 = StringHelper::join(tmpstringvec, del);
            
            output << std::setprecision(4) << bin->getMassShift() << "\t" <<
                bin->getCount() <<  "\t" <<
                bin->getCountDecoy() << "\t" <<
                bin->getCountTarget()<< "\t" <<
                bin->getLocalizeableTarget() << "\t" <<
                (bin->getCountTarget() - bin->getLocalizeableTarget()) << "\t" <<
                (bin->getCount() == 0 ? NAN : static_cast<double>(bin->getCountDecoy() / bin->getCount())) << "\t" <<
                Math::NormalDistribution ( 0, 1, bin->ComputeZ(0.01)) << "\t" <<
                Math::NormalDistribution (0, 1, bin->ComputeZ(0.255)) << "\t" <<
                (bin->getCountTarget() == 0 ? NAN : static_cast<double>(bin->getLocalizeableTarget() / bin->getCountTarget())) << "\t" <<
                bin->getMine() << "\t" << bin->getUnimodId() << "\t" <<
                bin->getUnimodFormulas() << "\t" <<
                bin->getUnimodDiffs() << "\t" << bin->AA << "\t" <<
                bin->getCombos() << "\t" <<
                tmpstring1 + "\t" <<
                tmpstring2 + "\t" <<
                tmpstring3 + "\t" <<
                (bin->getLocalizeableTarget() == 0 ? NAN : static_cast<double>(bin->getProtNlocCount() / bin->getLocalizeableTarget())) << "\t" <<
                (bin->getLocalizeableTarget() == 0 ? NAN : static_cast<double>(bin->getPepNlocCount() / bin->getLocalizeableTarget())) << "\t" <<
                (bin->getLocalizeableTarget() == 0 ? NAN : static_cast<double>(bin->getPepClocCount() / bin->getLocalizeableTarget())) << "\t" <<
                (bin->getLocalizeableTarget() == 0 ? NAN : static_cast<double>(bin->getProtClocCount() / bin->getLocalizeableTarget())) << "\t" <<
                bin->getFracWithSingle() << "\t" <<
                (static_cast<double>(bin->getOverlapping()) / bin->getCountTarget()) << "\t" <<
                bin->getMedianLength() << "\t" <<  bin->getUniprotID() << std::endl;
        }
        
        output.close();
    }
    
    void PostSearchAnalysisTask::WritePsmsForPercolator(std::vector<PeptideSpectralMatch*> &psmList,
                                                        const std::string &writtenFileForPercolator,
                                                        double qValueCutoff)
    {
        std::ofstream output(writtenFileForPercolator);
        output << "SpecId\tLabel\tScanNr\tF1\tF2\tPeptide\tProteins" << std::endl;
        output << "DefaultDirection\t-\t-\t1\t1\t\t" << std::endl;
        for (int i = 0; i < (int)psmList.size(); i++)
        {
            auto psm = psmList[i];
            
            if (psm->getFdrInfo()->getQValue() > qValueCutoff || psm->getFdrInfo()->getQValueNotch() > qValueCutoff)
            {
                continue;
            }
            
            output << std::to_string(i);
            output << '\t' << std::to_string(psm->getIsDecoy() ? -1 : 1);
            output << '\t' << psm->getScanNumber();
            
            // Features
            std::string del = "\t";
            auto vec = psm->getFeatures();
            std::vector<std::string> stringvec;
            for ( auto v: vec ) {
                stringvec.push_back(std::to_string(v));
            }
            output << StringHelper::toString('\t') << StringHelper::join(stringvec, del );
            
            // HACKY: Ignores all ambiguity
            auto pwsm = std::get<1>(psm->getBestMatchingPeptides().front());
            
            output << '\t' << pwsm->getPreviousAminoAcid() << "." << pwsm->getFullSequence() << "." << pwsm->getNextAminoAcid();
            output << '\t' << pwsm->getProtein()->getAccession();
            output << std::endl;
        }
        output.close();
    }
    
    
    void PostSearchAnalysisTask::WriteProteinGroupsToTsv(std::vector<EngineLayer::ProteinGroup*> &proteinGroups,
                                                         const std::string &filePath,
                                                         std::vector<std::string> &nestedIds,
                                                         double qValueCutoff)
    {
        if ( proteinGroups.size() > 0 )
        {
            std::ofstream output(filePath);
            output << proteinGroups.front()->GetTabSeparatedHeader() << std::endl;
            for (int i = 0; i < (int)proteinGroups.size(); i++)
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
                    output << proteinGroups[i] << std::endl;
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
        
        std::string s1="", s2="";
        flashLFQResults->WriteResults(s1, fullSeqPath, s2);
        
        FinishedWritingFile(fullSeqPath, nestedIds);
    }
    
    void PostSearchAnalysisTask::WritePeakQuantificationResultsToTsv(FlashLfqResults *flashLFQResults,
                                                                     const std::string &outputFolder,
                                                                     const std::string &fileName,
                                                                     std::vector<std::string> &nestedIds)
    {
        auto peaksPath = outputFolder + "/" +  fileName + ".tsv";
        
        std::string s1="", s2="";
        flashLFQResults->WriteResults(peaksPath, s1, s2);
        
        FinishedWritingFile(peaksPath, nestedIds);
    }
}
