#include "PeptideSpectralMatch.h"
#include "IScan.h"
#include "bankersrounding.h"
#include "stringhelper.h"
#include "GlobalVariables.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <tuple>
#include <experimental/filesystem>
#include <optional>
#include <numeric>
#include <iomanip>
#include <sstream>

#include "Group.h"

using namespace Chemistry;
using namespace EngineLayer::FdrAnalysis;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

    PeptideSpectralMatch::PeptideSpectralMatch(PeptideWithSetModifications *peptide, int notch, double score, int scanIndex, IScan *scan, DigestionParams *digestionParams, std::vector<MatchedFragmentIon*> &matchedFragmentIons) 
    {
        _bestMatchingPeptides = std::vector<std::tuple<int, PeptideWithSetModifications*>>();
        privateScanIndex = scanIndex;
        privateFullFilePath = scan->getFullFilePath();
        privateScanNumber = scan->getOneBasedScanNumber();
        privatePrecursorScanNumber = scan->getOneBasedPrecursorScanNumber();
        privateScanRetentionTime = scan->getRetentionTime();
        privateScanExperimentalPeaks = scan->getNumPeaks();
        privateTotalIonCurrent = scan->getTotalIonCurrent();
        privateScanPrecursorCharge = scan->getPrecursorCharge();
        privateScanPrecursorMonoisotopicPeakMz = scan->getPrecursorMonoisotopicPeakMz();
        privateScanPrecursorMass = scan->getPrecursorMass();
        setAllScores(std::vector<double>());
        setPeptidesToMatchingFragments(std::unordered_map<PeptideWithSetModifications*, std::vector<MatchedFragmentIon*>>());
        this->digestionParams = new DigestionParams(*digestionParams);
        
        AddOrReplace(peptide, score, notch, true, matchedFragmentIons);
    }
    
    std::vector<std::tuple<int, PeptideWithSetModifications *>> PeptideSpectralMatch::getBestMatchingPeptides() const
    {
        return _bestMatchingPeptides;
    }

    ChemicalFormula *PeptideSpectralMatch::getModsChemicalFormula() const
    {
        return privateModsChemicalFormula;
    }
    
    void PeptideSpectralMatch::setModsChemicalFormula(ChemicalFormula *value)
    {
        privateModsChemicalFormula = value;
    }
    
    std::string PeptideSpectralMatch::getFullSequence() const
    {
        return privateFullSequence;
    }
    
    void PeptideSpectralMatch::setFullSequence(const std::string &value)
    {
        privateFullSequence = value;
    }
    
    std::optional<int> PeptideSpectralMatch::getNotch() const
    {
        return privateNotch;
    }
    
    void PeptideSpectralMatch::setNotch(const std::optional<int> &value)
    {
        privateNotch = value;
    }
    
    std::string PeptideSpectralMatch::getBaseSequence() const
    {
        return privateBaseSequence;
    }
    
    void PeptideSpectralMatch::setBaseSequence(const std::string &value)
    {
        privateBaseSequence = value;
    }
    
    std::optional<int> PeptideSpectralMatch::getPeptideLength() const
    {
        return privatePeptideLength;
    }
    
    void PeptideSpectralMatch::setPeptideLength(const std::optional<int> &value)
    {
        privatePeptideLength = value;
    }
    
    std::optional<int> PeptideSpectralMatch::getOneBasedStartResidueInProtein() const
    {
        return privateOneBasedStartResidueInProtein;
    }
    
    void PeptideSpectralMatch::setOneBasedStartResidueInProtein(const std::optional<int> &value)
    {
        privateOneBasedStartResidueInProtein = value;
    }
    
    std::optional<int> PeptideSpectralMatch::getOneBasedEndResidueInProtein() const
    {
        return privateOneBasedEndResidueInProtein;
    }
    
    void PeptideSpectralMatch::setOneBasedEndResidueInProtein(const std::optional<int> &value)
    {
        privateOneBasedEndResidueInProtein = value;
    }
    
    std::optional<double> PeptideSpectralMatch::getPeptideMonisotopicMass() const
    {
        return privatePeptideMonisotopicMass;
    }
    
    void PeptideSpectralMatch::setPeptideMonisotopicMass(const std::optional<double> &value)
    {
        privatePeptideMonisotopicMass = value;
    }
    
    std::optional<int> PeptideSpectralMatch::getProteinLength() const
    {
        return privateProteinLength;
    }
    
    void PeptideSpectralMatch::setProteinLength(const std::optional<int> &value)
    {
        privateProteinLength = value;
    }
    
    std::string PeptideSpectralMatch::getProteinAccession() const
    {
        return privateProteinAccession;
    }
    
    void PeptideSpectralMatch::setProteinAccession(const std::string &value)
    {
        privateProteinAccession = value;
    }
    
    std::string PeptideSpectralMatch::getOrganism() const
    {
        return privateOrganism;
    }
    
    void PeptideSpectralMatch::setOrganism(const std::string &value)
    {
        privateOrganism = value;
    }
    
    std::vector<MatchedFragmentIon*> PeptideSpectralMatch::getMatchedFragmentIons() const
    {
        return privateMatchedFragmentIons;
    }
    
    void PeptideSpectralMatch::setMatchedFragmentIons(const std::vector<MatchedFragmentIon*> &value)
    {
        privateMatchedFragmentIons = value;
    }
    
    std::unordered_map<std::string, int> PeptideSpectralMatch::getModsIdentified() const
    {
        return privateModsIdentified;
    }
    
    void PeptideSpectralMatch::setModsIdentified(const std::unordered_map<std::string, int> &value)
    {
        privateModsIdentified = value;
    }
    
    std::vector<double> PeptideSpectralMatch::getLocalizedScores() const
    {
        return privateLocalizedScores;
    }
    
    void PeptideSpectralMatch::setLocalizedScores(const std::vector<double> &value)
    {
        privateLocalizedScores = value;
    }
    
    int PeptideSpectralMatch::getScanNumber() const
    {
        return privateScanNumber;
    }
    
    std::optional<int> PeptideSpectralMatch::getPrecursorScanNumber() const
    {
        return privatePrecursorScanNumber;
    }
    
    double PeptideSpectralMatch::getScanRetentionTime() const
    {
        return privateScanRetentionTime;
    }
    
    int PeptideSpectralMatch::getScanExperimentalPeaks() const
    {
        return privateScanExperimentalPeaks;
    }
    
    double PeptideSpectralMatch::getTotalIonCurrent() const
    {
        return privateTotalIonCurrent;
    }
    
    int PeptideSpectralMatch::getScanPrecursorCharge() const
    {
        return privateScanPrecursorCharge;
    }
    
    double PeptideSpectralMatch::getScanPrecursorMonoisotopicPeakMz() const
    {
        return privateScanPrecursorMonoisotopicPeakMz;
    }
    
    double PeptideSpectralMatch::getScanPrecursorMass() const
    {
        return privateScanPrecursorMass;
    }
    
    std::string PeptideSpectralMatch::getFullFilePath() const
    {
        return privateFullFilePath;
    }
    
    int PeptideSpectralMatch::getScanIndex() const
    {
        return privateScanIndex;
    }
    
    int PeptideSpectralMatch::getNumDifferentMatchingPeptides() const
    {
        return _bestMatchingPeptides.size();
    }
    
    FdrInfo *PeptideSpectralMatch::getFdrInfo() const
    {
        return privateFdrInfo;
    }
    
    void PeptideSpectralMatch::setFdrInfo(FdrInfo *value)
    {
        privateFdrInfo = value;
    }
    
    double PeptideSpectralMatch::getScore() const
    {
        return privateScore;
    }
    
    void PeptideSpectralMatch::setScore(double value)
    {
        privateScore = value;
    }
    
    double PeptideSpectralMatch::getDeltaScore() const
    {
        return privateDeltaScore;
    }
    
    void PeptideSpectralMatch::setDeltaScore(double value)
    {
        privateDeltaScore = value;
    }
    
    double PeptideSpectralMatch::getRunnerUpScore() const
    {
        return privateRunnerUpScore;
    }
    
    void PeptideSpectralMatch::setRunnerUpScore(double value)
    {
        privateRunnerUpScore = value;
    }
    
    bool PeptideSpectralMatch::getIsDecoy() const
    {
        return privateIsDecoy;
    }
    
    void PeptideSpectralMatch::setIsDecoy(bool value)
    {
        privateIsDecoy = value;
    }
    
    bool PeptideSpectralMatch::getIsContaminant() const
    {
        return privateIsContaminant;
    }
    
    void PeptideSpectralMatch::setIsContaminant(bool value)
    {
        privateIsContaminant = value;
    }
    
    std::vector<double> PeptideSpectralMatch::getAllScores() const
    {
        return privateAllScores;
    }
    
    void PeptideSpectralMatch::setAllScores(const std::vector<double> &value)
    {
        privateAllScores = value;
    }
    
    std::unordered_map<PeptideWithSetModifications*, std::vector<MatchedFragmentIon*>> PeptideSpectralMatch::getPeptidesToMatchingFragments() const
    {
        return privatePeptidesToMatchingFragments;
    }
    
    void PeptideSpectralMatch::setPeptidesToMatchingFragments(const std::unordered_map<PeptideWithSetModifications*, std::vector<MatchedFragmentIon*>> &value)
    {
        privatePeptidesToMatchingFragments = value;
    }
    
#ifdef ORIG
    // Edgar: not entirely sure what this is.
    private *IEnumerable < PeptideSpectralMatch::(int Notch, PeptideWithSetModifications *Peptide)
    {
        get;
        
	get
	{
            return _bestMatchingPeptides.OrderBy([&] (std::any p) {
                    p::Item2->FullSequence;
                }).ThenBy([&] (std::any p)  {
                        p::Item2->Protein.Accession;
                    }).ThenBy([&] (std::any p)  {
                            p::Item2->OneBasedStartResidueInProtein;
			});
	}
    }
#endif
    
    std::vector<double> PeptideSpectralMatch::getFeatures() const
    {
        return std::vector<double> {BankersRounding::round(getScore()), getScore() - BankersRounding::round(getScore())};
    }
    
    std::string PeptideSpectralMatch::GetTabSeparatedHeader()
    {
        // return std::string::Join("\t", DataDictionary(nullptr, nullptr).Keys);
        //std::unordered_map<std::string, std::string> dict = DataDictionary(nullptr, nullptr);
        std::vector<std::tuple<std::string, std::string>> dict = DataDictionary(nullptr, nullptr);
        std::stringstream ss;
        for ( auto it = dict.begin(); it != dict.end(); it++ ) {
            ss << std::get<0>(*it) << "\t";
         }
        return ss.str();
    }
    
    void PeptideSpectralMatch::AddOrReplace(PeptideWithSetModifications *pwsm, double newScore, int notch,
                                            bool reportAllAmbiguity, std::vector<MatchedFragmentIon*> &matchedFragmentIons)
    {
        if (newScore - getScore() > ToleranceForScoreDifferentiation) //if new score beat the old score, overwrite it
        {
            _bestMatchingPeptides.clear();
            _bestMatchingPeptides.push_back(std::make_tuple(notch, pwsm));
            
            if (getScore() - getRunnerUpScore() > ToleranceForScoreDifferentiation)
            {
                setRunnerUpScore(getScore());
            }
            
            setScore(newScore);
            
            privatePeptidesToMatchingFragments.clear();
            privatePeptidesToMatchingFragments.emplace(pwsm, matchedFragmentIons);
        }
        else if (newScore - getScore() > -ToleranceForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
        {
            _bestMatchingPeptides.push_back(std::make_tuple(notch, pwsm));
            
            if (privatePeptidesToMatchingFragments.find(pwsm) == privatePeptidesToMatchingFragments.end())
            {
                privatePeptidesToMatchingFragments.emplace(pwsm, matchedFragmentIons);
            }
        }
        else if (getScore() - getRunnerUpScore() > ToleranceForScoreDifferentiation)
        {
            setRunnerUpScore(newScore);
        }
    }
    
    std::string PeptideSpectralMatch::ToString()
    {
        std::unordered_map<std::string, int> mymap;
        return ToString(&mymap);
    }
    
    std::string PeptideSpectralMatch::ToString(std::unordered_map<std::string, int> *ModstoWritePruned)
    {
        // return std::string::Join("\t", DataDictionary(this, ModstoWritePruned).Values);
        //std::unordered_map<std::string, std::string> dict = DataDictionary(this, ModstoWritePruned);
        std::vector<std::tuple<std::string, std::string>> dict = DataDictionary(this, ModstoWritePruned);
        std::stringstream ss;
        for ( auto it = dict.begin(); it != dict.end(); it++ ) {
            ss << std::get<1>(*it)  <<  "\t";
        }
        return ss.str();
    }
    
    std::vector<std::tuple<std::string, std::string>> PeptideSpectralMatch::DataDictionary(PeptideSpectralMatch *psm,
                                                             std::unordered_map<std::string, int> *ModsToWritePruned)
    {
        std::vector<std::tuple<std::string, std::string>> s;
        AddBasicMatchData(s, psm);
        AddPeptideSequenceData(s, psm, ModsToWritePruned);
        AddMatchedIonsData(s, psm);
        AddMatchScoreData(s, psm);
        return s;
    }
    
    void PeptideSpectralMatch::CalculateDeltaScore(double scoreCutoff)
    {
        setDeltaScore(getScore() - std::max(getRunnerUpScore(), scoreCutoff));
    }
    
    void PeptideSpectralMatch::SetFdrValues(double cumulativeTarget, double cumulativeDecoy, double qValue, double cumulativeTargetNotch,
                                            double cumulativeDecoyNotch, double qValueNotch, double maximumLikelihood, double eValue,
                                            double eScore, bool calculateEValue)
    {
        FdrInfo* tempVar= new FdrInfo();
        setFdrInfo(tempVar);
        getFdrInfo()->setCumulativeTarget(cumulativeTarget);
        getFdrInfo()->setCumulativeDecoy(cumulativeDecoy);
        getFdrInfo()->setQValue(qValue);
        getFdrInfo()->setCumulativeTargetNotch(cumulativeTargetNotch);
        getFdrInfo()->setCumulativeDecoyNotch(cumulativeDecoyNotch);
        getFdrInfo()->setQValueNotch(qValueNotch);
        getFdrInfo()->setMaximumLikelihood(maximumLikelihood);
        getFdrInfo()->setEScore(eScore);
        getFdrInfo()->setEValue(eValue);
        getFdrInfo()->setCalculateEValue(calculateEValue);
    }
    
    void PeptideSpectralMatch::ResolveAllAmbiguities()
    {
#ifdef ORIG
        // EDGAR Note: Pwsm is the name of the second element in the tuple
        // in the C# version. The C++ version does not use 'names' for
        // the elements of the tuple. Its still the second element,
        // so we access it using std::get<1>
        setIsDecoy(_bestMatchingPeptides.Any([&] (std::any p)  {
                    p::Pwsm::Protein::IsDecoy;
                }));
#endif
        bool found = false;
        for ( auto p: _bestMatchingPeptides ) {
            if ( std::get<1>(p)->getProtein()->getIsDecoy()) {
                found = true;
                break;
            }
        }
        setIsDecoy(found);
        
#ifdef ORIG
        setIsContaminant(_bestMatchingPeptides.Any([&] (std::any p) {
                    p::Pwsm::Protein::IsContaminant;
                }));
#endif
        found = false;
        for ( auto p: _bestMatchingPeptides ) {
            if ( std::get<1>(p)->getProtein()->getIsContaminant()) {
                found = true;
                break;
            }
        }
        setIsContaminant(found);
        
#ifdef ORIG
        setFullSequence(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                        b::Pwsm::FullSequence;
                    }))->ResolvedValue);
#endif
        std::vector<std::string> vs;
        for ( auto b: _bestMatchingPeptides ) {
            vs.push_back( std::get<1>(b)->getFullSequence() );
        }
        std::tuple<std::string, std::string> res = Resolve(vs);
        setFullSequence(std::get<1>(res));
        
        
#ifdef ORIG
        setBaseSequence(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                        b::Pwsm::BaseSequence;
                    }))->ResolvedValue);
#endif
        vs.clear();
        for ( auto b: _bestMatchingPeptides ) {
            vs.push_back( std::get<1>(b)->getBaseSequence() );
        }
        std::tuple<std::string, std::string> res2 = Resolve(vs);
        setBaseSequence(std::get<1>(res2));
        
        
#ifdef ORIG
        setPeptideLength(Resolve(_bestMatchingPeptides.Select([&] (std::any b){
                        b::Pwsm->Length;
                    }))->ResolvedValue);
#endif
        std::vector<int> vi;
        for ( auto b: _bestMatchingPeptides ) {
            vi.push_back( std::get<1>(b)->getLength() );
        }
        std::tuple<std::string, std::optional<int>> res3 = Resolve(vi);
        setPeptideLength(std::get<1>(res3));
        
#ifdef ORIG            
        setOneBasedStartResidueInProtein(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                        b::Pwsm::OneBasedStartResidueInProtein;
                    }))->ResolvedValue);
#endif
        vi.clear();
        for ( auto b: _bestMatchingPeptides ) {
            vi.push_back( std::get<1>(b)->getOneBasedStartResidueInProtein() );
        }
        std::tuple<std::string, std::optional<int>> res4 = Resolve(vi);
        setOneBasedStartResidueInProtein(std::get<1>(res4));
        
#ifdef ORIG
        setOneBasedEndResidueInProtein(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                        b::Pwsm::OneBasedEndResidueInProtein;
                    }))->ResolvedValue);
#endif
        vi.clear();
        for ( auto b: _bestMatchingPeptides ) {
            vi.push_back( std::get<1>(b)->getOneBasedEndResidueInProtein() );
        }
        std::tuple<std::string, std::optional<int>> res5 = Resolve(vi);
        setOneBasedEndResidueInProtein(std::get<1>(res5));
        
#ifdef ORIG
        setProteinLength(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                        b::Pwsm::Protein->Length;
                    }))->ResolvedValue);
#endif
        vi.clear();
        for ( auto b: _bestMatchingPeptides ) {
            vi.push_back( std::get<1>(b)->getProtein()->getLength() );
        }
        std::tuple<std::string, std::optional<int>> res6 = Resolve(vi);
        setProteinLength(std::get<1>(res6));
        
#ifdef ORIG
        setPeptideMonisotopicMass(Resolve(_bestMatchingPeptides.Select([&] (std::any b)  {
                        b::Pwsm::MonoisotopicMass;
                    }))->ResolvedValue);
#endif
        std::vector <double> vo;
        for ( auto b: _bestMatchingPeptides ) {
            vo.push_back( std::get<1>(b)->getMonoisotopicMass() );
        }
        std::tuple<std::string, std::optional<double>> res7 = Resolve(vo);
        setPeptideMonisotopicMass(std::get<1>(res7));
        
#ifdef ORIG
        setProteinAccession(Resolve(_bestMatchingPeptides.Select([&] (std::any b)	{
                        b::Pwsm::Protein::Accession;
                    }))->ResolvedValue);
#endif
        vs.clear();
        for ( auto b: _bestMatchingPeptides ) {
            vs.push_back( std::get<1>(b)->getProtein()->getAccession() );
        }
        std::tuple<std::string, std::string> res8 = Resolve(vs);
        setProteinAccession(std::get<1>(res8));
        
#ifdef ORIG
        setOrganism(Resolve(_bestMatchingPeptides.Select([&] (std::any b)	{
                        b::Pwsm::Protein::Organism;
                    }))->ResolvedValue);
#endif
        vs.clear();
        for ( auto b: _bestMatchingPeptides ) {
            vs.push_back( std::get<1>(b)->getProtein()->getOrganism() );
        }
        std::tuple<std::string, std::string> res9 = Resolve(vs);
        setOrganism(std::get<1>(res9));
        
#ifdef ORIG
        setModsIdentified(Resolve(_bestMatchingPeptides.Select([&] (std::any b)	{
                        b::Pwsm::AllModsOneIsNterminus;
                    }))->ResolvedValue);
#endif
        std::vector<std::unordered_map<int, Modification*>> vmap;
        for ( auto b: _bestMatchingPeptides ) {
            vmap.push_back( std::get<1>(b)->getAllModsOneIsNterminus() );
        }
        std::tuple<std::string, std::unordered_map<std::string, int>> res10 = Resolve(vmap);
        setModsIdentified(std::get<1>(res10));
        
#ifdef ORIG
        setModsChemicalFormula(Resolve(_bestMatchingPeptides.Select([&] (std::any b){
                        b::Pwsm::AllModsOneIsNterminus->Select([&] (std::any c)  {
                                (c->Value);
                            });
                    }))->ResolvedValue);
#endif
        std::vector<std::vector<Modification *>> vmod(_bestMatchingPeptides.size() );
        int i=0;
        for ( auto b: _bestMatchingPeptides ) {
            for ( auto c: std::get<1>(b)->getAllModsOneIsNterminus() ) {
                vmod[i].push_back(std::get<1>(c));
            }
            i++;
        }
        std::tuple<std::string, ChemicalFormula *> res11 = Resolve (vmod);
        setModsChemicalFormula(std::get<1>(res11));
        
#ifdef ORIG
        setNotch(Resolve(_bestMatchingPeptides.Select([&] (std::any b)  {
                        b::Notch;
                    }))->ResolvedValue);
#endif
        vi.clear();
        for ( auto b: _bestMatchingPeptides ) {
            vi.push_back(std::get<0>(b));
        }
        std::tuple<std::string, std::optional<int>> res12 = Resolve (vi);
        setNotch(std::get<1>(res12));
        
        
        // if the PSM matches a target and a decoy and they are the SAME SEQUENCE, remove the decoy
        if (getIsDecoy())
        {
            bool removedPeptides = false;
#ifdef ORIG
            auto hits = _bestMatchingPeptides.GroupBy([&] (std::any p)  {
                    p::Pwsm::FullSequence;
                });
#endif     
            std::function<bool(std::tuple<int, PeptideWithSetModifications*>,std::tuple<int, PeptideWithSetModifications*>)> f1 = [&](std::tuple<int, PeptideWithSetModifications*>l, std::tuple<int, PeptideWithSetModifications*>r) {
                return  std::get<1>(l)->getFullSequence() < std::get<1>(r)->getFullSequence(); };
            std::function<bool(std::tuple<int, PeptideWithSetModifications*>,std::tuple<int, PeptideWithSetModifications*>)> f2 = [&](std::tuple<int, PeptideWithSetModifications*>l, std::tuple<int, PeptideWithSetModifications*>r) {
                return  std::get<1>(l)->getFullSequence() != std::get<1>(r)->getFullSequence(); } ;
            std::vector<std::vector<std::tuple<int, PeptideWithSetModifications*>>> hits = Group::GroupBy(_bestMatchingPeptides, f1, f2);
            
            for (auto hit : hits)
            {
# ifdef ORIG
                if (hit->Any([&] (std::any p) {
                            p::Pwsm::Protein::IsDecoy;
                        })                   &&
                    hit->Any([&] (std::any p) {
                            !p::Pwsm::Protein::IsDecoy;
                        })) {
                    _bestMatchingPeptides.RemoveAll([&] (std::any p)   {
                            return p::Pwsm->FullSequence == hit->Key && p::Pwsm::Protein::IsDecoy;
                        });
                    {
                        removedPeptides = true;
                    }
                }
#endif
                bool someisdecoy=false, someisnotdecoy=false; 
                for ( auto p: hit ) {
                    if ( std::get<1>(p)->getProtein()->getIsDecoy() ) {
                        someisdecoy = true;
                        break;
                    }
                }
                for ( auto p: hit ) {
                    if ( !std::get<1>(p)->getProtein()->getIsDecoy() ) {
                        someisnotdecoy = true;
                        break;
                    }
                }
                
                if ( someisdecoy && someisnotdecoy ) {
                    for ( auto p = _bestMatchingPeptides.begin(); p != _bestMatchingPeptides.end(); p++ ) {
                        if ( std::get<1>(*p)->getFullSequence() == std::get<1>(hit[0])->getFullSequence() &&
                             std::get<1>(*p)->getProtein()->getIsDecoy()             ) {
                            _bestMatchingPeptides.erase(p);
                            removedPeptides = true;
                        }
                    }
                }                                
            }
            
            if (removedPeptides)
            {
                ResolveAllAmbiguities();
            }
        }
        // TODO: technically, different peptide options for this PSM can have different matched ions
        // we can write a Resolve method for this if we want...
        //MatchedFragmentIons = PeptidesToMatchingFragments.First().Value;
        privateMatchedFragmentIons = std::get<1>(*privatePeptidesToMatchingFragments.begin());
    }
    
    void PeptideSpectralMatch::TrimProteinMatches(std::vector<Protein*> parsimoniousProteins)
    {
        if (privateIsDecoy)
        {
#ifdef ORIG
            if (_bestMatchingPeptides::Any([&] (std::any p) {
                        return std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(),
                                         p::Pwsm::Protein) != parsimoniousProtei ns.end() &&
                            p::Pwsm::Protein::IsDecoy;
                    }))
            {
                _bestMatchingPeptides::RemoveAll([&] (std::any p) {
                        std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(), p::Item2->Protein) ==
                            parsimoniousProteins.end ();
                    });
            }
#endif
            bool found =false;
            for ( auto p = _bestMatchingPeptides.begin(); p !=_bestMatchingPeptides.end(); p++ ) {
                if ( std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(),
                               std::get<1>(*p)->getProtein() ) != parsimoniousProteins.end() &&
                     std::get<1>(*p)->getProtein()->getIsDecoy() ) {
                    found = true;
                    break;
                }
            }
            if ( found ) {
                for ( auto p = _bestMatchingPeptides.begin(); p !=_bestMatchingPeptides.end(); p++ ) {
                    if ( std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(),
                                   std::get<1>(*p)->getProtein() ) == parsimoniousProteins.end()) {
                        _bestMatchingPeptides.erase (p);
                    }
                }
            }
            
            
            // else do nothing
        }
        else
        {
#ifdef ORIG
            _bestMatchingPeptides::RemoveAll([&] (std::any p){
                    std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(), p::Item2->Protein) == parsimoniousProteins.end();
                });
#endif
            for ( auto p = _bestMatchingPeptides.begin(); p !=_bestMatchingPeptides.end(); p++ ) {
                if ( std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(),
                               std::get<1>(*p)->getProtein() ) == parsimoniousProteins.end()) {
                    _bestMatchingPeptides.erase (p);
                }
            }
        }
        
        ResolveAllAmbiguities();
    }
    
    void PeptideSpectralMatch::AddProteinMatch(std::tuple<int, PeptideWithSetModifications*>  peptideWithNotch)
    {
        _bestMatchingPeptides.push_back(peptideWithNotch);
        ResolveAllAmbiguities();
    }
    
    void PeptideSpectralMatch::AddBasicMatchData(std::vector<std::tuple<std::string, std::string>> &s, PeptideSpectralMatch *psm)
    {
        if ( psm == nullptr ) {
            s.push_back(std::make_tuple("File Name", " "));
            s.push_back(std::make_tuple("Scan Number"  , " "));
            s.push_back(std::make_tuple("Scan Retention Time" , " "));
            s.push_back(std::make_tuple("Num Experimental Peaks"  , " "));
            s.push_back(std::make_tuple("Total Ion Current"  , " "));
            s.push_back(std::make_tuple("Precursor Scan Number", " "));
            s.push_back(std::make_tuple("Precursor Charge"  , " "));
            s.push_back(std::make_tuple("Precursor MZ"  , " "));
            s.push_back(std::make_tuple("Precursor Mass" , " "));
            s.push_back(std::make_tuple("Score" , " "));
            s.push_back(std::make_tuple("Delta Score" , " "));
            s.push_back(std::make_tuple("Notch",  " "));
            s.push_back(std::make_tuple("Different Peak Matches" , " "));

        }
        else {
            s.push_back(std::make_tuple("File Name", std::experimental::filesystem::path(psm->getFullFilePath()).stem() ));            
            s.push_back(std::make_tuple("Scan Number" ,std::to_string(psm->getScanNumber() ) ));
            s.push_back(std::make_tuple("Scan Retention Time" ,std::to_string(psm->getScanRetentionTime())));
            s.push_back(std::make_tuple("Num Experimental Peaks" ,std::to_string(psm->getScanExperimentalPeaks())));
            s.push_back(std::make_tuple("Total Ion Current" ,std::to_string(psm->getTotalIonCurrent())));
            s.push_back(std::make_tuple("Precursor Scan Number" ,psm->getPrecursorScanNumber().has_value() ? std::to_string(psm->getPrecursorScanNumber().value()) : "unknown"));
            s.push_back(std::make_tuple("Precursor Charge" ,std::to_string(psm->getScanPrecursorCharge())));
            s.push_back(std::make_tuple("Precursor MZ" ,std::to_string(psm->getScanPrecursorMonoisotopicPeakMz())));
            s.push_back(std::make_tuple("Precursor Mass" ,std::to_string(psm->getScanPrecursorMass())));
            s.push_back(std::make_tuple("Score" ,std::to_string(psm->getScore())));
            s.push_back(std::make_tuple("Delta Score" ,std::to_string(psm->getDeltaScore())));
            std::vector<int> v;
            for ( auto p: psm->_bestMatchingPeptides ) {
                v.push_back (std::get<0>(p));
            }
            std::tuple<std::string, std::optional<int>> res = Resolve ( v);
            s.push_back(std::make_tuple("Notch", std::get<0>(res)));            
            s.push_back(std::make_tuple("Different Peak Matches", std::to_string(psm->getNumDifferentMatchingPeptides())));
        }
                
    }
    
    void PeptideSpectralMatch::AddPeptideSequenceData(std::vector<std::tuple<std::string, std::string>>& s,
                                                      PeptideSpectralMatch *psm,
                                                      std::unordered_map<std::string, int> *ModsToWritePruned)
    {
        bool pepWithModsIsNull = psm == nullptr || psm->_bestMatchingPeptides.empty() == true;
        
        std::vector<PeptideWithSetModifications*> pepsWithMods;
        if ( !pepWithModsIsNull ) {
            for ( auto p: psm->_bestMatchingPeptides ) {
                pepsWithMods.push_back(std::get<1>(p));
            }
        }

        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Base Sequence", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                vs.push_back(b->getBaseSequence());
            }
            std::tuple<std::string, std::string>res = Resolve(vs);
            s.push_back(std::make_tuple("Base Sequence", std::get<0>(res) ));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Full Sequence", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                vs.push_back(b->getFullSequence());
            }
            std::tuple<std::string, std::string>res = Resolve(vs);
            s.push_back(std::make_tuple("Full Sequence", std::get<0>(res) ));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Essential Sequence", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                vs.push_back(b->EssentialSequence(ModsToWritePruned));
            }
            std::tuple<std::string, std::string>res = Resolve(vs);
            s.push_back(std::make_tuple("Essential Sequence", std::get<0>(res) ));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Mods", " "));
        }
        else {
            std::vector<std::unordered_map<int, Modification*>> vc;
            for ( auto b : pepsWithMods ) {
                vc.push_back( b->getAllModsOneIsNterminus() );
            }
            std::tuple<std::string, std::unordered_map<std::string, int>>res = Resolve(vc);
            s.push_back(std::make_tuple("Mods", std::get<0>(res) ));
        }
        
#ifdef ORIG
        s["Mods Chemical Formulas"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any p)
        {
            p::AllModsOneIsNterminus->Select([&] (std::any v)
        {
            v->Value;
        });
        }))->ResolvedString;
#endif
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Mods Chemical Formulas", " "));
        }
        else {
            std::vector<std::vector<Modification *>> vm(pepsWithMods.size() );
            int i=0;
            for ( auto b : pepsWithMods ) {
                for ( auto v : b->getAllModsOneIsNterminus() ) {
                    vm[i].push_back(std::get<1>(v));
                }
                i++;
            }
            std::tuple<std::string, ChemicalFormula *>res = Resolve(vm);
            s.push_back(std::make_tuple("Mods Chemical Formulas", std::get<0>(res) ));
        }
        
#ifdef ORIG
        s["Mods Combined Chemical Formula"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
        {
            b::AllModsOneIsNterminus->Select([&] (std::any c)
        {
            (dynamic_cast<Modification*>(c->Value));
        });
        }))->ResolvedString;
#endif
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Mods Combined Chemical Formula", " "));
        }
        else {
            std::vector<std::vector<Modification *>> vm(pepsWithMods.size() );
            int i=0;
            for ( auto b : pepsWithMods ) {
                for ( auto v : b->getAllModsOneIsNterminus() ) {
                    vm[i].push_back(std::get<1>(v));
                }
                i++;
            }
            std::tuple<std::string, ChemicalFormula *>res = Resolve(vm);
            s.push_back(std::make_tuple("Mods Combined Chemical Formula", std::get<0>(res) ));
        }
        
#ifdef ORIG
        s["Num Variable Mods"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
        {
            b::NumVariableMods;
        }))->Item1;
#endif
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Num Variable Mods", " "));
        }
        else {
            std::vector<int> vi;
            for ( auto b : pepsWithMods ) {
                vi.push_back(b->getNumVariableMods());
            }
            std::tuple<std::string, std::optional<int>>res = Resolve(vi);
            s.push_back(std::make_tuple("Num Variable Mods", std::get<0>(res) ));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Missed Cleavages", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                vs.push_back(std::to_string(b->getMissedCleavages()));
            }
            std::tuple<std::string, std::string>res = Resolve(vs);
            s.push_back(std::make_tuple("Missed Cleavages", std::get<0>(res) ));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Peptide Monoisotopic Mass", " "));
        }
        else {
            std::vector<double> vd;
            for ( auto b : pepsWithMods ) {
                vd.push_back(b->getMonoisotopicMass());
            }
            std::tuple<std::string, std::optional<double>>res = Resolve(vd);
            s.push_back(std::make_tuple("Peptide Monoisotopic Mass", std::get<0>(res) ));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Mass Diff (Da)", " "));
        }
        else {
            std::vector<double> vd;
            for ( auto b : pepsWithMods ) {
                vd.push_back((psm->getScanPrecursorMass() - b->getMonoisotopicMass() ));
            }
            std::tuple<std::string, std::optional<double>>res = Resolve(vd);
            s.push_back(std::make_tuple("Mass Diff (Da)", std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Mass Diff (ppm)", " "));
        }
        else {
            std::vector<double> vd;
            for ( auto b : pepsWithMods ) {
                vd.push_back((psm->getScanPrecursorMass() - b->getMonoisotopicMass())/ b->getMonoisotopicMass() * 1e6 );
            }
            std::tuple<std::string, std::optional<double>>res = Resolve(vd);
            s.push_back(std::make_tuple("Mass Diff (ppm)", std::get<0>(res) ));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Protein Accession", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                vs.push_back(b->getProtein()->getAccession());
            }
            std::tuple<std::string, std::string>res = Resolve(vs, psm->getFullSequence());
            s.push_back(std::make_tuple("Protein Accession", std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Protein Name", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                vs.push_back(b->getProtein()->getFullName());
            }
            std::tuple<std::string, std::string>res = Resolve(vs, psm->getFullSequence());
            s.push_back(std::make_tuple("Protein Name", std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Gene Name", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                std::string s;
                std::vector<std::tuple<std::string, std::string>> genes = b->getProtein()->getGeneNames();
                for ( auto d = genes.begin(); d != genes.end(); d++  ) {
                    if ( d == (genes.end() - 1 )){ 
                        s += std::get<0>(*d) + ":" + std::get<1>(*d) ;
                    }
                    else {
                        s += std::get<0>(*d) + ":" + std::get<1>(*d) + ", ";
                    }
                }
                vs.push_back(s);
            }
            std::tuple<std::string, std::string>res = Resolve(vs, psm->getFullSequence() );
            s.push_back(std::make_tuple("Gene Name", std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Intersecting Sequence Variations", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                std::string s;
                std::vector<SequenceVariation*> variations = b->getProtein()->getAppliedSequenceVariations();
                for ( auto av = variations.begin(); av != variations.end(); av++ ) {
                    if ( IntersectsWithVariation(b, *av, false ) ) {
                        s += SequenceVariantString(b, *av ) + ", ";
                    }
                    vs.push_back(s);
                }
            }
            std::tuple<std::string, std::string>res = Resolve(vs );
            s.push_back(std::make_tuple("Intersecting Sequence Variations", std::get<0>(res) ));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Identified Sequence Variations", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                std::string s;
                std::vector<SequenceVariation*> variations = b->getProtein()->getAppliedSequenceVariations();
                for ( auto av = variations.begin(); av != variations.end(); av++ ) {
                    if ( IntersectsWithVariation(b, *av, true ) ) {
                        s += SequenceVariantString(b, *av ) + ", ";
                    }
                    vs.push_back(s);
                }
            }
            std::tuple<std::string, std::string>res = Resolve(vs );
            s.push_back(std::make_tuple("Identified Sequence Variations", std::get<0>(res) ));
        }
        
#ifdef ORIG
        s["Splice Sites"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
        {
            std::string::Join(", ", b::Protein::SpliceSites::Where([&] (std::any d)
        {
            Includes(b, d);
        })->Select([&] (std::any d)
        {
            StringHelper::formatSimple("{0}-{1}", d::OneBasedBeginPosition.ToString(), d::OneBasedEndPosition.ToString());
        }));
        }))->ResolvedString;
#endif
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Splice Sites", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                std::string s;
                std::vector<SpliceSite*> sites = b->getProtein()->getSpliceSites();
                for ( auto d = sites.begin(); d != sites.end(); d++ ) {
                    if ( Includes (b, *d) ) {
                        s+= std::to_string((*d)->getOneBasedBeginPosition()) + "-" + std::to_string((*d)->getOneBasedEndPosition());
                    }
                    vs.push_back(s);                        
                }
            }
            std::tuple<std::string, std::string>res = Resolve(vs );
            s.push_back(std::make_tuple("Splice Sites", std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Organism Name", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                vs.push_back(b->getProtein()->getOrganism());
            }
            std::tuple<std::string, std::string>res = Resolve(vs);
            s.push_back(std::make_tuple("Organism Name", std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Contaminant", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                std::string yes = "Y";
                std::string no = "N";
                if ( b->getProtein()->getIsContaminant() ) {
                    vs.push_back(yes);
                }
                else {
                    vs.push_back(no);
                }
            }
            std::tuple<std::string, std::string>res = Resolve(vs);
            s.push_back(std::make_tuple("Contaminant", std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Decoy", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                std::string yes = "Y";
                std::string no = "N";
                if ( b->getProtein()->getIsDecoy() ) {
                    vs.push_back(yes);
                }
                else {
                    vs.push_back(no);
                }
            }
            std::tuple<std::string, std::string>res = Resolve(vs );
            s.push_back(std::make_tuple("Decoy", std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Peptide Description", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                vs.push_back(b->getPeptideDescription() );
            }
            std::tuple<std::string, std::string>res = Resolve(vs);
            s.push_back(std::make_tuple("Peptide Description", std::get<0>(res) ));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Start and End Residues In Protein", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                std::string s;
                s = "[" + std::to_string(b->getOneBasedStartResidueInProtein()) + " to " +
                    std::to_string(b->getOneBasedEndResidueInProtein()) + "]";
                vs.push_back(s);
            }
            std::tuple<std::string, std::string>res = Resolve(vs, psm->getFullSequence() );
            s.push_back(std::make_tuple("Start and End Residues In Protein",  std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Previous Amino Acid", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                std::stringstream ss;
                ss << b->getPreviousAminoAcid();
                vs.push_back(ss.str() );
            }
            std::tuple<std::string, std::string>res = Resolve(vs);
            s.push_back(std::make_tuple("Previous Amino Acid", std::get<0>(res)));
        }
        
        if ( pepWithModsIsNull ) {
            s.push_back(std::make_tuple("Next Amino Acid", " "));
        }
        else {
            std::vector<std::string> vs;
            for ( auto b : pepsWithMods ) {
                std::stringstream ss;
                ss << b->getNextAminoAcid();
                vs.push_back(ss.str() );
            }
            std::tuple<std::string, std::string>res = Resolve(vs );
            s.push_back(std::make_tuple("Next Amino Acid", std::get<0>(res)));
        }
        
        std::string allScores = " ";
        std::string theoreticalsSearched = " ";
        if (!pepWithModsIsNull && psm->getFdrInfo() != nullptr && psm->getFdrInfo()->getCalculateEValue())
        {
            std::vector<double> scores = psm->getAllScores();
            for ( auto p = scores.begin(); p != scores.end(); p++ ) {
                allScores += std::to_string(*p) + ";";
            }
            theoreticalsSearched = std::to_string(psm->getAllScores().size());
        }
        
        s.push_back(std::make_tuple("All Scores", allScores));
        s.push_back(std::make_tuple("Theoreticals Searched", theoreticalsSearched));
        s.push_back(std::make_tuple("Decoy/Contaminant/Target", pepWithModsIsNull ? " " : psm->getIsDecoy() ? "D" : psm->getIsContaminant() ? "C" : "T" ));
    }
    
    
    bool PeptideSpectralMatch::Includes(PeptideWithSetModifications *pep, SpliceSite *site)
    {
        return pep->getOneBasedStartResidueInProtein() <= site->getOneBasedBeginPosition() &&
            pep->getOneBasedEndResidueInProtein() >= site->getOneBasedEndPosition();
    }
    
    bool PeptideSpectralMatch::IntersectsWithVariation(PeptideWithSetModifications *pep, SequenceVariation *appliedVariation,
                                                       bool checkUnique)
    {
        // does it intersect? 
        int intersectOneBasedStart = std::max(pep->getOneBasedStartResidueInProtein(), appliedVariation->getOneBasedBeginPosition());
        int intersectOneBasedEnd = std::min(pep->getOneBasedEndResidueInProtein(), appliedVariation->getOneBasedEndPosition());
        if (intersectOneBasedEnd < intersectOneBasedStart)
        {
            return false;
        }
        else if (!checkUnique)
        {
            return true;
        }
        else
        {
            // if the original sequence is too short or long, the intersect of the peptide and variant is unique
            int intersectSize = intersectOneBasedEnd - intersectOneBasedStart + 1;
            int variantZeroBasedStart = intersectOneBasedStart - appliedVariation->getOneBasedBeginPosition();
            bool origSeqIsShort = (int)appliedVariation->getOriginalSequence().size() - variantZeroBasedStart < intersectSize;
            bool origSeqIsLong = (int)appliedVariation->getOriginalSequence().size() > intersectSize &&
                pep->getOneBasedEndResidueInProtein() > intersectOneBasedEnd;
            if (origSeqIsShort || origSeqIsLong)
            {
                return true;
            }
            
            // is the variant sequence intersecting the peptide different than the original sequence?
            std::string originalAtIntersect = appliedVariation->getOriginalSequence().substr(intersectOneBasedStart -
                                                           appliedVariation->getOneBasedBeginPosition(), intersectSize);
                std::string variantAtIntersect = appliedVariation->getVariantSequence().substr(intersectOneBasedStart -
                                                              appliedVariation->getOneBasedBeginPosition(), intersectSize);
                return originalAtIntersect != variantAtIntersect;
        }
    }
    
    std::string PeptideSpectralMatch::SequenceVariantString(PeptideWithSetModifications *p, SequenceVariation *applied)
    {
#ifdef ORIG
        auto modsOnVariantOneIsNTerm = p->AllModsOneIsNterminus.Where([&] (std::any kv)  {
                return kv->Key == 1 && applied->OneBasedBeginPosition == 1 ||
                applied->OneBasedBeginPosition <= kv::Key - 2 + p->OneBasedStartResidueInProtein &&
                        kv::Key - 2 + p->OneBasedStartResidueInProtein <= applied->OneBasedEndPosition;
            }).ToDictionary([&] (std::any kv)   {
                    return kv::Key - applied->OneBasedBeginPosition + 1;
                }, [&] (std::any kv)    {
                    kv->Value;
                });
#endif
        std::unordered_map<int, Modification*> modsOnVariantOneIsNTerm;
        std::unordered_map<int, Modification*> tmp = p->getAllModsOneIsNterminus();
        for ( auto kv = tmp.begin(); kv != tmp.end(); kv++ ) {
            if ( (kv->first == 1 && applied->getOneBasedBeginPosition() == 1)  ||
                 (applied->getOneBasedBeginPosition() <= kv->first - 2 + p->getOneBasedStartResidueInProtein() &&
                  (kv->first -2 + p->getOneBasedStartResidueInProtein()) <= applied ->getOneBasedEndPosition() )) {
                modsOnVariantOneIsNTerm.emplace((kv->first - applied->getOneBasedBeginPosition() + 1),kv->second);
            }
        }
        
        PeptideWithSetModifications *variantWithAnyMods = new PeptideWithSetModifications(p->getProtein(),
                                                                                          p->getDigestionParams(),
                                                                                          applied->getOneBasedBeginPosition(),
                                                                                          applied->getOneBasedEndPosition(),
                                                                                          p->getCleavageSpecificityForFdrCategory(),
                                                                                          p->getPeptideDescription(),
                                                                                          p->getMissedCleavages(),
                                                                                          modsOnVariantOneIsNTerm,
                                                                                          p->NumFixedMods);
        
        delete variantWithAnyMods;
        return StringHelper::formatSimple("{0}{1}{2}", applied->getOriginalSequence(), applied->getOneBasedBeginPosition(),
                                          variantWithAnyMods->getFullSequence());
    }
    
    void PeptideSpectralMatch::AddMatchedIonsData(std::vector<std::tuple<std::string, std::string>> &s, PeptideSpectralMatch *psm)
    {
        bool nullPsm = (psm == nullptr);
        
        StringBuilder *seriesStringBuilder = new StringBuilder();
        StringBuilder *mzStringBuilder = new StringBuilder();
        StringBuilder *fragmentDaErrorStringBuilder = new StringBuilder();
        StringBuilder *fragmentPpmErrorStringBuilder = new StringBuilder();
        StringBuilder *fragmentIntensityStringBuilder = new StringBuilder();
        std::vector<StringBuilder*> stringBuilders = {seriesStringBuilder, mzStringBuilder, fragmentDaErrorStringBuilder,
                                                      fragmentPpmErrorStringBuilder, fragmentIntensityStringBuilder};
        
        if (!nullPsm)
        {
            std::vector<MatchedFragmentIon*> matchedIons = psm->getMatchedFragmentIons();
            if (matchedIons.empty())
            {
                matchedIons = psm->getPeptidesToMatchingFragments().begin()->second;
            }
            // using ", " instead of "," improves human readability
            const std::string delimiter = ", ";
            
#ifdef ORIG
            auto matchedIonsGroupedByProductType = matchedIons->GroupBy([&] (std::any i) {
                    i::NeutralTheoreticalProduct::ProductType;
                }).OrderBy([&] (std::any i)  {
                        i::Key;
                    }).ToList();
#endif
            std::function<bool(MatchedFragmentIon *,MatchedFragmentIon *)> f1 = [&](MatchedFragmentIon *l, MatchedFragmentIon *r) {
                return  l->NeutralTheoreticalProduct->productType < r->NeutralTheoreticalProduct->productType; };
            std::function<bool(MatchedFragmentIon *,MatchedFragmentIon *)> f2 = [&](MatchedFragmentIon *l, MatchedFragmentIon *r) {
                return  l->NeutralTheoreticalProduct->productType != r->NeutralTheoreticalProduct->productType; } ;
            std::vector<std::vector<MatchedFragmentIon*>> matchedIonsGroupedByProductType = Group::GroupBy(matchedIons, f1, f2);
            
            
            for (auto productType : matchedIonsGroupedByProductType)
            {
#ifdef ORIG
                auto products = productType.OrderBy([&] (std::any p) {
                        p::NeutralTheoreticalProduct::TerminusFragment::FragmentNumber;
                    }).ToList();
#endif
                std::vector<MatchedFragmentIon *> products (productType.begin(), productType.end());
                std::sort(products.begin(), products.end(), [&] (MatchedFragmentIon *l, MatchedFragmentIon *r) {
                        return l->NeutralTheoreticalProduct->TerminusFragment->FragmentNumber <
                            r->NeutralTheoreticalProduct->TerminusFragment->FragmentNumber; } );
                
                std::for_each(stringBuilders.begin(), stringBuilders.end(), [&] (StringBuilder *p) {
                        p->append("[");
                    });
                
                for (int i = 0; i < (int)products.size(); i++)
                {
                    MatchedFragmentIon *ion = products[i];
                    std::string ionLabel;
                    
                    double massError = Chemistry::ClassExtensions::ToMass(ion->Mz, ion->Charge) -
                        ion->NeutralTheoreticalProduct->NeutralMass;
                    double ppmMassError = massError / ion->NeutralTheoreticalProduct->NeutralMass * 1e6;
                    
                    if (ion->NeutralTheoreticalProduct->NeutralLoss == 0)
                    {
                        // no neutral loss
                        auto pty = ion->NeutralTheoreticalProduct->productType;
                        ionLabel = ProductTypeToString(pty) + 
                            std::to_string(ion->NeutralTheoreticalProduct->TerminusFragment->FragmentNumber) + "+" +
                            std::to_string(ion->Charge);
                    }
                    else
                    {
                        // ion label with neutral loss
                        auto pty = ion->NeutralTheoreticalProduct->productType;                        
                        ionLabel = "(" + ProductTypeToString(pty) + 
                            std::to_string(ion->NeutralTheoreticalProduct->TerminusFragment->FragmentNumber) + "-" +
                            std::to_string(ion->NeutralTheoreticalProduct->NeutralLoss) + ")" + "+" +
                            std::to_string(ion->Charge);
                    }
                    
                    // append ion label
                    seriesStringBuilder->append(ionLabel);
                    
                    // append experimental m/z
                    mzStringBuilder->append(ionLabel + ":" + std::to_string(ion->Mz));
                    
                    // append absolute mass error
                    fragmentDaErrorStringBuilder->append(ionLabel + ":" + std::to_string(massError));
                    
                    // append ppm mass error                        
                    fragmentPpmErrorStringBuilder->append(ionLabel + ":" + std::to_string(ppmMassError));
                    
                    // append fragment ion intensity
                    fragmentIntensityStringBuilder->append(ionLabel + ":" + std::to_string(ion->Intensity));
                    
                    // append delimiter ", "
                    if (i < (int)products.size() - 1)
                    {
                        std::for_each(stringBuilders.begin(), stringBuilders.end(), [&] (StringBuilder* p)  {
                                p->append(delimiter);
                            });
                    }
                }
                
                // append product type delimiter                           
                std::for_each (stringBuilders.begin(), stringBuilders.end(), [&] (StringBuilder *p)  {
                        p->append("];");
                    });
            }
        }
        
        std::string sep=";";
        s.push_back(std::make_tuple("Matched Ion Series", nullPsm ? " " : StringHelper::trimEnd(seriesStringBuilder->toString(), sep) ));
        
        //C# TO C++ CONVERTER TODO TASK: The following lambda expression could not be converted:
        //stringBuilders.ForEach(p => TangibleLambdaToken86)s["Matched Ion Series"]);
        s.push_back(std::make_tuple("Matched Ion Mass-To-Charge Ratios",  nullPsm ? " " : StringHelper::trimEnd(mzStringBuilder->toString(), sep )));
        s.push_back(std::make_tuple("Matched Ion Mass Diff (Da)", nullPsm ? " " : StringHelper::trimEnd(fragmentDaErrorStringBuilder->toString(), sep)));
        s.push_back(std::make_tuple("Matched Ion Mass Diff (Ppm)", nullPsm ? " " : StringHelper::trimEnd(fragmentPpmErrorStringBuilder->toString(), sep )));
        s.push_back(std::make_tuple("Matched Ion Intensities", nullPsm ? " " : StringHelper::trimEnd(fragmentIntensityStringBuilder->toString(), sep)));
        
        // number of matched ions
        s.push_back(std::make_tuple("Matched Ion Counts", nullPsm ? " " : std::to_string(psm->getMatchedFragmentIons().size() )));
    }
    
    void PeptideSpectralMatch::AddMatchScoreData(std::vector<std::tuple<std::string, std::string>> &s, PeptideSpectralMatch *peptide)
    {
        std::string localizedScores = " ";
        std::string improvementPossible = " ";
	if (peptide != nullptr && !peptide->getLocalizedScores().empty() ) {
            //	localizedScores = GlobalVariables.CheckLengthOfOutput(("[" + string.Join(",",
            //                  peptide.LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]"));
            std::vector<std::string>vs;
            for ( auto b : peptide->getLocalizedScores() ) {
                std::stringstream ss;
                ss << std::setprecision(6) << b;
                vs.push_back(ss.str());
            }
            std::string del=",";
            localizedScores = "[" + StringHelper::join(vs, del) + "]";
            //	improvementPossible = (peptide.LocalizedScores.Max() - peptide.Score).ToString("F3", CultureInfo.InvariantCulture);
            //auto max = std::max_element(peptide->getLocalizedScores().begin(), peptide->getLocalizedScores().end()); 
            auto max = peptide->getLocalizedScores()[0];
            for ( auto p : peptide->getLocalizedScores() ) {
                if ( p > max ) {
                    max = p;
                }
            }
            std::stringstream ss2;
            ss2 << (max - peptide->getScore());
            improvementPossible = ss2.str();
	}
	s.push_back(std::make_tuple("Localized Scores", localizedScores));
        s.push_back(std::make_tuple("Improvement Possible", improvementPossible));
        std::string cumulativeTarget = " ";
        std::string cumulativeDecoy = " ";
        std::string qValue = " ";
        std::string cumulativeTargetNotch = " ";
        std::string cumulativeDecoyNotch = " ";
        std::string qValueNotch = " ";
        std::string eValue = " ";
        std::string eScore = " ";
	if (peptide != nullptr && peptide->getFdrInfo() != nullptr) {
            std::stringstream ss;
            ss << peptide->getFdrInfo()->getCumulativeTarget();
            cumulativeTarget = ss.str();

            ss.str("");
            ss << peptide->getFdrInfo()->getCumulativeDecoy();
            cumulativeDecoy = ss.str();

            ss.str("");
            ss << peptide->getFdrInfo()->getQValue();
            qValue = ss.str();

            ss.str("");
            ss << peptide->getFdrInfo()->getCumulativeTargetNotch();
            cumulativeTargetNotch = ss.str();

            ss.str("");
            ss << peptide->getFdrInfo()->getCumulativeDecoyNotch();
            cumulativeDecoyNotch = ss.str();

            ss.str("");
            ss << peptide->getFdrInfo()->getQValueNotch();
            qValueNotch = ss.str();
            

            if (peptide->getFdrInfo()->getCalculateEValue() ) {
                ss.str("");
                ss << peptide->getFdrInfo()->getEValue();
                eValue = ss.str();

                ss.str("");
                ss << peptide->getFdrInfo()->getEScore();
                eScore = ss.str();
            }
	}
	s.push_back(std::make_tuple("Cumulative Target", cumulativeTarget));
	s.push_back(std::make_tuple("Cumulative Decoy", cumulativeDecoy));
	s.push_back(std::make_tuple("QValue", qValue));
	s.push_back(std::make_tuple("Cumulative Target Notch", cumulativeTargetNotch));
	s.push_back(std::make_tuple("Cumulative Decoy Notch", cumulativeDecoyNotch));
	s.push_back(std::make_tuple("QValue Notch", qValueNotch));
	s.push_back(std::make_tuple("eValue", eValue));
	s.push_back(std::make_tuple("eScore", eScore));
		        
    }
    
    //static(string ResolvedString, ChemicalFormula ResolvedValue) Resolve(IEnumerable<IEnumerable<Modification>> enumerable);
    std::tuple<std::string, ChemicalFormula*> PeptideSpectralMatch::Resolve(std::vector<std::vector<Modification *>> enumerable)
    {
        ChemicalFormula *firstChemFormula = new ChemicalFormula();
        bool equals = true;

        if ( enumerable.empty() ) {
            return (std::make_tuple("unknown", nullptr));
        }
        
        for ( auto firstMods : enumerable.front() ) {
            if (firstMods == nullptr || firstMods->getChemicalFormula() == nullptr) {
                delete firstChemFormula;
                return (std::make_tuple("unknown", nullptr));
            }
            firstChemFormula->Add(firstMods->getChemicalFormula() );
        }

        std::vector<ChemicalFormula*> formulas;
        for (auto anEnum : enumerable ) {
            ChemicalFormula *fhere = new ChemicalFormula();
            for ( auto mod : anEnum ) {
                if (mod == nullptr || mod->getChemicalFormula() == nullptr)	{
                    delete firstChemFormula;
                    delete fhere;
                    for ( auto p: formulas ) {
                        if ( p != nullptr ) {
                            delete p;
                        }
                    }
                    return (std::make_tuple("unknown", nullptr));
                }
                fhere->Add(mod->getChemicalFormula());
            }
            if (!firstChemFormula->Equals(fhere)) {
                equals = false;
            }
            formulas.push_back(fhere);
        }


        if (!equals) {
            std::vector<std::string> vs;
            for ( auto f: formulas ) {
                vs.push_back ( f->getFormula() );
            }
            std::string returnString = GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, '|'));
            delete firstChemFormula;
            for ( auto p: formulas ) {
                if ( p != nullptr ) {
                    delete p;
                }
            }
            return (std::make_tuple(returnString, nullptr));
        }
        else {
            for ( auto p: formulas ) {
                if ( p != nullptr ) {
                    delete p;
                }
            }
            return (std::make_tuple(firstChemFormula->getFormula(), firstChemFormula));
        }
    }  

    std::tuple<std::string, std::unordered_map<std::string, int>> PeptideSpectralMatch::Resolve(std::vector<std::unordered_map<int, Modification*>> enumerable)
    {
        //var list = enumerable.ToList();
        //Dictionary<string, int> firstDict = list[0].Values.OrderBy(b => b.IdWithMotif).GroupBy(b =>
        //                                    b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());
        std::vector<Modification*> fvec;
        for ( auto b = enumerable.front().begin(); b != enumerable.front().end(); b++ ) {
            fvec.push_back(b->second);
        }
        std::sort ( fvec.begin(), fvec.end(), [&] (Modification *r, Modification *l) {
                return r->getIdWithMotif() < l->getIdWithMotif(); 
            });
        std::vector<std::vector<Modification *>> firstVec;
        int current = 0;
        for ( auto o = fvec.begin(); o != fvec.end(); o++ ) {
            if ( o == fvec.begin() ) {
                std::vector<Modification *> *v = new std::vector<Modification *>;
                firstVec.push_back(*v);
            }
            else {
                auto p = o - 1;
                if ( (*o)->getIdWithMotif() != (*p)->getIdWithMotif() ) {
                    std::vector<Modification *> *v = new std::vector<Modification *>;
                    firstVec.push_back(*v);
                    current++;
                }                
            }
            firstVec[current].push_back(*o);
        }

        std::unordered_map<std::string, int> firstDict; 
        for ( std::vector<Modification *> b : firstVec ) {
            firstDict.emplace ( b[0]->getIdWithMotif(), b.size() );
        }

        bool equals = true;
        //foreach (var dict in list) {
        //	Dictionary<string, int> okTest = dict.Values.OrderBy(b => b.IdWithMotif).GroupBy(b =>
        //                                   b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());
        //	if (!firstDict.SequenceEqual(okTest)) {
        //		equals = false;
        //		break;
        //	}
        //}
        for ( auto dict : enumerable ) {
            if ( dict == *enumerable.begin() ){
                continue;
            }
            std::vector<Modification*> mvec;
            for ( auto b = dict.begin(); b != dict.end(); b++ ) {
                mvec.push_back(b->second);
            }
            std::sort ( mvec.begin(), mvec.end(), [&] (Modification *r, Modification *l) {
                    return r->getIdWithMotif() < l->getIdWithMotif(); 
                });
            std::vector<std::vector<Modification *>> tempVec;
            int current = 0;
            for ( auto o = mvec.begin(); o != mvec.end(); o++ ) {
                if ( o == mvec.begin() ) {
                    std::vector<Modification *> *v = new std::vector<Modification *>;
                    tempVec.push_back(*v);
                }
                else {
                    auto p = o - 1;
                    if ( (*o)->getIdWithMotif() != (*p)->getIdWithMotif() ) {
                        std::vector<Modification *> *v = new std::vector<Modification *>;
                        tempVec.push_back(*v);
                        current++;
                    }                
                }
                tempVec[current].push_back(*o);
            }
            
            std::unordered_map<std::string, int> tempDict; 
            for ( std::vector<Modification *> b : tempVec ) {
                tempDict.emplace ( b[0]->getIdWithMotif(), b.size() );
            }
            // Compare tempDict to firstDict
            if ( tempDict.size() != firstDict.size()) {
                equals = false;
                break;
            }
            std::unordered_map<std::string, int>::iterator j = firstDict.begin();
            for ( auto i : tempDict ) {
                if ( i.first != j->first ||i.second != j->second ) {
                    equals = false;
                    break;
                }
                j++;
            }
        }
        //if (!equals) {
        //	var returnString = string.Join("|", list.Select(b => string.Join(" ", b.Values.Select(c => c.IdWithMotif).OrderBy(c => c))));
        //	returnString = GlobalVariables.CheckLengthOfOutput(returnString);
        //	return (returnString, nullptr);
        //}
        //else {
        //	return (string.Join(" ", list[0].Values.Select(c => c.IdWithMotif).OrderBy(c => c)), firstDict);
        //}
        if ( !equals ) {
            std::vector<std::string>vs;
            for ( auto i: enumerable) {
                for ( auto j = i.begin(); j != i.end(); j++ ) {
                    vs.push_back(j->second->getIdWithMotif());
                }
            }
            std::sort( vs.begin(), vs.end() );
            std::string returnString = GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, ' '));
            std::unordered_map<std::string, int> emptyMap;
            return ( std::make_tuple ( returnString, emptyMap ));            
        }
        else {
            std::vector<std::string>vs;
            for ( auto i = fvec.begin(); i != fvec.end(); i++ ) {
                vs.push_back( (*i)->getIdWithMotif());
            }
            std::sort( vs.begin(), vs.end() );
            std::string returnString = GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, ' '));
            return ( std::make_tuple ( returnString, firstDict ));            
        }
            
    }
    
    //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
    //static(string ResolvedString, Nullable<double> ResolvedValue) ResolveF2(IEnumerable<double> enumerable);
    std::tuple<std::string, std::optional<double>> PeptideSpectralMatch::ResolveF2(std::vector<double> enumerable)
    {
        //var list = enumerable.ToList();
        //if (list.Max() - list.Min() < ToleranceForDoubleResolutionF2) {
        //	return (list.Average().ToString("F2", CultureInfo.InvariantCulture), list.Average());
        //}
        //else {
        //	var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b =>
        //                   b.ToString("F2", CultreInfo.InvariantCulture))));
        //	return (returnString, nullptr);
        //}
        std::string emptyString = "";
        if ( enumerable.size() == 0 ) {
            return ( std::make_tuple ( emptyString, std::nullopt));
        }

        double max = *std::max_element(enumerable.begin(), enumerable.end());
        double min = *std::min_element(enumerable.begin(), enumerable.end());
        if ( max - min < ToleranceForDoubleResolutionF2 ) {
            double avg = std::accumulate(enumerable.begin(), enumerable.end(), 0.0 )/ enumerable.size();
            return (std::make_tuple(std::to_string(avg),std::make_optional(avg)));
        }
        else {
            std::vector<std::string>vs;
            for ( auto i: enumerable) {
                vs.push_back(std::to_string(i));
            }
            std::string returnString = GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, '|'));
            return ( std::make_tuple ( returnString, std::nullopt ));            
        }
    }
    
    //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
    //static(string ResolvedString, Nullable<double> ResolvedValue) Resolve(IEnumerable<double> enumerable);
    std::tuple<std::string, std::optional<double>> PeptideSpectralMatch::Resolve(std::vector<double> enumerable)
    {
        std::string emptyString = "";
        if ( enumerable.size() == 0 ) {
            return ( std::make_tuple ( emptyString, std::nullopt));
        }

        double max = *std::max_element(enumerable.begin(), enumerable.end());
        double min = *std::min_element(enumerable.begin(), enumerable.end());
        if ( max - min < ToleranceForDoubleResolutionF5 ) {
            double avg = std::accumulate(enumerable.begin(), enumerable.end(), 0.0 )/ enumerable.size();
            return (std::make_tuple(std::to_string(avg),std::make_optional(avg)));
        }
        else {
            std::vector<std::string>vs;
            for ( auto i: enumerable) {
                vs.push_back(std::to_string(i));
            }
            std::string returnString = GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, '|'));
            return ( std::make_tuple ( returnString, std::nullopt ));            
        }
    }
    
    //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
    //static(string ResolvedString, Nullable<int> ResolvedValue) Resolve(IEnumerable<int> enumerable);
    std::tuple<std::string, std::optional<int>>  PeptideSpectralMatch::Resolve(std::vector<int> enumerable)
    {
        //var list = enumerable.ToList();
        //var first = list[0];
        //if (list.All(b => first.Equals(b))) {
        //	return (first.ToString(CultureInfo.InvariantCulture), first);
        //}
        //else{
        //	var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b =>
        //                                                  b.ToString(CultureInfo.InvariantCulture))));
        //	return (returnString, nullptr);
        //}
        std::string emptyString = "";
        if ( enumerable.size() == 0 ) {
            return ( std::make_tuple ( emptyString, std::nullopt));
        }

        int first = enumerable.front();
        bool allfirst=true;
        for ( auto b : enumerable )  {
            if ( b != first ) {
                allfirst = false;
            }            
        }
        if ( allfirst ) {
            return ( std::make_tuple(std::to_string(first), std::make_optional(first)));
        }
        else {
            std::vector<std::string>vs;
            for ( auto i: enumerable) {
                vs.push_back(std::to_string(i));
            }
            std::string returnString = GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, '|'));
            return ( std::make_tuple ( returnString, std::nullopt ));            
        }
    }
    
    //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
    //static(string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable);
    std::tuple<std::string, std::string> PeptideSpectralMatch::Resolve(std::vector<std::string> enumerable)
    {
        std::string emptyString = "";
        if ( enumerable.size() == 0 ) {
            return ( std::make_tuple ( emptyString, emptyString ));
        }
        
        std::string first= enumerable.front() ;
        bool allnull=true, allfirst=true;
        for ( auto b : enumerable )  {
            if ( b != "" ) {
                allnull = false;
            }
            if ( b != first ) {
                allfirst = false;
            }            
        }
        if ( allnull|| allfirst ) {
            // Only first if list is either all null or all equal to the first
            return ( std::make_tuple(first, first) );
        }
        else {
            std::string returnString = GlobalVariables::CheckLengthOfOutput(StringHelper::join(enumerable, '|'));
            return ( std::make_tuple ( returnString, emptyString ));
        }

    }
    
    //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
    //static(string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable, string ambiguousIfNull);
    std::tuple<std::string, std::string> PeptideSpectralMatch::Resolve(std::vector<std::string> enumerable,
                                                                       std::string ambiguousIfNull)
    {

        std::string emptyString = "";
        if ( enumerable.size() == 0 ) {
            return ( std::make_tuple ( emptyString, emptyString ));
        }
        
        std::string first= enumerable.front() ;
        bool allnull=true, allfirst=true;
        for ( auto b : enumerable )  {
            if ( b != "" ) {
                allnull = false;
            }
            if ( b != first ) {
                allfirst = false;
            }            
        }
        if ( allnull|| allfirst ) {
            // Only first if list is either all null or all equal to the first
            return ( std::make_tuple(first, first) );
        }
        else if ( ambiguousIfNull != "" ) {
            std::vector<std::string> uniqueStrings (enumerable.begin(), enumerable.end() );
            std::sort(uniqueStrings.begin(), uniqueStrings.end() );
            auto last = std::unique ( uniqueStrings.begin(), uniqueStrings.end());
            uniqueStrings.erase ( last, uniqueStrings.end());
            std::string returnString = GlobalVariables::CheckLengthOfOutput(StringHelper::join(uniqueStrings, '|'));
            return ( std::make_tuple ( returnString, emptyString ));
        }
        else {
            std::string returnString = GlobalVariables::CheckLengthOfOutput(StringHelper::join(enumerable, '|'));
            return ( std::make_tuple ( returnString, emptyString ));
        }
    }
}



