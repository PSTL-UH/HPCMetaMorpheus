#include "ProteinGroup.h"
#include "../PeptideSpectralMatch.h"
#include "../GlobalVariables.h"
#include "Group.h"

#include <bits/stdc++.h> 
#include <string>

using namespace FlashLFQ;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

    ProteinGroup::ProteinGroup(std::unordered_set<Protein*> proteins,
                               std::unordered_set<PeptideWithSetModifications*> peptides,
                               std::unordered_set<PeptideWithSetModifications*> uniquePeptides)
    {
        setProteins(proteins);
#ifdef ORIG        
        ListOfProteinsOrderedByAccession = getProteins().OrderBy([&] (std::any p)   {
                p::Accession;
            }).ToList();
#endif
        for ( auto p : proteins ) {
            ListOfProteinsOrderedByAccession.push_back(p);
        }
        std::sort(ListOfProteinsOrderedByAccession.begin(), ListOfProteinsOrderedByAccession.end(), [&]
                  (Protein *l, Protein *r ) {
                      return l->getAccession() < r->getAccession();
                  });
#ifdef ORIG
        setProteinGroupName(std::string::Join("|", ListOfProteinsOrderedByAccession.Select([&] (std::any p) {
			p::Accession;
                    })));
#endif
        for ( auto p= ListOfProteinsOrderedByAccession.begin(); p!= ListOfProteinsOrderedByAccession.end(); p++  ) {
            privateProteinGroupName += (*p)->getAccession();
            if ( p+1 != ListOfProteinsOrderedByAccession.end() ) {
                privateProteinGroupName += "|";
            }
        }
        
        setAllPeptides(peptides);
        setUniquePeptides(uniquePeptides);
        setAllPsmsBelowOnePercentFDR(std::unordered_set<PeptideSpectralMatch*>());
        setSequenceCoveragePercent(std::vector<double>());
        setSequenceCoverageDisplayList(std::vector<std::string>());
        setSequenceCoverageDisplayListWithMods(std::vector<std::string>());
        setProteinGroupScore(0);
        setBestPeptideScore(0);
        setQValue(0);
        privateIsDecoy = false;
        privateIsContaminant = false;
        setModsInfo(std::vector<std::string>());
        
        // if any of the proteins in the protein group are decoys, the protein group is a decoy
        for (auto protein : proteins)
        {
            if (protein->getIsDecoy())
            {
                privateIsDecoy = true;
                break;
            }
            if (protein->getIsContaminant())
            {
                privateIsContaminant = true;
                break;
            }
        }
    }
    
    bool ProteinGroup::getIsDecoy() const
    {
        return privateIsDecoy;
    }
    
    bool ProteinGroup::getIsContaminant() const
    {
        return privateIsContaminant;
    }
    
    std::vector<SpectraFileInfo*> ProteinGroup::getFilesForQuantification() const
    {
        return privateFilesForQuantification;
    }
    
    void ProteinGroup::setFilesForQuantification(const std::vector<SpectraFileInfo*> &value)
    {
        privateFilesForQuantification = value;
    }
    
    std::unordered_set<Protein*> ProteinGroup::getProteins() const
    {
        return privateProteins;
    }
    
    void ProteinGroup::setProteins(const std::unordered_set<Protein*> &value)
    {
        privateProteins = value;
    }
    
    std::string ProteinGroup::getProteinGroupName() const
    {
        return privateProteinGroupName;
    }
    
    void ProteinGroup::setProteinGroupName(const std::string &value)
    {
        privateProteinGroupName = value;
    }
    
    double ProteinGroup::getProteinGroupScore() const
    {
        return privateProteinGroupScore;
    }
    
    void ProteinGroup::setProteinGroupScore(double value)
    {
        privateProteinGroupScore = value;
    }
    
    std::unordered_set<PeptideWithSetModifications*> ProteinGroup::getAllPeptides() const
    {
        return privateAllPeptides;
    }
    
    void ProteinGroup::setAllPeptides(const std::unordered_set<PeptideWithSetModifications*> &value)
    {
        privateAllPeptides = value;
    }
    
    std::unordered_set<PeptideWithSetModifications*> ProteinGroup::getUniquePeptides() const
    {
        return privateUniquePeptides;
    }
    
    void ProteinGroup::setUniquePeptides(const std::unordered_set<PeptideWithSetModifications*> &value)
    {
        privateUniquePeptides = value;
    }
    
    std::unordered_set<PeptideSpectralMatch*> ProteinGroup::getAllPsmsBelowOnePercentFDR() const
    {
        return privateAllPsmsBelowOnePercentFDR;
    }
    
    void ProteinGroup::setAllPsmsBelowOnePercentFDR(const std::unordered_set<PeptideSpectralMatch*> &value)
    {
        privateAllPsmsBelowOnePercentFDR = value;
    }
    
    std::vector<double> ProteinGroup::getSequenceCoveragePercent() const
    {
        return privateSequenceCoveragePercent;
    }
    
    void ProteinGroup::setSequenceCoveragePercent(const std::vector<double> &value)
    {
        privateSequenceCoveragePercent = value;
    }
    
    std::vector<std::string> ProteinGroup::getSequenceCoverageDisplayList() const
    {
        return privateSequenceCoverageDisplayList;
    }
    
    void ProteinGroup::setSequenceCoverageDisplayList(const std::vector<std::string> &value)
    {
        privateSequenceCoverageDisplayList = value;
    }
    
    std::vector<std::string> ProteinGroup::getSequenceCoverageDisplayListWithMods() const
    {
        return privateSequenceCoverageDisplayListWithMods;
    }
    
    void ProteinGroup::setSequenceCoverageDisplayListWithMods(const std::vector<std::string> &value)
    {
        privateSequenceCoverageDisplayListWithMods = value;
    }
    
    double ProteinGroup::getQValue() const
    {
        return privateQValue;
    }
    
    void ProteinGroup::setQValue(double value)
    {
        privateQValue = value;
    }
    
    double ProteinGroup::getBestPeptideQValue() const
    {
        return privateBestPeptideQValue;
    }
    
    void ProteinGroup::setBestPeptideQValue(double value)
    {
        privateBestPeptideQValue = value;
    }
    
    double ProteinGroup::getBestPeptideScore() const
    {
        return privateBestPeptideScore;
    }
    
    void ProteinGroup::setBestPeptideScore(double value)
    {
        privateBestPeptideScore = value;
    }
    
    int ProteinGroup::getCumulativeTarget() const
    {
        return privateCumulativeTarget;
    }
    
    void ProteinGroup::setCumulativeTarget(int value)
    {
        privateCumulativeTarget = value;
    }
    
    int ProteinGroup::getCumulativeDecoy() const
    {
        return privateCumulativeDecoy;
    }
    
    void ProteinGroup::setCumulativeDecoy(int value)
    {
        privateCumulativeDecoy = value;
    }
    
    bool ProteinGroup::getDisplayModsOnPeptides() const
    {
        return privateDisplayModsOnPeptides;
    }
    
    void ProteinGroup::setDisplayModsOnPeptides(bool value)
    {
        privateDisplayModsOnPeptides = value;
    }
    
    std::vector<std::string> ProteinGroup::getModsInfo() const
    {
        return privateModsInfo;
    }
    
    void ProteinGroup::setModsInfo(const std::vector<std::string> &value)
    {
        privateModsInfo = value;
    }
    
    std::unordered_map<SpectraFileInfo*, double> ProteinGroup::getIntensitiesByFile() const
    {
        return privateIntensitiesByFile;
    }
    
    void ProteinGroup::setIntensitiesByFile(const std::unordered_map<SpectraFileInfo*, double> &value)
    {
        privateIntensitiesByFile = value;
    }
    
    std::string ProteinGroup::GetTabSeparatedHeader()
    {
        auto sb = new StringBuilder();
        sb->append("Protein Accession" + StringHelper::toString('\t'));
        sb->append("Gene" + StringHelper::toString('\t'));
        sb->append("Organism" + StringHelper::toString('\t'));
        sb->append("Protein Full Name" + StringHelper::toString('\t'));
        sb->append("Protein Unmodified Mass" + StringHelper::toString('\t'));
        sb->append("Number of Proteins in Group" + StringHelper::toString('\t'));
        sb->append("Unique Peptides" + StringHelper::toString('\t'));
        sb->append("Shared Peptides" + StringHelper::toString('\t'));
        sb->append("Number of Peptides" + StringHelper::toString('\t'));
        sb->append("Number of Unique Peptides" + StringHelper::toString('\t'));
        sb->append("Sequence Coverage %" + StringHelper::toString('\t'));
        sb->append("Sequence Coverage" + StringHelper::toString('\t'));
        sb->append("Sequence Coverage with Mods" + StringHelper::toString('\t'));
        sb->append(std::string("Modification Info List") + "\t");
        if (getFilesForQuantification().size() > 0)
        {
            for (int i = 0; i < (int)getFilesForQuantification().size(); i++)
            {
                sb->append("Intensity_" + getFilesForQuantification()[i]->FilenameWithoutExtension + '\t');
            }
        }
        sb->append("Number of PSMs" + StringHelper::toString('\t'));
        sb->append("Protein Decoy/Contaminant/Target" + StringHelper::toString('\t'));
        sb->append("Protein Cumulative Target" + StringHelper::toString('\t'));
        sb->append("Protein Cumulative Decoy" + StringHelper::toString('\t'));
        sb->append("Protein QValue" + StringHelper::toString('\t'));
        sb->append("Best Peptide Score" + StringHelper::toString('\t'));
        sb->append("Best Peptide Notch QValue");
        
        std::string s =  sb->toString();
        delete sb;
        return s;
    }
    
    std::string ProteinGroup::ToString()
    {
        auto sb = new StringBuilder();
        
        // list of protein accession numbers
        sb->append(getProteinGroupName());
        sb->append("\t");
        
        // genes
#ifdef ORIG
        sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", ListOfProteinsOrderedByAccession.Select([&]
                                                                                                                 (std::any p) {
                            p::GeneNames->Select([&] (std::any x) {
                                    x::Item2;
                                }).FirstOrDefault();
                                                                                                                       }))
                       ));
#endif
        std::string s;
        for ( auto p= ListOfProteinsOrderedByAccession.begin(); p != ListOfProteinsOrderedByAccession.end(); p++ ){
            s += std::get<1>((*p)->getGeneNames()[0]);
            if ( p+1 != ListOfProteinsOrderedByAccession.end() ){
                s+= "|";
            }
        }
        sb->append(GlobalVariables::CheckLengthOfOutput(s));
        sb->append("\t");
        
        // organisms
#ifdef ORIG
        sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", ListOfProteinsOrderedByAccession.Select([&]
                                                                                                             ( std::any p)	{
                 p::Organism;
        }).Distinct())));
#endif
        std::vector<std::string> vs;
        for ( auto p= ListOfProteinsOrderedByAccession.begin(); p != ListOfProteinsOrderedByAccession.end(); p++ ){
            std::string tmps = (*p)->getOrganism();
            bool found = false;
            for ( auto x: vs ) {
                if ( tmps == x ) {
                    found = true;
                    break;
                }
            }
            if ( !found  ) {
                vs.push_back(tmps);
            }
        }
        std::string ps = "|";
        sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, ps)));
        sb->append("\t");
        
        // list of protein names
#ifdef ORIG
        sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", ListOfProteinsOrderedByAccession.Select([&] (std::any p)	{
                            p->FullName;
                        }).Distinct())));
#endif
        vs.clear();
        for ( auto p : ListOfProteinsOrderedByAccession ){
            std::string tmps = p->getFullName();
            bool found = false;
            for ( auto x: vs ) {
                if ( tmps == x ) {
                    found = true;
                    break;
                }
            }
            if ( !found  ) {
                vs.push_back(tmps);
            }
        }
        sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, ps)));
        sb->append("\t");
        
        // list of masses
#ifdef ORIG
        auto sequences = ListOfProteinsOrderedByAccession.Select([&] (std::any p)   {
                p::BaseSequence;
            }).Distinct();
#endif
        std::vector<std::string> sequences;
        for ( auto p : ListOfProteinsOrderedByAccession ){
            std::string tmps = p->getBaseSequence();
            bool found = false;
            for ( auto x: sequences ) {
                if ( tmps == x ) {
                    found = true;
                    break;
                }
            }
            if ( !found  ) {
                sequences.push_back(tmps);
            }
        }
        

        std::vector<double> masses;
        std::vector<std::string> massesString;
        for (auto sequence : sequences)
        {
            try
            {
                auto  tempVar = new Proteomics::AminoAcidPolymer::Peptide(sequence);
                masses.push_back(tempVar->getMonoisotopicMass());
                massesString.push_back(std::to_string(tempVar->getMonoisotopicMass()));
            }
            catch (const std::runtime_error &e1)
            {
                masses.push_back(NAN);
                massesString.push_back(std::to_string(NAN));
            }
        }
        sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join( massesString, ps)));
        sb->append("\t");
        
        // number of proteins in group
        sb->append("" + std::to_string(getProteins().size()));
        sb->append("\t");
        
        // list of unique peptides
        int uniquePeptidescount=0;
        if (!getDisplayModsOnPeptides())
        {
#ifdef ORIG
            sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", getUniquePeptides().Select([&] (std::any p) {
				p::BaseSequence;
                            }).Distinct())));
#endif
            vs.clear();
            for ( auto p : getUniquePeptides() ){
                std::string tmps = p->getBaseSequence();
                bool found = false;
                for ( auto x: vs ) {
                    if ( tmps == x ) {
                        found = true;
                        break;
                    }
                }
                if ( !found  ) {
                    vs.push_back(tmps);
                }
            }
            uniquePeptidescount = vs.size();
            sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, ps )));
        }
        else
        {
#ifdef ORIG
            sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", getUniquePeptides().Select([&] (std::any p)  {
				p::FullSequence;
                            }).Distinct())));
#endif
            vs.clear();
            for ( auto p : getUniquePeptides() ){
                std::string tmps = p->getFullSequence();
                bool found = false;
                for ( auto x: vs ) {
                    if ( tmps == x ) {
                        found = true;
                        break;
                    }
                }
                if ( !found  ) {
                    vs.push_back(tmps);
                }
            }
            uniquePeptidescount = vs.size();
            sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, ps )));            
        }
        sb->append("\t");
        
        // list of shared peptides
#ifdef ORIG
        auto SharedPeptides = getAllPeptides().Except(getUniquePeptides());
#endif
        auto SharedPeptides = getAllPeptides();
        for ( auto p : getUniquePeptides() ) {
            for ( auto v: SharedPeptides ) {
                if ( p == v ) {
                    SharedPeptides.erase(v);
                    break;
                }
            }
        }

        if (!getDisplayModsOnPeptides())
        {
#ifdef ORIG
            sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", SharedPeptides->Select([&] (std::any p) {
				p::BaseSequence;
                            }).Distinct())));

#endif
            vs.clear();
            for ( auto p : SharedPeptides ){
                std::string tmps = p->getBaseSequence();
                bool found = false;
                for ( auto x: vs ) {
                    if ( tmps == x ) {
                        found = true;
                        break;
                    }
                }
                if ( !found  ) {
                    vs.push_back(tmps);
                }
            }
            sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, ps )));
        }
        else
        {
#ifdef ORIG
            sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", SharedPeptides->Select([&] (std::any p) {
				p::FullSequence;
                            }).Distinct())));
#endif
            vs.clear();
            for ( auto p : SharedPeptides ){
                std::string tmps = p->getFullSequence();
                bool found = false;
                for ( auto x: vs ) {
                    if ( tmps == x ) {
                        found = true;
                        break;
                    }
                }
                if ( !found  ) {
                    vs.push_back(tmps);
                }
            }
            sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, ps )));
            
        }
        sb->append("\t");
        
        // number of peptides
        if (!getDisplayModsOnPeptides())
        {
#ifdef ORIG
            sb->append("" + getAllPeptides().Select([&] (std::any p) {
                        p::BaseSequence;
                    }).Distinct()->Count());
#endif
            vs.clear();
            for ( auto p : getAllPeptides() ){
                std::string tmps = p->getBaseSequence();
                bool found = false;
                for ( auto x: vs ) {
                    if ( tmps == x ) {
                        found = true;
                        break;
                    }
                }
                if ( !found  ) {
                    vs.push_back(tmps);
                }
            }
            sb->append(std::to_string(vs.size()));            
        }
        else
        {
#ifdef ORIG
            sb->append("" + getAllPeptides().Select([&] (std::any p){
                        p::FullSequence;
                    }).Distinct()->Count());
#endif
            vs.clear();
            for ( auto p : getAllPeptides() ){
                std::string tmps = p->getFullSequence();
                bool found = false;
                for ( auto x: vs ) {
                    if ( tmps == x ) {
                        found = true;
                        break;
                    }
                }
                if ( !found  ) {
                    vs.push_back(tmps);
                }
            }
            sb->append(std::to_string(vs.size()));            
        }
        sb->append("\t");
        
        // number of unique peptides
#ifdef ORIG
        if (!getDisplayModsOnPeptides())
        {
            sb->append("" + getUniquePeptides().Select([&] (std::any p) {
                        p::BaseSequence;
                    }).Distinct()->Count());
        }
        else
        {
            sb->append("" + getUniquePeptides().Select([&] (std::any p)	{
                        p::FullSequence;
                    }).Distinct()->Count());
        }
#endif
        // determined the count already in the previous handling of uniquePeptides
        sb->append(std::to_string(uniquePeptidescount));
        sb->append("\t");
        
        // sequence coverage percent
#ifdef ORIG
        sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", getSequenceCoveragePercent().Select([&] (std::any p) {
                            std::string::Format(std::string("{0:0}") + "%", (p * 100));
                        }))));
#endif
        vs.clear();
        for ( auto p: getSequenceCoveragePercent() ) {
            std::string tmps = StringHelper::formatSimple("{0:0}", (p*100)) + "%";
            vs.push_back(tmps);
        }
        sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(vs, ps )));
        sb->append("\t");
        
        
        // sequence coverage
#ifdef ORIG
        sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", getSequenceCoverageDisplayList())));
#endif
        sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(getSequenceCoverageDisplayList(), ps)));        
	sb->append("\t");

        // sequence coverage with mods
#ifdef ORIG
        sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", getSequenceCoverageDisplayListWithMods())));
#endif
        sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(getSequenceCoverageDisplayListWithMods(), ps)));
        sb->append("\t");

        //Detailed mods information list
#ifdef ORIG
        sb->append(GlobalVariables::CheckLengthOfOutput(std::string::Join("|", getModsInfo())));
#endif
        sb->append(GlobalVariables::CheckLengthOfOutput(StringHelper::join(getModsInfo(), ps)));
        sb->append("\t");
        
        // MS1 intensity (retrieved from FlashLFQ in the SearchTask)
        if (getIntensitiesByFile().size() > 0 && getFilesForQuantification().size() > 0)
        {
            for (auto file : getFilesForQuantification())
            {
                if (getIntensitiesByFile()[file] > 0)
                {
                    sb->append(getIntensitiesByFile()[file]);
                }
                else
                {
                    sb->append("");
                }
                sb->append("\t");
            }
        }
        
        // number of PSMs for listed peptides
        sb->append("" + std::to_string(getAllPsmsBelowOnePercentFDR().size()));
        sb->append("\t");
        
        // isDecoy
        if (getIsDecoy())
        {
            sb->append("D");
        }
        else if (getIsContaminant())
        {
            sb->append("C");
        }
        else
        {
            sb->append("T");
        }
        sb->append("\t");
        
        // cumulative target
        sb->append(getCumulativeTarget());
        sb->append("\t");
        
        // cumulative decoy
        sb->append(getCumulativeDecoy());
        sb->append("\t");
        
        // q value
        sb->append(getQValue());
        sb->append("\t");
        
        // best peptide score
        sb->append(getBestPeptideScore());
        sb->append("\t");
        
        // best peptide q value
        sb->append(getBestPeptideQValue());
        sb->append("\t");
        
        delete sb;
        return sb->toString();
    }
    
    void ProteinGroup::Score()
    {
        // sum the scores of the best PSM per base sequence
#ifdef ORIG
        setProteinGroupScore(getAllPsmsBelowOnePercentFDR().GroupBy([&] (std::any p)   {
                    p::BaseSequence;
                })->Select([&] (std::any p)  {
                        p->Select([&] (std::any x)  {
                                x::Score;
                            }).Max();
                    }).Sum());
#endif
        std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f1 = [&](PeptideSpectralMatch *l, PeptideSpectralMatch *r) {return l->getBaseSequence() < r->getBaseSequence(); } ;
        std::function<bool(PeptideSpectralMatch*,PeptideSpectralMatch*)> f2 = [&](PeptideSpectralMatch *l, PeptideSpectralMatch *r) {return l->getBaseSequence() != r->getBaseSequence(); } ;
        std::vector<PeptideSpectralMatch*> tmpinput;
        for ( auto t : getAllPsmsBelowOnePercentFDR() ) {
            tmpinput.push_back(t);
        }
        std::vector<std::vector<PeptideSpectralMatch*>> tmpres = Group::GroupBy ( tmpinput, f1, f2);
        double sum = 0.0;
        for ( auto t : tmpres ) {
            double max = t[0]->getScore();
            for ( auto p: t ) {
                if ( p->getScore() > max ) max = p->getScore();
            }
            sum += max;
        }
        setProteinGroupScore(sum);
    }
    
    void ProteinGroup::CalculateSequenceCoverage()
    {
        auto proteinsWithUnambigSeqPsms = std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>();
        auto proteinsWithPsmsWithLocalizedMods = std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>();
        
        for (auto protein : getProteins())
        {
            proteinsWithUnambigSeqPsms.emplace(protein, std::vector<PeptideWithSetModifications*>());
            proteinsWithPsmsWithLocalizedMods.emplace(protein, std::vector<PeptideWithSetModifications*>());
        }
        
        for (auto psm : getAllPsmsBelowOnePercentFDR())
        {
            // null BaseSequence means that the amino acid sequence is ambiguous; do not use these to calculate sequence coverage
            if (psm->getBaseSequence().length() != 0)
            {
#ifdef ORIG
                auto peptides = psm->BestMatchingPeptides->Select([&] (std::any p) {
                        p::Peptide;
                    });
#endif
                std::vector<PeptideWithSetModifications*> peptides;
                for ( auto p : psm->getBestMatchingPeptides() ) {
                    peptides.push_back (std::get<1>(p));                    
                }
                
                for (auto peptide : peptides)
                {
                    // might be unambiguous but also shared; make sure this protein group contains this peptide+protein combo
                    auto Proteins = getProteins();
                    if (std::find(Proteins.begin(), Proteins.end(), peptide->getProtein()) != Proteins.end())
                    {
                        proteinsWithUnambigSeqPsms[peptide->getProtein()].push_back(peptide);
                        
                        // null FullSequence means that mods were not successfully localized; do not display them on
                        // the sequence coverage mods info
                        if (psm->getFullSequence().length() !=  0 )
                        {
                            proteinsWithPsmsWithLocalizedMods[peptide->getProtein()].push_back(peptide);
                        }
                    }
                }
            }
        }
        
        for (auto protein : ListOfProteinsOrderedByAccession)
        {
            bool errorResult = false;
            auto sequenceCoverageDisplay = StringHelper::toLower(protein->getBaseSequence());
            std::unordered_set<int> coveredOneBasedResidues;
            
            // get residue numbers of each peptide in the protein and identify them as observed if the sequence is unambiguous
            for (auto peptide : proteinsWithUnambigSeqPsms[protein])
            {
                std::string sequenceExtractedFromProtein = "";
                for (int i = peptide->getOneBasedStartResidueInProtein(); i <= peptide->getOneBasedEndResidueInProtein(); i++)
                {
                    // check for bugs in sequence coverage; make sure we have the right amino acids!
                    sequenceExtractedFromProtein += sequenceCoverageDisplay[i - 1];
                    coveredOneBasedResidues.insert(i);
                }
                
                if (StringHelper::toUpper(sequenceExtractedFromProtein) != peptide->getBaseSequence())
                {
                    errorResult = true;
                }
            }
            
            // calculate sequence coverage percent
            double seqCoveragePercent = static_cast<double>(coveredOneBasedResidues.size()) / protein->getLength();
            if (seqCoveragePercent > 1)
            {
                errorResult = true;
            }
            
            // add the percent coverage or NaN if there was an error
            if (!errorResult)
            {
                getSequenceCoveragePercent().push_back(seqCoveragePercent);
            }
            else
            {
                getSequenceCoveragePercent().push_back(NAN);
            }
            
            // convert the observed amino acids to upper case if they are unambiguously observed
#ifdef ORIG
            auto coverageArray = sequenceCoverageDisplay->ToCharArray();
#endif
            char coverageArray[sequenceCoverageDisplay.length() + 1 ];
            strcpy (coverageArray, sequenceCoverageDisplay.c_str() );
            
            for (auto obsResidueLocation : coveredOneBasedResidues)
            {
                coverageArray[obsResidueLocation - 1] = std::toupper(coverageArray[obsResidueLocation - 1]);
            }
            sequenceCoverageDisplay = std::string(coverageArray);
            
            // check to see if there was an errored result; if not, add the coverage display
            if (!errorResult)
            {
                getSequenceCoverageDisplayList().push_back(sequenceCoverageDisplay);
            }
            else
            {
                getSequenceCoverageDisplayList().push_back("Error calculating sequence coverage");
            }
            
            // put mods in the sequence coverage display
            if (!errorResult)
            {
                // get mods to display in sequence (only unambiguously identified mods)
#ifdef ORIG
                auto modsOnThisProtein = std::unordered_set<KeyValuePair<int, Modification*>*>();
#endif
                PGroupTuple_set modsOnThisProtein;
                for (auto pep : proteinsWithPsmsWithLocalizedMods[protein])
                {
                    for (auto mod : pep->getAllModsOneIsNterminus())
                    {
                        std::string modType = std::get<1>(mod)->getModificationType();
                        if (modType.find("PeptideTermMod")  == std::string::npos  &&
                            modType.find("Common Variable") == std::string::npos &&
                            modType.find("Common Fixed") == std::string::npos )
                        {
                            modsOnThisProtein.insert(std::make_tuple(pep->getOneBasedStartResidueInProtein()
                                                                     + std::get<0>(mod)-2,
                                                                     std::get<1>(mod)));
                        }
                    }
                }
                
#ifdef ORIG
                auto temp1 = modsOnThisProtein.OrderBy([&] (std::any p)   {
                        p::Key;
                    }).ToList();
#endif
                std::vector<PGroupTuple> temp1;
                for ( auto p = modsOnThisProtein.begin(); p != modsOnThisProtein.end(); p++ ) {
                    temp1.push_back(*p);
                }
                std::sort(temp1.begin(), temp1.end(), [&] (PGroupTuple l, PGroupTuple r) {
                        return std::get<0>(l) < std::get<0>(r);
                    });
                
                for (auto mod : temp1)
                {
                    if (std::get<1>(mod)->getLocationRestriction() == "N-terminal.")
                    {
#ifdef ORIG
                        sequenceCoverageDisplay = sequenceCoverageDisplay->Insert(0, "[" +
                                                                                  std::get<1>(mod)->getIdWithMotif() + "]-");
#endif
                        sequenceCoverageDisplay = "[" + std::get<1>(mod)->getIdWithMotif() + "]-" + sequenceCoverageDisplay;
                    }
                    else if (std::get<1>(mod)->getLocationRestriction() == "Anywhere.")
                    {
                        int modStringIndex = sequenceCoverageDisplay.length() - (protein->getLength() - std::get<0>(mod));
#ifdef ORIG
                        sequenceCoverageDisplay = sequenceCoverageDisplay->Insert(modStringIndex, "[" +
                                                                                  std::get<1>(mod)->getIdWithMotif() + "]");
#endif
                        sequenceCoverageDisplay == sequenceCoverageDisplay.substr(0, modStringIndex) +
                            "[" + std::get<1>(mod)->getIdWithMotif() + "]" + sequenceCoverageDisplay.substr(modStringIndex);
                    }
                    else if (std::get<1>(mod)->getLocationRestriction() == "C-terminal." )
                    {
#ifdef ORIG
                        sequenceCoverageDisplay = sequenceCoverageDisplay->Insert(sequenceCoverageDisplay.length(), "-[" +
                                                                                  std::get<1>(mod)->getIdWithMotif() + "]");
#endif
                        sequenceCoverageDisplay += "-[" + std::get<1>(mod)->getIdWithMotif() + "]";
                    }
                }
                
                getSequenceCoverageDisplayListWithMods().push_back(sequenceCoverageDisplay);
                
                if (!modsOnThisProtein.empty())
                {
                    // calculate spectral count percentage of modified observation
                    std::string tempModStrings = ""; //The whole string
                    std::vector<int> tempPepModTotals; //The List of (For one mod, The Modified Pep Num)
                    std::vector<int> tempPepTotals; //The List of (For one mod, The total Pep Num)
                    std::vector<std::string> tempPepModValues; //The List of (For one mod, the Modified Name)
                    std::vector<int> tempModIndex; //The Index of the modified position.
                    
                    for (auto pep : proteinsWithPsmsWithLocalizedMods[protein])
                    {
                        for (auto mod : pep->getAllModsOneIsNterminus())
                        {
                            int tempPepNumTotal = 0; //For one mod, The total Pep Num
                            std::string modType = std::get<1>(mod)->getModificationType();
                            auto  locRestr = std::get<1>(mod)->getLocationRestriction();
                            if (modType.find("Common Variable") == std::string::npos   &&
                                modType.find("Common Fixed") == std::string::npos &&
                                //locRestr != ModLocationOnPeptideOrProtein::PepC  &&
                                //locRestr != ModLocationOnPeptideOrProtein::NPep )
                                locRestr != "PepC"  &&
                                locRestr != "NPep" )

                                {
                                int tempIndexInProtein;
                                if (locRestr == "N-terminal.")
                                {
                                    tempIndexInProtein = 1;
                                }
                                else if (locRestr == "Anywhere." )
                                {
                                    tempIndexInProtein = pep->getOneBasedStartResidueInProtein() + std::get<0>(mod) - 2;
                                }
                                else if (locRestr == "C-terminal.")
                                {
                                    tempIndexInProtein = protein->getLength();
                                }
                                else
                                {
                                    // In case it's a peptide terminal mod, skip!
                                    // we don't want this annotated in the protein's modifications
                                    continue;
                                }

                                auto tempi = std::find(tempModIndex.begin(), tempModIndex.end(), tempIndexInProtein);
                                if ( tempi != tempModIndex.end()                                &&
                                    tempPepModValues[*tempi] == std::get<1>(mod)->getIdWithMotif())
                                {
                                    tempPepModTotals[*tempi] += 1;
                                }
                                else
                                {
                                    tempModIndex.push_back(tempIndexInProtein);
                                    for (auto pept : proteinsWithPsmsWithLocalizedMods[protein])
                                    {
                                        if (tempIndexInProtein >= pept->getOneBasedStartResidueInProtein() - (tempIndexInProtein == 1 ? 1 : 0) && tempIndexInProtein <= pept->getOneBasedEndResidueInProtein())
                                        {
                                            tempPepNumTotal += 1;
                                        }
                                    }
                                    tempPepTotals.push_back(tempPepNumTotal);
                                    tempPepModValues.push_back(std::get<1>(mod)->getIdWithMotif());
                                    tempPepModTotals.push_back(1);
                                }
                            }
                        }
                    }
                    for (int i = 0; i < (int)tempPepModTotals.size(); i++)
                    {
                        
                        std::string tempString = ("#aa" + std::to_string(tempModIndex[i]) + "[" +
                                                  tempPepModValues[i] + ",info:occupancy=" +
                                                  std::to_string(static_cast<double>(tempPepModTotals[i]) /
                                                                 static_cast<double>(tempPepTotals[i])) + "(" +
                                                  std::to_string(tempPepModTotals[i]) + "/" +
                                                  std::to_string(tempPepTotals[i]) + ")" + "];");
                        tempModStrings += tempString;
                    }
                    
                    if (!tempModStrings.empty())
                    {
                        getModsInfo().push_back(tempModStrings);
                    }
                }
            }
            
            //C# TO C++ CONVERTER TODO TASK: A 'delete sequenceCoverageDisplay' statement was not added since
            //sequenceCoverageDisplay was passed to a method or constructor. Handle memory management manually.
        }
    }
    
    void ProteinGroup::MergeProteinGroupWith(ProteinGroup *other)
    {
#ifdef ORIG
        this->getProteins().UnionWith(other->getProteins());
#endif
        auto tmp1 = this->getProteins();
        for ( auto p : other->getProteins() ) {
            bool found = false;
            for ( auto v: tmp1 ) {
                if ( p == v ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                tmp1.insert(p);
            }
        }
        this->setProteins(tmp1);
        
#ifdef ORIG
        this->getAllPeptides().UnionWith(other->getAllPeptides());
#endif
        auto tmp2 = this->getAllPeptides();
        for ( auto p : other->getAllPeptides() ) {
            bool found = false;
            for ( auto v: tmp2 ) {
                if ( p == v ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                tmp2.insert(p);
            }
        }
        this->setAllPeptides(tmp2);
        
#ifdef ORIG
        this->getUniquePeptides().UnionWith(other->getUniquePeptides());
#endif
        auto tmp3 = this->getUniquePeptides();
        for ( auto p : other->getUniquePeptides() ) {
            bool found = false;
            for ( auto v: tmp3 ) {
                if ( p == v ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                tmp3.insert(p);
            }
        }
        this->setUniquePeptides(tmp3);

#ifdef ORIG        
        this->getAllPsmsBelowOnePercentFDR().UnionWith(other->getAllPsmsBelowOnePercentFDR());
#endif
        auto tmp4 = this->getAllPsmsBelowOnePercentFDR();
        for ( auto p : other->getAllPsmsBelowOnePercentFDR() ) {
            bool found = false;
            for ( auto v: tmp4 ) {
                if ( p == v ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                tmp4.insert(p);
            }
        }
        this->setAllPsmsBelowOnePercentFDR(tmp4);
        
        other->setProteinGroupScore(0);
        
#ifdef ORIG
        ListOfProteinsOrderedByAccession = getProteins().OrderBy([&] (std::any p)       {
                p::Accession;
            }).ToList();
#endif
        ListOfProteinsOrderedByAccession.clear();
        for ( auto p: getProteins() ){
            ListOfProteinsOrderedByAccession.push_back(p);
        }
        std::sort(ListOfProteinsOrderedByAccession.begin(), ListOfProteinsOrderedByAccession.end(), [&]
                  (Protein *l, Protein *r) {
                      return l->getAccession() < r->getAccession();
                  });
        
#ifdef ORIG
        setProteinGroupName(std::string::Join("|", ListOfProteinsOrderedByAccession.Select([&] (std::any p)  {
                        p::Accession;
                    })));
#endif
        std::vector<std::string>vs;
        for ( auto p: ListOfProteinsOrderedByAccession ) {
            vs.push_back(p->getAccession());
        }
        std::string ps = "|";
        setProteinGroupName(StringHelper::join(vs, ps));
    }

    ProteinGroup *ProteinGroup::ConstructSubsetProteinGroup(const std::string &fullFilePath)
    {
#ifdef ORIg
        auto allPsmsForThisFile = std::unordered_set<PeptideSpectralMatch*>(this->getAllPsmsBelowOnePercentFDR().Where([&] (std::any p)	{
                    p::FullFilePath->Equals(fullFilePath);
		}));
#endif
        std::unordered_set<PeptideSpectralMatch*> allPsmsForThisFile;
        for ( auto p: this->getAllPsmsBelowOnePercentFDR() ) {
            if (p->getFullFilePath() == fullFilePath  ) {
                allPsmsForThisFile.insert(p);
            }
        }

#ifdef ORIG        
        auto allPeptidesForThisFile = std::unordered_set<PeptideWithSetModifications*>(allPsmsForThisFile.SelectMany([&] (std::any p){
                    p::BestMatchingPeptides->Select([&] (std::any v)  {
                            v::Peptide;
			});
		}));
#endif
        std::unordered_set<PeptideWithSetModifications*>  allPeptidesForThisFile;
        for ( auto p: allPsmsForThisFile ) {
            for ( auto v: p->getBestMatchingPeptides() ) {
                allPeptidesForThisFile.insert(std::get<1>(v) );                
            }
        }

#ifdef ORIG
        auto allUniquePeptidesForThisFile = std::unordered_set<PeptideWithSetModifications*>(this->getUniquePeptides().Intersect(allPeptidesForThisFile));
#endif
        std::unordered_set<PeptideWithSetModifications*> allUniquePeptidesForThisFile;
        for ( auto p: this->getUniquePeptides() ) {
            bool found = false;
            for ( auto v: allPeptidesForThisFile ) {
                if ( p == v ) {
                    found = true;
                    break;
                }
            }
            if ( found ) {
                allUniquePeptidesForThisFile.insert(p);
            }
        }

        auto gP = this->getProteins();
        ProteinGroup *subsetPg = new ProteinGroup(gP, allPeptidesForThisFile, allUniquePeptidesForThisFile);
        subsetPg->setAllPsmsBelowOnePercentFDR(allPsmsForThisFile);
        subsetPg->setDisplayModsOnPeptides(this->getDisplayModsOnPeptides());
        
        SpectraFileInfo *spectraFileInfo = nullptr;
        if (getFilesForQuantification().size() > 0)
        {
#ifdef ORIG
            spectraFileInfo = getFilesForQuantification().Where([&] (std::any p) {
                    return p->FullFilePathWithExtension == fullFilePath;
                }).First();
#endif
            for ( auto p : getFilesForQuantification() ) {
                if ( p->FullFilePathWithExtension == fullFilePath ) {
                    spectraFileInfo = p;
                    break;
                }
            }
            subsetPg->setFilesForQuantification(std::vector<SpectraFileInfo*> {spectraFileInfo});
        }
        
        if (getIntensitiesByFile().empty())
        {
            subsetPg->getIntensitiesByFile().clear();
        }
        else
        {
            subsetPg->setIntensitiesByFile(std::unordered_map<SpectraFileInfo*, double>  {
                    {spectraFileInfo, getIntensitiesByFile()[spectraFileInfo]}
                });
        }
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete subsetPg' statement was not added since subsetPg was used in
        //a 'return' or 'throw' statement.
        return subsetPg;
    }
}
