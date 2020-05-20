#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include "stringhelper.h"
#include "stringbuilder.h"

#include "../PeptideSpectralMatch.h"
#include "FlashLFQ/FlashLFQ.h"
using namespace FlashLFQ;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

    typedef std::tuple<int, Modification*> PGroupTuple;

    struct PGroupTuple_hash: public std::unary_function<PGroupTuple, std::size_t>{
        std::size_t operator() (const PGroupTuple& k ) const
        {
            size_t h1= std::hash<int>{}(std::get<0>(k));
            size_t h2= std::hash<int>{}((long)std::get<1>(k));
            return h1 ^ (h2 << 1);
        }
    };

    struct PGroupTuple_equal: public std::binary_function<PGroupTuple, PGroupTuple, bool>{
        bool operator() (const PGroupTuple& lhs, const PGroupTuple& rhs) const
        {
            return std::get<0>(lhs) == std::get<0>(rhs) &&
                std::get<1>(lhs) == std::get<1>(rhs);
        }
    };

    typedef std::unordered_set<PGroupTuple, PGroupTuple_hash, PGroupTuple_equal> PGroupTuple_set;
    

    class ProteinGroup
    {
    private:
        bool privateIsDecoy = false;
        bool privateIsContaminant = false;
        std::vector<SpectraFileInfo*> privateFilesForQuantification;
        std::unordered_set<Protein*> privateProteins;
        std::string privateProteinGroupName;
        double privateProteinGroupScore = 0;
        std::unordered_set<PeptideWithSetModifications*> privateAllPeptides;
        std::unordered_set<PeptideWithSetModifications*> privateUniquePeptides;
        std::unordered_set<PeptideSpectralMatch*> privateAllPsmsBelowOnePercentFDR;
        std::vector<double> privateSequenceCoveragePercent;
        std::vector<std::string> privateSequenceCoverageDisplayList;
        std::vector<std::string> privateSequenceCoverageDisplayListWithMods;
        double privateQValue = 0;
        double privateBestPeptideQValue = 0;
        double privateBestPeptideScore = 0;
        int privateCumulativeTarget = 0;
        int privateCumulativeDecoy = 0;
        bool privateDisplayModsOnPeptides = false;
        std::vector<std::string> privateModsInfo;
        std::unordered_map<SpectraFileInfo*, double> privateIntensitiesByFile;
        
    public:
        ProteinGroup(std::unordered_set<Protein*> &proteins, std::unordered_set<PeptideWithSetModifications*> &peptides,
                     std::unordered_set<PeptideWithSetModifications*> &uniquePeptides);
        
        bool getIsDecoy() const;
        
        bool getIsContaminant() const;
        
        std::vector<SpectraFileInfo*> getFilesForQuantification() const;
        void setFilesForQuantification(const std::vector<SpectraFileInfo*> &value);
        
        std::unordered_set<Protein*> getProteins() const;
        void setProteins(const std::unordered_set<Protein*> &value);
        
        std::string getProteinGroupName() const;
        void setProteinGroupName(const std::string &value);
        
        double getProteinGroupScore() const;
        void setProteinGroupScore(double value);
        
        std::unordered_set<PeptideWithSetModifications*> getAllPeptides() const;
        void setAllPeptides(const std::unordered_set<PeptideWithSetModifications*> &value);
        
        std::unordered_set<PeptideWithSetModifications*> getUniquePeptides() const;
        void setUniquePeptides(const std::unordered_set<PeptideWithSetModifications*> &value);
        
        std::unordered_set<PeptideSpectralMatch*> getAllPsmsBelowOnePercentFDR() const;
        void setAllPsmsBelowOnePercentFDR(const std::unordered_set<PeptideSpectralMatch*> &value);
        
        std::vector<double> getSequenceCoveragePercent() const;
        void setSequenceCoveragePercent(const std::vector<double> &value);
        
        std::vector<std::string> getSequenceCoverageDisplayList() const;
        void setSequenceCoverageDisplayList(const std::vector<std::string> &value);
        
        std::vector<std::string> getSequenceCoverageDisplayListWithMods() const;
        void setSequenceCoverageDisplayListWithMods(const std::vector<std::string> &value);
        
        double getQValue() const;
        void setQValue(double value);
        
        double getBestPeptideQValue() const;
        void setBestPeptideQValue(double value);
        
        double getBestPeptideScore() const;
        void setBestPeptideScore(double value);
        
        int getCumulativeTarget() const;
        void setCumulativeTarget(int value);
        
        int getCumulativeDecoy() const;
        void setCumulativeDecoy(int value);
        
        bool getDisplayModsOnPeptides() const;
        void setDisplayModsOnPeptides(bool value);
        
        std::vector<std::string> getModsInfo() const;
        void setModsInfo(const std::vector<std::string> &value);
        
        std::unordered_map<SpectraFileInfo*, double> getIntensitiesByFile() const;
        void setIntensitiesByFile(const std::unordered_map<SpectraFileInfo*, double> &value);
        
    private:
        std::vector<Protein*> ListOfProteinsOrderedByAccession;
        
    public:
        std::string GetTabSeparatedHeader();
        
        std::string ToString();
        
        // this method is only used internally, to make protein grouping faster
        // this is NOT an output and is NOT used for protein FDR calculations
        void Score();
        
        void CalculateSequenceCoverage();
        
        void MergeProteinGroupWith(ProteinGroup *other);
        
        ProteinGroup *ConstructSubsetProteinGroup(const std::string &fullFilePath);
    };
}
