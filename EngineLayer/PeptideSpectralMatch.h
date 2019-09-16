#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <any>
#include <optional>
#include "stringhelper.h"
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
//namespace EngineLayer { class IScan; }
#include "IScan.h"

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

#include "FdrAnalysis/FdrInfo.h"
using namespace EngineLayer::FdrAnalysis;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    class PeptideSpectralMatch
    {
    private:
        ChemicalFormula *privateModsChemicalFormula;
        std::string privateFullSequence;
        std::optional<int> privateNotch;
        std::string privateBaseSequence;
        std::optional<int> privatePeptideLength;
        std::optional<int> privateOneBasedStartResidueInProtein;
        std::optional<int> privateOneBasedEndResidueInProtein;
        std::optional<double> privatePeptideMonisotopicMass;
        std::optional<int> privateProteinLength;
        std::string privateProteinAccession;
        std::string privateOrganism;
        std::vector<MatchedFragmentIon*> privateMatchedFragmentIons;
        std::unordered_map<std::string, int> privateModsIdentified;
        std::vector<double> privateLocalizedScores;
        int privateScanNumber = 0;
        std::optional<int> privatePrecursorScanNumber;
        double privateScanRetentionTime = 0;
        int privateScanExperimentalPeaks = 0;
        double privateTotalIonCurrent = 0;
        int privateScanPrecursorCharge = 0;
        double privateScanPrecursorMonoisotopicPeakMz = 0;
        double privateScanPrecursorMass = 0;
        std::string privateFullFilePath;
        int privateScanIndex = 0;
        FdrInfo *privateFdrInfo;
        double privateScore = 0;
        double privateDeltaScore = 0;
        double privateRunnerUpScore = 0;
        bool privateIsDecoy = false;
        bool privateIsContaminant = false;
        std::vector<double> privateAllScores;
        std::unordered_map<PeptideWithSetModifications*, std::vector<MatchedFragmentIon*>> privatePeptidesToMatchingFragments;
        
        static constexpr double ToleranceForDoubleResolutionF5 = 1e-6;
        static constexpr double ToleranceForDoubleResolutionF2 = 1e-3;
        std::vector<std::tuple<int, PeptideWithSetModifications*>> _bestMatchingPeptides;
        
    public:
        static constexpr double ToleranceForScoreDifferentiation = 1e-9;
        
        virtual ~PeptideSpectralMatch()
        {
            delete digestionParams;
        }
        
        PeptideSpectralMatch(PeptideWithSetModifications *peptide, int notch, double score, int scanIndex, IScan *scan, DigestionParams *digestionParams, std::vector<MatchedFragmentIon*> &matchedFragmentIons);
        
        // these fields will be null if they are ambiguous
        ChemicalFormula *getModsChemicalFormula() const;
        void setModsChemicalFormula(ChemicalFormula *value);
        std::string getFullSequence() const;
        void setFullSequence(const std::string &value);
        std::optional<int> getNotch() const;
        void setNotch(const std::optional<int> &value);
        std::string getBaseSequence() const;
        void setBaseSequence(const std::string &value);
        std::optional<int> getPeptideLength() const;
        void setPeptideLength(const std::optional<int> &value);
        std::optional<int> getOneBasedStartResidueInProtein() const;
        void setOneBasedStartResidueInProtein(const std::optional<int> &value);
        std::optional<int> getOneBasedEndResidueInProtein() const;
        void setOneBasedEndResidueInProtein(const std::optional<int> &value);
        std::optional<double> getPeptideMonisotopicMass() const;
        void setPeptideMonisotopicMass(const std::optional<double> &value);
        std::optional<int> getProteinLength() const;
        void setProteinLength(const std::optional<int> &value);
        std::string getProteinAccession() const;
        void setProteinAccession(const std::string &value);
        std::string getOrganism() const;
        void setOrganism(const std::string &value);
        std::vector<MatchedFragmentIon*> getMatchedFragmentIons() const;
        void setMatchedFragmentIons(const std::vector<MatchedFragmentIon*> &value);
        
        // these should never be null under normal circumstances
        std::unordered_map<std::string, int> getModsIdentified() const;
        void setModsIdentified(const std::unordered_map<std::string, int> &value);
        std::vector<double> getLocalizedScores() const;
        void setLocalizedScores(const std::vector<double> &value);
        int getScanNumber() const;
        std::optional<int> getPrecursorScanNumber() const;
        double getScanRetentionTime() const;
        int getScanExperimentalPeaks() const;
        double getTotalIonCurrent() const;
        int getScanPrecursorCharge() const;
        double getScanPrecursorMonoisotopicPeakMz() const;
        double getScanPrecursorMass() const;
        std::string getFullFilePath() const;
        int getScanIndex() const;
        int getNumDifferentMatchingPeptides() const;
        FdrInfo *getFdrInfo() const;
        void setFdrInfo(FdrInfo *value);
        double getScore() const;
        void setScore(double value);
        double getDeltaScore() const;
        void setDeltaScore(double value);
        double getRunnerUpScore() const;
        void setRunnerUpScore(double value);
        bool getIsDecoy() const;
        void setIsDecoy(bool value);
        bool getIsContaminant() const;
        void setIsContaminant(bool value);
        DigestionParams *digestionParams;
        std::vector<double> getAllScores() const;
        void setAllScores(const std::vector<double> &value);
        std::unordered_map<PeptideWithSetModifications*, std::vector<MatchedFragmentIon*>> getPeptidesToMatchingFragments() const;
        void setPeptidesToMatchingFragments(const std::unordered_map<PeptideWithSetModifications*, std::vector<MatchedFragmentIon*>> &value);
        
    private:
        //private *IEnumerable < (int Notch, PeptideWithSetModifications *Peptide) BestMatchingPeptides;
        std::vector<std::tuple<int, PeptideWithSetModifications *>> BestMatchingPeptides;

        /// <summary>
        /// Used for Percolator output
        /// </summary>
        std::vector<double> getFeatures() const;
        
        static std::string GetTabSeparatedHeader();
        
        void AddOrReplace(PeptideWithSetModifications *pwsm, double newScore, int notch, bool reportAllAmbiguity,
                          std::vector<MatchedFragmentIon*> &matchedFragmentIons);
        
        std::string ToString();
        //std::string ToString(IReadOnlyDictionary<std::string, int> *ModstoWritePruned);
        std::string ToString(std::unordered_map<std::string, int> *ModstoWritePruned);
        
        
        static std::unordered_map<std::string, std::string> DataDictionary(PeptideSpectralMatch *psm,
                                                                         std::unordered_map<std::string, int> *ModsToWritePruned);
        
        void CalculateDeltaScore(double scoreCutoff);
        
        void SetFdrValues(double cumulativeTarget, double cumulativeDecoy, double qValue, double cumulativeTargetNotch,
                          double cumulativeDecoyNotch, double qValueNotch, double maximumLikelihood,
                          double eValue, double eScore, bool calculateEValue);

        /// <summary>
        /// This method saves properties of this PSM for internal use. It is NOT used for any output.
        /// These resolved fields are (usually) null if there is more than one option.
        /// e.g., if this PSM can be explained by more than one base sequence, the BaseSequence property will be null
        /// </summary>
        void ResolveAllAmbiguities();
        
        
    public:
        /// <summary>
        /// This method is used by protein parsimony to remove PeptideWithSetModifications objects
        /// that have non-parsimonious protein associations
        /// </summary>
        //C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
        void TrimProteinMatches(std::vector<Protein*> parsimoniousProteins);
    
        /// <summary>
        /// This method is used by protein parsimony to add PeptideWithSetModifications objects for
        /// modification-agnostic parsimony
        /// </summary>
        //C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
        void AddProteinMatch(std::tuple<int, PeptideWithSetModifications*> peptideWithNotch);
    
    private:
        /// <summary>
        /// Resolve Methods()
        /// if all 'values' are the same this returns the one value otherwise you get a separated list
        /// of all values in their original order.
        /// for example:
        /// Notches 1,1,1,1 returns as 1
        /// Notches 1,0,1,0 returns as 1|0|1|0
        /// </summary>

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, ChemicalFormula ResolvedValue) Resolve(IEnumerable<IEnumerable<Modification>> enumerable);
        // static std::tuple<std::string, ChemicalFormula*> Resolve(std::vector<std::vector<Modification *>> enumerable);
        static std::tuple<std::string, ChemicalFormula*> Resolve(std::vector<Modification *> enumerable);
        
        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, Dictionary<string, int> ResolvedValue) Resolve(IEnumerable<Dictionary<int, Modification>> enumerable);
        static std::tuple<std::string, std::unordered_map<std::string, int>> Resolve(std::vector<std::unordered_map<int, Modification*>> enumerable);

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, Nullable<double> ResolvedValue) ResolveF2(IEnumerable<double> enumerable);
        static std::tuple<std::string, std::optional<double>> ResolveF2(std::vector<double> enumerable);

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, Nullable<double> ResolvedValue) Resolve(IEnumerable<double> enumerable);
        static std::tuple<std::string, std::optional<double>> Resolve(std::vector<double> enumerable);

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, Nullable<int> ResolvedValue) Resolve(IEnumerable<int> enumerable);
        static std::tuple<std::string, std::optional<int>>  Resolve(std::vector<int> enumerable);
        
        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable);
        static std::tuple<std::string, std::string> Resolve(std::vector<std::string> enumerable);

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable, string ambiguousIfNull);
        static std::tuple<std::string, std::string> Resolve(std::vector<std::string> enumerable, std::string ambiguousIfNull);

        //C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
        //static void AddBasicMatchData(Dictionary<string, string> s, PeptideSpectralMatch psm);
        static void AddBasicMatchData(std::unordered_map<std::string, std::string> &s, PeptideSpectralMatch *psm);

        //C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
        //static void AddPeptideSequenceData(Dictionary<string, string> s, PeptideSpectralMatch psm,
        //                                   IReadOnlyDictionary<string, int> ModsToWritePruned);
        static void AddPeptideSequenceData(std::unordered_map<std::string, std::string> s,
                                           PeptideSpectralMatch *psm,
                                           std::unordered_map<std::string, int> *ModsToWritePruned);

        /// <summary>
        /// Determines whether a peptide includes a splice site
        /// </summary>
        /// <param name="pep"></param>
        /// <param name="site"></param>
        /// <returns></returns>
        //C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
        static bool Includes(PeptideWithSetModifications pep, SpliceSite site);
            
        /// <summary>
        /// Checks for an intersection between a peptide and applied variant that shows a sequence change.
        /// </summary>
        /// <param name="pep"></param>
        /// <param name="appliedVariation"></param>
        /// <returns></returns>
        //C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
        static bool IntersectsWithVariation(PeptideWithSetModifications pep, SequenceVariation appliedVariation,
                                            bool checkUnique);

        /// <summary>
        /// Makes the string representing a detected sequence variation, including any modifications on a variant amino acid
        /// </summary>
        /// <param name="p"></param>
        /// <param name="d"></param>
        /// <returns></returns>
        //C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
        static std::string SequenceVariantString(PeptideWithSetModifications p, SequenceVariation applied);

        //C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
        //static void AddMatchedIonsData(Dictionary<string, string> s, PeptideSpectralMatch psm);
        static void AddMatchedIonsData(std::unordered_map<std::string, std::string> s, PeptideSpectralMatch *psm);

        //C# TO C++ CONVERTER TODO TASK: Local functions are not converted by C# to C++ Converter:
        //static void AddMatchScoreData(Dictionary<string, string> s, PeptideSpectralMatch peptide)
        static void AddMatchScoreData(std::unordered_map<std::string, std::string> s, PeptideSpectralMatch *peptide);

    };
}
