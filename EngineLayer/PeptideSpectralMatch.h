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
namespace EngineLayer { class IScan; }

using namespace Chemistry;
using namespace EngineLayer::FdrAnalysis;
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

		std::vector<(int Notch, PeptideWithSetModifications Pwsm)*> _bestMatchingPeptides;

	public:
		static constexpr double ToleranceForScoreDifferentiation = 1e-9;

		virtual ~PeptideSpectralMatch()
		{
			delete DigestionParams;
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
		DigestionParams *const DigestionParams;
		std::vector<double> getAllScores() const;
		void setAllScores(const std::vector<double> &value);
		std::unordered_map<PeptideWithSetModifications*, std::vector<MatchedFragmentIon*>> getPeptidesToMatchingFragments() const;
		void setPeptidesToMatchingFragments(const std::unordered_map<PeptideWithSetModifications*, std::vector<MatchedFragmentIon*>> &value);

		private *IEnumerable < (int Notch, PeptideWithSetModifications *Peptide);
	}

		/// <summary>
		/// Used for Percolator output
		/// </summary>
		std::vector<double> getFeatures() const;

		static std::string GetTabSeparatedHeader();

		void AddOrReplace(PeptideWithSetModifications *pwsm, double newScore, int notch, bool reportAllAmbiguity, std::vector<MatchedFragmentIon*> &matchedFragmentIons);

		std::string ToString() override;

		std::string ToString(IReadOnlyDictionary<std::string, int> *ModstoWritePruned);

		static std::unordered_map<std::string, std::string> DataDictionary(PeptideSpectralMatch *psm, IReadOnlyDictionary<std::string, int> *ModsToWritePruned);

		void CalculateDeltaScore(double scoreCutoff);

		void SetFdrValues(double cumulativeTarget, double cumulativeDecoy, double qValue, double cumulativeTargetNotch, double cumulativeDecoyNotch, double qValueNotch, double maximumLikelihood, double eValue, double eScore, bool calculateEValue);

		/// <summary>
		/// This method saves properties of this PSM for internal use. It is NOT used for any output.
		/// These resolved fields are (usually) null if there is more than one option.
		/// e.g., if this PSM can be explained by more than one base sequence, the BaseSequence property will be null
		/// </summary>
		void ResolveAllAmbiguities();

			// TODO: technically, different peptide options for this PSM can have different matched ions
			// we can write a Resolve method for this if we want...
			setMatchedFragmentIons(getPeptidesToMatchingFragments().First()->Value);
};

		/// <summary>
		/// This method is used by protein parsimony to remove PeptideWithSetModifications objects that have non-parsimonious protein associations
		/// </summary>
//C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
		public void TrimProteinMatches(HashSet<Protein> parsimoniousProteins)

		/// <summary>
		/// This method is used by protein parsimony to add PeptideWithSetModifications objects for modification-agnostic parsimony
		/// </summary>
//C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
		public void AddProteinMatch((int, PeptideWithSetModifications) peptideWithNotch)

		/// <summary>
		/// Resolve Methods()
		/// if all 'values' are the same this returns the one value otherwise you get a separated list of all values in their original order.
		/// for example:
		/// Notches 1,1,1,1 returns as 1
		/// Notches 1,0,1,0 returns as 1|0|1|0
		/// </summary>
//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//		private static(string ResolvedString, ChemicalFormula ResolvedValue) Resolve(IEnumerable<IEnumerable<Modification>> enumerable)
//		{
//			var list = enumerable.ToList();
//			ChemicalFormula firstChemFormula = new ChemicalFormula();
//			foreach (var firstMods in list[0])
//			{
//				if (firstMods == nullptr || firstMods.ChemicalFormula == nullptr)
//				{
//					return ("unknown", nullptr);
//				}
//				firstChemFormula.Add(firstMods.ChemicalFormula);
//			}
//
//			bool equals = true;
//			List<ChemicalFormula> formulas = new List<ChemicalFormula>();
//			foreach (var anEnum in list)
//			{
//				ChemicalFormula fhere = new ChemicalFormula();
//				foreach (var mod in anEnum)
//				{
//					if (mod == nullptr || mod.ChemicalFormula == nullptr)
//					{
//						return ("unknown", nullptr);
//					}
//					fhere.Add(mod.ChemicalFormula);
//				}
//				if (!firstChemFormula.Equals(fhere))
//				{
//					equals = false;
//				}
//				formulas.Add(fhere);
//			}
//			if (!equals)
//			{
//				var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", formulas.Select(b => b.Formula)));
//				return (returnString, nullptr);
//			}
//			else
//			{
//				return (firstChemFormula.Formula, firstChemFormula);
//			}
//		}

//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//		private static(string ResolvedString, Dictionary<string, int> ResolvedValue) Resolve(IEnumerable<Dictionary<int, Modification>> enumerable)
//		{
//			var list = enumerable.ToList();
//			Dictionary<string, int> firstDict = list[0].Values.OrderBy(b => b.IdWithMotif).GroupBy(b => b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());
//
//			bool equals = true;
//			foreach (var dict in list)
//			{
//				Dictionary<string, int> okTest = dict.Values.OrderBy(b => b.IdWithMotif).GroupBy(b => b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());
//				if (!firstDict.SequenceEqual(okTest))
//				{
//					equals = false;
//					break;
//				}
//			}
//			if (!equals)
//			{
//				var returnString = string.Join("|", list.Select(b => string.Join(" ", b.Values.Select(c => c.IdWithMotif).OrderBy(c => c))));
//				returnString = GlobalVariables.CheckLengthOfOutput(returnString);
//				return (returnString, nullptr);
//			}
//			else
//			{
//				return (string.Join(" ", list[0].Values.Select(c => c.IdWithMotif).OrderBy(c => c)), firstDict);
//			}
//		}

//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//		private static(string ResolvedString, Nullable<double> ResolvedValue) ResolveF2(IEnumerable<double> enumerable)
//		{
//			var list = enumerable.ToList();
//			if (list.Max() - list.Min() < ToleranceForDoubleResolutionF2)
//			{
//				return (list.Average().ToString("F2", CultureInfo.InvariantCulture), list.Average());
//			}
//			else
//			{
//				var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F2", CultureInfo.InvariantCulture))));
//				return (returnString, nullptr);
//			}
//		}

//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//		private static(string ResolvedString, Nullable<double> ResolvedValue) Resolve(IEnumerable<double> enumerable)
//		{
//			var list = enumerable.ToList();
//			if (list.Max() - list.Min() < ToleranceForDoubleResolutionF5)
//			{
//				return (list.Average().ToString("F5", CultureInfo.InvariantCulture), list.Average());
//			}
//			else
//			{
//				var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))));
//				return (returnString, nullptr);
//			}
//		}

//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//		private static(string ResolvedString, Nullable<int> ResolvedValue) Resolve(IEnumerable<int> enumerable)
//		{
//			var list = enumerable.ToList();
//			var first = list[0];
//			if (list.All(b => first.Equals(b)))
//			{
//				return (first.ToString(CultureInfo.InvariantCulture), first);
//			}
//			else
//			{
//				var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString(CultureInfo.InvariantCulture))));
//				return (returnString, nullptr);
//			}
//		}

//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//		private static(string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable)
//		{
//			var list = enumerable.ToList();
//			string first = list.FirstOrDefault(b => b != nullptr);
//			// Only first if list is either all null or all equal to the first
//			if (list.All(b => b == nullptr) || list.All(b => first.Equals(b)))
//			{
//				return (first, first);
//			}
//			else
//			{
//				var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list));
//				return (returnString, nullptr);
//			}
//		}

//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//		private static(string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable, string ambiguousIfNull)
//		{
//			var list = enumerable.ToList();
//			string first = list.FirstOrDefault(b => b != nullptr);
//			// Only first if list is either all null or all equal to the first
//			if (list.All(b => b == nullptr) || list.All(b => first.Equals(b)))
//			{
//				return (first, first);
//			}
//			// use only distinct names if all of the base sequences are the same
//			else if (ambiguousIfNull != nullptr)
//			{
//				var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Distinct()));
//				return (returnString, nullptr);
//			}
//			else
//			{
//				var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list));
//				return (returnString, nullptr);
//			}
//		}

//C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
		private static void AddBasicMatchData(Dictionary<string, string> s, PeptideSpectralMatch psm)

//C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
		private static void AddPeptideSequenceData(Dictionary<string, string> s, PeptideSpectralMatch psm, IReadOnlyDictionary<string, int> ModsToWritePruned)

		/// <summary>
		/// Determines whether a peptide includes a splice site
		/// </summary>
		/// <param name="pep"></param>
		/// <param name="site"></param>
		/// <returns></returns>
//C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
		private static bool Includes(PeptideWithSetModifications pep, SpliceSite site)

		/// <summary>
		/// Checks for an intersection between a peptide and applied variant that shows a sequence change.
		/// </summary>
		/// <param name="pep"></param>
		/// <param name="appliedVariation"></param>
		/// <returns></returns>
//C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
		private static bool IntersectsWithVariation(PeptideWithSetModifications pep, SequenceVariation appliedVariation, bool checkUnique)

		/// <summary>
		/// Makes the string representing a detected sequence variation, including any modifications on a variant amino acid
		/// </summary>
		/// <param name="p"></param>
		/// <param name="d"></param>
		/// <returns></returns>
//C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
		private static string SequenceVariantString(PeptideWithSetModifications p, SequenceVariation applied)

//C# TO C++ CONVERTER TODO TASK: The following line could not be converted:
		public static void AddMatchedIonsData(Dictionary<string, string> s, PeptideSpectralMatch psm)
