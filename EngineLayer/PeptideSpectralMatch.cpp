#include "PeptideSpectralMatch.h"
#include "IScan.h"

using namespace Chemistry;
using namespace EngineLayer::FdrAnalysis;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

	PeptideSpectralMatch::PeptideSpectralMatch(PeptideWithSetModifications *peptide, int notch, double score, int scanIndex, IScan *scan, DigestionParams *digestionParams, std::vector<MatchedFragmentIon*> &matchedFragmentIons) : DigestionParams(digestionParams)
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

            AddOrReplace(peptide, score, notch, true, matchedFragmentIons);
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
#endif
        
	std::vector<double> PeptideSpectralMatch::getFeatures() const
	{
            return std::vector<std::any> {BankersRounding::round(getScore()), getScore() - BankersRounding::round(getScore())};
	}

	std::string PeptideSpectralMatch::GetTabSeparatedHeader()
	{
            return std::string::Join("\t", DataDictionary(nullptr, nullptr).Keys);
	}

	void PeptideSpectralMatch::AddOrReplace(PeptideWithSetModifications *pwsm, double newScore, int notch,
                                                bool reportAllAmbiguity, std::vector<MatchedFragmentIon*> &matchedFragmentIons)
	{
            if (newScore - getScore() > ToleranceForScoreDifferentiation) //if new score beat the old score, overwrite it
            {
                _bestMatchingPeptides.clear();
                _bestMatchingPeptides.push_back((notch, pwsm));
                
                if (getScore() - getRunnerUpScore() > ToleranceForScoreDifferentiation)
                {
                    setRunnerUpScore(getScore());
                }
                
                setScore(newScore);
                
                getPeptidesToMatchingFragments().clear();
                getPeptidesToMatchingFragments().emplace(pwsm, matchedFragmentIons);
            }
            else if (newScore - getScore() > -ToleranceForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
            {
                _bestMatchingPeptides.push_back((notch, pwsm));
                
                if (getPeptidesToMatchingFragments().find(pwsm) == getPeptidesToMatchingFragments().end())
                {
                    getPeptidesToMatchingFragments().emplace(pwsm, matchedFragmentIons);
                }
            }
            else if (getScore() - getRunnerUpScore() > ToleranceForScoreDifferentiation)
            {
                setRunnerUpScore(newScore);
            }
	}
        
	std::string PeptideSpectralMatch::ToString()
	{
            return ToString(std::unordered_map<std::string, int>());
	}
        
	std::string PeptideSpectralMatch::ToString(IReadOnlyDictionary<std::string, int> *ModstoWritePruned)
	{
            return std::string::Join("\t", DataDictionary(this, ModstoWritePruned).Values);
	}
        
	std::unordered_map<std::string, std::string> PeptideSpectralMatch::DataDictionary(PeptideSpectralMatch *psm,
                                                                                          std::unordered_map<std::string, int> *ModsToWritePruned)
	{
            std::unordered_map<std::string, std::string> s;
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
            FdrInfo tempVar();
            setFdrInfo(&tempVar);
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
            setIsDecoy(_bestMatchingPeptides.Any([&] (std::any p)  {
                        p::Pwsm::Protein::IsDecoy;
                    }));

            setIsContaminant(_bestMatchingPeptides.Any([&] (std::any p) {
			p::Pwsm::Protein::IsContaminant;
                    }));
            
            setFullSequence(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                            b::Pwsm::FullSequence;
                        }))->ResolvedValue);

            setBaseSequence(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                            b::Pwsm::BaseSequence;
                        }))->ResolvedValue);

            setPeptideLength(Resolve(_bestMatchingPeptides.Select([&] (std::any b){
                            b::Pwsm->Length;
                        }))->ResolvedValue);

            setOneBasedStartResidueInProtein(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                            b::Pwsm::OneBasedStartResidueInProtein;
                        }))->ResolvedValue);

            setOneBasedEndResidueInProtein(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                            b::Pwsm::OneBasedEndResidueInProtein;
                        }))->ResolvedValue);
            
            setProteinLength(Resolve(_bestMatchingPeptides.Select([&] (std::any b) {
                            b::Pwsm::Protein->Length;
                        }))->ResolvedValue);

            setPeptideMonisotopicMass(Resolve(_bestMatchingPeptides.Select([&] (std::any b)  {
                            b::Pwsm::MonoisotopicMass;
                        }))->ResolvedValue);

            setProteinAccession(Resolve(_bestMatchingPeptides.Select([&] (std::any b)	{
                            b::Pwsm::Protein::Accession;
                        }))->ResolvedValue);

            setOrganism(Resolve(_bestMatchingPeptides.Select([&] (std::any b)	{
                            b::Pwsm::Protein::Organism;
                        }))->ResolvedValue);

            setModsIdentified(Resolve(_bestMatchingPeptides.Select([&] (std::any b)	{
                            b::Pwsm::AllModsOneIsNterminus;
                        }))->ResolvedValue);

            setModsChemicalFormula(Resolve(_bestMatchingPeptides.Select([&] (std::any b){
                            b::Pwsm::AllModsOneIsNterminus->Select([&] (std::any c)  {
                                    (c->Value);
                                });
                        }))->ResolvedValue);

            setNotch(Resolve(_bestMatchingPeptides.Select([&] (std::any b)  {
                            b::Notch;
                        }))->ResolvedValue);

            // if the PSM matches a target and a decoy and they are the SAME SEQUENCE, remove the decoy
            if (getIsDecoy())
            {
                bool removedPeptides = false;
                auto hits = _bestMatchingPeptides.GroupBy([&] (std::any p)  {
                        p::Pwsm::FullSequence;
                    });
                
                for (auto hit : hits)
                {
                    if (hit->Any([&] (std::any p) {
                                p::Pwsm::Protein::IsDecoy;
                            }) && hit->Any([&] (std::any p) {
                                    !p::Pwsm::Protein::IsDecoy;
				})){_bestMatchingPeptides.RemoveAll([&] (std::any p)   {
                                return p::Pwsm->FullSequence == hit->Key && p::Pwsm::Protein::IsDecoy;
                            });
                        {
                            removedPeptides = true;
                        }
                    }
		}

                if (removedPeptides)
                {
                    ResolveAllAmbiguities();
                }
            }
        }
        
	void PeptideSpectralMatch::TrimProteinMatches(std::unordered_set<Protein*> &parsimoniousProteins)
	{
            if (IsDecoy)
            {
                if (_bestMatchingPeptides::Any([&] (std::any p) {
                     return std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(), p::Pwsm::Protein) != parsimoniousProteins.end() &&
                         p::Pwsm::Protein::IsDecoy;
			}))
                {
                    _bestMatchingPeptides::RemoveAll([&] (std::any p) {
                            std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(), p::Item2->Protein) == parsimoniousProteins.end();
                        });
                }
                // else do nothing
            }
            else
            {
                _bestMatchingPeptides::RemoveAll([&] (std::any p){
                        std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(), p::Item2->Protein) == parsimoniousProteins.end();
                    });
            }
            
            ResolveAllAmbiguities();
	}

	public void PeptideSpectralMatch::AddProteinMatch((int, PeptideWithSetModifications)  *peptideWithNotch)
	{
            _bestMatchingPeptides->Add(peptideWithNotch);
            ResolveAllAmbiguities();
	}

	private void PeptideSpectralMatch::AddBasicMatchData(std::unordered_map<std::string, std::string> &s, PeptideSpectralMatch *psm)
	{
            s["File Name"] = psm == nullptr ? " " : Path::GetFileNameWithoutExtension(psm->getFullFilePath());
            s["Scan Number"] = psm == nullptr ? " " : psm->getScanNumber().ToString(CultureInfo::InvariantCulture);
            s["Scan Retention Time"] = psm == nullptr ? " " : psm->getScanRetentionTime().ToString("F5", CultureInfo::InvariantCulture);
            s["Num Experimental Peaks"] = psm == nullptr ? " " : psm->getScanExperimentalPeaks().ToString("F5", CultureInfo::InvariantCulture);
            s["Total Ion Current"] = psm == nullptr ? " " : psm->getTotalIonCurrent().ToString("F5", CultureInfo::InvariantCulture);
            s["Precursor Scan Number"] = psm == nullptr ? " " : psm->getPrecursorScanNumber().HasValue ? psm->getPrecursorScanNumber().Value.ToString(CultureInfo::InvariantCulture) : "unknown";
            s["Precursor Charge"] = psm == nullptr ? " " : psm->getScanPrecursorCharge().ToString("F5", CultureInfo::InvariantCulture);
            s["Precursor MZ"] = psm == nullptr ? " " : psm->getScanPrecursorMonoisotopicPeakMz().ToString("F5", CultureInfo::InvariantCulture);
            s["Precursor Mass"] = psm == nullptr ? " " : psm->getScanPrecursorMass().ToString("F5", CultureInfo::InvariantCulture);
            s["Score"] = psm == nullptr ? " " : psm->getScore().ToString("F3", CultureInfo::InvariantCulture);
            s["Delta Score"] = psm == nullptr ? " " : psm->getDeltaScore().ToString("F3", CultureInfo::InvariantCulture);
            s["Notch"] = psm == nullptr ? " " : Resolve(psm->BestMatchingPeptides->Select([&] (std::any p)
            {
                p::Notch;
            }))->ResolvedString;
            
            s["Different Peak Matches"] = psm == nullptr ? " " : psm->getNumDifferentMatchingPeptides().ToString("F5", CultureInfo::InvariantCulture);
	}

	void PeptideSpectralMatch::AddPeptideSequenceData(std::unordered_map<std::string, std::string> s, PeptideSpectralMatch *psm,
                                                          std::unordered_map<std::string, int> *ModsToWritePruned)
	{
            bool pepWithModsIsNull = psm == nullptr || psm->BestMatchingPeptides == nullptr || !psm->BestMatchingPeptides.Any();

            std::vector<PeptideWithSetModifications*> pepsWithMods = pepWithModsIsNull ? nullptr : psm->BestMatchingPeptides->Select([&] (std::any p)
		{
                    p::Peptide;
		}).ToList();

            s["Base Sequence"] = pepWithModsIsNull ? " " : Resolve(pepWithModsIsNull ? nullptr : pepsWithMods.Select([&] (std::any b)
		{
                    b::BaseSequence;
		}))->ResolvedString;
            s["Full Sequence"] = pepWithModsIsNull ? " " : Resolve(pepWithModsIsNull ? nullptr : pepsWithMods.Select([&] (std::any b)
		{
                    b::FullSequence;
		}))->ResolvedString;
            s["Essential Sequence"] = pepWithModsIsNull ? " " : Resolve(pepWithModsIsNull ? nullptr : pepsWithMods.Select([&] (std::any b)
		{
                    b::EssentialSequence(ModsToWritePruned);
		}))->ResolvedString;
            s["Mods"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::AllModsOneIsNterminus;
		}))->ResolvedString;
            s["Mods Chemical Formulas"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any p)
		{
                    p::AllModsOneIsNterminus->Select([&] (std::any v)
			{
                            v->Value;
			});
		}))->ResolvedString;
            s["Mods Combined Chemical Formula"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::AllModsOneIsNterminus->Select([&] (std::any c)
			{
                            (dynamic_cast<Modification*>(c->Value));
			});
		}))->ResolvedString;
            s["Num Variable Mods"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::NumVariableMods;
		}))->Item1;
            s["Missed Cleavages"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::MissedCleavages.ToString(CultureInfo::InvariantCulture);
		}))->ResolvedString;
            s["Peptide Monoisotopic Mass"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::MonoisotopicMass;
		}))->ResolvedString;
            s["Mass Diff (Da)"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    return psm->getScanPrecursorMass() - b::MonoisotopicMass;
		}))->ResolvedString;
            s["Mass Diff (ppm)"] = pepWithModsIsNull ? " " : ResolveF2(pepsWithMods.Select([&] (std::any b)
		{
                    ((psm->getScanPrecursorMass() - b::MonoisotopicMass) / b::MonoisotopicMass * 1e6);
		})).ResolvedString;
            s["Protein Accession"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::Protein::Accession;
		}), psm->getFullSequence())->ResolvedString;
            s["Protein Name"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::Protein->FullName;
		}), psm->getFullSequence())->ResolvedString;
            s["Gene Name"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    std::string::Join(", ", b::Protein::GeneNames->Select([&] (std::any d)
			{
                            StringHelper::formatSimple("{0}:{1}", d::Item1, d::Item2);
			}));
		}), psm->getFullSequence())->ResolvedString;
            s["Intersecting Sequence Variations"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    std::string::Join(", ", b::Protein::AppliedSequenceVariations::Where([&] (std::any av)
			{
                            IntersectsWithVariation(b, av, false);
			})->Select([&] (std::any av)  {
                                SequenceVariantString(b, av);
                            }));
		}))->ResolvedString;
            s["Identified Sequence Variations"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    std::string::Join(", ", b::Protein::AppliedSequenceVariations::Where([&] (std::any av)
			{
                            IntersectsWithVariation(b, av, true);
			})->Select([&] (std::any av)
			{
                            SequenceVariantString(b, av);
			}));
		}))->ResolvedString;
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
            s["Organism Name"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::Protein::Organism;
		}))->Item1;
            s["Contaminant"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::Protein::IsContaminant ? "Y" : "N";
		}))->Item1;
            s["Decoy"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::Protein::IsDecoy ? "Y" : "N";
		}))->Item1;
            s["Peptide Description"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::PeptideDescription;
		}))->Item1;
            s["Start and End Residues In Protein"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    (StringHelper::formatSimple("[{0} to {1}]", b::OneBasedStartResidueInProtein.ToString(CultureInfo::InvariantCulture),
                                                b::OneBasedEndResidueInProtein.ToString(CultureInfo::InvariantCulture)));
		}), psm->getFullSequence())->ResolvedString;
            s["Previous Amino Acid"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::PreviousAminoAcid.ToString();
		}))->ResolvedString;
            s["Next Amino Acid"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
                    b::NextAminoAcid.ToString();
		}))->ResolvedString;
            
            std::string allScores = " ";
            std::string theoreticalsSearched = " ";
            if (!pepWithModsIsNull && psm->getFdrInfo() != nullptr && psm->getFdrInfo()->getCalculateEValue())
            {
                allScores = std::string::Join(";", psm->getAllScores().Select([&] (std::any p)
                {
                    p.ToString("F2", CultureInfo::InvariantCulture);
                }));
                theoreticalsSearched = std::to_string(psm->getAllScores().size());
            }

            s["All Scores"] = allScores;
            s["Theoreticals Searched"] = theoreticalsSearched;
            s["Decoy/Contaminant/Target"] = pepWithModsIsNull ? " " : psm->getIsDecoy() ? "D" : psm->getIsContaminant() ? "C" : "T";
	}

	bool PeptideSpectralMatch::Includes(PeptideWithSetModifications *pep, SpliceSite *site)
	{
            return pep->OneBasedStartResidueInProtein <= site->OneBasedBeginPosition &&
                pep->OneBasedEndResidueInProtein >= site->OneBasedEndPosition;
	}

	bool PeptideSpectralMatch::IntersectsWithVariation(PeptideWithSetModifications *pep, SequenceVariation *appliedVariation, bool checkUnique)
	{
            // does it intersect? 
            int intersectOneBasedStart = std::max(pep->OneBasedStartResidueInProtein, appliedVariation->OneBasedBeginPosition);
            int intersectOneBasedEnd = std::min(pep->OneBasedEndResidueInProtein, appliedVariation->OneBasedEndPosition);
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
                int variantZeroBasedStart = intersectOneBasedStart - appliedVariation->OneBasedBeginPosition;
                bool origSeqIsShort = appliedVariation->OriginalSequence->Length - variantZeroBasedStart < intersectSize;
                bool origSeqIsLong = appliedVariation->OriginalSequence->Length > intersectSize && pep->OneBasedEndResidueInProtein > intersectOneBasedEnd;
                if (origSeqIsShort || origSeqIsLong)
                {
                    return true;
                }
                
                // is the variant sequence intersecting the peptide different than the original sequence?
                std::string originalAtIntersect = appliedVariation->OriginalSequence->substr(intersectOneBasedStart - appliedVariation->OneBasedBeginPosition, intersectSize);
                std::string variantAtIntersect = appliedVariation->VariantSequence->substr(intersectOneBasedStart - appliedVariation->OneBasedBeginPosition, intersectSize);
                return originalAtIntersect != variantAtIntersect;
            }
	}

	std::string PeptideSpectralMatch::SequenceVariantString(PeptideWithSetModifications *p, SequenceVariation *applied)
	{
            auto modsOnVariantOneIsNTerm = p->AllModsOneIsNterminus.Where([&] (std::any kv)
		{
                    return kv->Key == 1 && applied->OneBasedBeginPosition == 1 || applied->OneBasedBeginPosition <= kv::Key - 2 + p->OneBasedStartResidueInProtein && kv::Key - 2 + p->OneBasedStartResidueInProtein <= applied->OneBasedEndPosition;
		}).ToDictionary([&] (std::any kv)
		{
                    return kv::Key - applied->OneBasedBeginPosition + 1;
		}, [&] (std::any kv)
		{
                    kv->Value;
		});
            PeptideWithSetModifications *variantWithAnyMods = new PeptideWithSetModifications(p->Protein, p->DigestionParams,
                                                                                              applied->OneBasedBeginPosition,
                                                                                              applied->OneBasedEndPosition,
                                                                                              p->CleavageSpecificityForFdrCategory,
                                                                                              p->PeptideDescription,
                                                                                              p->MissedCleavages,
                                                                                              modsOnVariantOneIsNTerm,
                                                                                              p->NumFixedMods);

            delete variantWithAnyMods;
            return StringHelper::formatSimple("{0}{1}{2}", applied->OriginalSequence, applied->OneBasedBeginPosition,
                                              variantWithAnyMods->FullSequence);
	}

	void PeptideSpectralMatch::AddMatchedIonsData(std::unordered_map<std::string, std::string> &s, PeptideSpectralMatch *psm)
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
                auto matchedIons = psm->MatchedFragmentIons;
                if (matchedIons == nullptr)
                {
                    matchedIons = psm->getPeptidesToMatchingFragments().First()->Value;
                }
                
                // using ", " instead of "," improves human readability
                const std::string delimiter = ", ";
                
                auto matchedIonsGroupedByProductType = matchedIons->GroupBy([&] (std::any i) {
                        i::NeutralTheoreticalProduct::ProductType;
                    }).OrderBy([&] (std::any i)  {
                            i::Key;
			}).ToList();

			for (auto productType : matchedIonsGroupedByProductType)
			{
                            auto products = productType.OrderBy([&] (std::any p)
				{
                                    p::NeutralTheoreticalProduct::TerminusFragment::FragmentNumber;
				}).ToList();

                            std::for_each(stringBuilders.begin(), stringBuilders.end(), [&] (std::any p)
				{
                                    p->Append("[");
				});

                            for (int i = 0; i < products.size(); i++)
                            {
                                MatchedFragmentIon *ion = products[i];
                                std::string ionLabel;
                                
                                double massError = ion->Mz.ToMass(ion->Charge) - ion->NeutralTheoreticalProduct.NeutralMass;
                                double ppmMassError = massError / ion->NeutralTheoreticalProduct.NeutralMass * 1e6;
                                
                                if (ion->NeutralTheoreticalProduct->NeutralLoss == 0)
                                {
                                    // no neutral loss
                                    ionLabel = ion->NeutralTheoreticalProduct.ProductType + "" + ion->NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + "+" + ion->Charge;
                                }
                                else
                                {
                                    // ion label with neutral loss
                                    
                                    ionLabel = "(" + ion->NeutralTheoreticalProduct.ProductType + "" + ion->NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + "-" + ion->NeutralTheoreticalProduct.NeutralLoss.ToString("F2") + ")" + "+" + ion->Charge;
                                }
                                
                                // append ion label
                                seriesStringBuilder->append(ionLabel);
                                
                                // append experimental m/z
                                
                                mzStringBuilder->append(ionLabel + ":" + ion->Mz.ToString("F5"));
                                
                                // append absolute mass error
                                
                                fragmentDaErrorStringBuilder->append(ionLabel + ":" + massError.ToString("F5"));
                                
                                // append ppm mass error
                                
                                fragmentPpmErrorStringBuilder->append(ionLabel + ":" + ppmMassError.ToString("F2"));
                                
                                // append fragment ion intensity
                                
                                fragmentIntensityStringBuilder->append(ionLabel + ":" + ion->Intensity.ToString("F0"));
                                
                                // append delimiter ", "
                                if (i < products.size() - 1)
                                {
                                    std::for_each(stringBuilders.begin(), stringBuilders.end(), [&] (std::any p)  {
                                            p->Append(delimiter);
                                        });
                                }
                            }
                            
                            // append product type delimiter
                            
                            >
                            {
                                p->Append("];
				}"));
			}}s["Matched Ion Series"] = nullPsm ? " " : seriesStringBuilder->ToString()->TrimEnd(L';');
            //C# TO C++ CONVERTER TODO TASK: The following lambda expression could not be converted:
            stringBuilders.ForEach(p => TangibleLambdaToken86}s["Matched Ion Series"];
            s["Matched Ion Mass-To-Charge Ratios"] = nullPsm ? " " : mzStringBuilder->toString()->TrimEnd(L';');
            s["Matched Ion Mass Diff (Da)"] = nullPsm ? " " : fragmentDaErrorStringBuilder->toString()->TrimEnd(L';');
            s["Matched Ion Mass Diff (Ppm)"] = nullPsm ? " " : fragmentPpmErrorStringBuilder->toString()->TrimEnd(L';');
            s["Matched Ion Intensities"] = nullPsm ? " " : fragmentIntensityStringBuilder->toString()->TrimEnd(L';');
            
            // number of matched ions
            s["Matched Ion Counts"] = nullPsm ? " " : std::to_string(psm->MatchedFragmentIons->Count);
        }
        
        //C# TO C++ CONVERTER TODO TASK: Local functions are not converted by C# to C++ Converter:
        static void AddMatchScoreData(std::unordered_map<std::string, std::string> s, PeptideSpectralMatch peptide)
        {
	//		{
	//			string localizedScores = " ";
	//			string improvementPossible = " ";
	//			if (peptide != nullptr && peptide.LocalizedScores != nullptr)
	//			{
	//				localizedScores = GlobalVariables.CheckLengthOfOutput(("[" + string.Join(",", peptide.LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]"));
	//				improvementPossible = (peptide.LocalizedScores.Max() - peptide.Score).ToString("F3", CultureInfo.InvariantCulture);
	//			}
	//			s["Localized Scores"] = localizedScores;
	//			s["Improvement Possible"] = improvementPossible;
	//
	//			string cumulativeTarget = " ";
	//			string cumulativeDecoy = " ";
	//			string qValue = " ";
	//			string cumulativeTargetNotch = " ";
	//			string cumulativeDecoyNotch = " ";
	//			string qValueNotch = " ";
	//			string eValue = " ";
	//			string eScore = " ";
	//			if (peptide != nullptr && peptide.FdrInfo != nullptr)
	//			{
	//				cumulativeTarget = peptide.FdrInfo.CumulativeTarget.ToString(CultureInfo.InvariantCulture);
	//				cumulativeDecoy = peptide.FdrInfo.CumulativeDecoy.ToString(CultureInfo.InvariantCulture);
	//				qValue = peptide.FdrInfo.QValue.ToString("F6", CultureInfo.InvariantCulture);
	//				cumulativeTargetNotch = peptide.FdrInfo.CumulativeTargetNotch.ToString(CultureInfo.InvariantCulture);
	//				cumulativeDecoyNotch = peptide.FdrInfo.CumulativeDecoyNotch.ToString(CultureInfo.InvariantCulture);
	//				qValueNotch = peptide.FdrInfo.QValueNotch.ToString("F6", CultureInfo.InvariantCulture);
	//				if (peptide.FdrInfo.CalculateEValue)
	//				{
	//					eValue = peptide.FdrInfo.EValue.ToString("F6", CultureInfo.InvariantCulture);
	//					eScore = peptide.FdrInfo.EScore.ToString("F6", CultureInfo.InvariantCulture);
	//				}
	//			}
	//			s["Cumulative Target"] = cumulativeTarget;
	//			s["Cumulative Decoy"] = cumulativeDecoy;
	//			s["QValue"] = qValue;
	//			s["Cumulative Target Notch"] = cumulativeTargetNotch;
	//			s["Cumulative Decoy Notch"] = cumulativeDecoyNotch;
	//			s["QValue Notch"] = qValueNotch;
	//			s["eValue"] = eValue;
	//			s["eScore"] = eScore;
	//		}
	//}

		delete fragmentIntensityStringBuilder;
		delete fragmentPpmErrorStringBuilder;
		delete fragmentDaErrorStringBuilder;
		delete mzStringBuilder;
		delete seriesStringBuilder;
	}

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, ChemicalFormula ResolvedValue) Resolve(IEnumerable<IEnumerable<Modification>> enumerable);
        static std::tuple<std::string, ChemicalFormula*> Resolve(std::vector<std::vector<Modification *>> enumerable)
        {
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
        }


        static std::tuple<std::string, std::unordered_map<std::string, int>> Resolve(std::vector<std::unoredered_map<int, Modification*>> enumerable)
        {
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
        }

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, Nullable<double> ResolvedValue) ResolveF2(IEnumerable<double> enumerable);
        static std::tuple<std::string, std::optional<double>> ResolveF2(std::vector<double> enumerable)
        {
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
        }

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, Nullable<double> ResolvedValue) Resolve(IEnumerable<double> enumerable);
        static std::tuple<std::string, std::optional<double>> Resolve(std::vector<double> enumerable)
        {
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
        }

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, Nullable<int> ResolvedValue) Resolve(IEnumerable<int> enumerable);
        static std::tuple<std::string, std::optional<int>>  Resolve(IEnumerable<int> enumerable)
        {
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
        }

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable);
        static std::tuple<std::string, std::string> Resolve(std::vector<std::string> enumerable)
        {
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
        }

        //C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
        //static(string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable, string ambiguousIfNull);
        static std::tuple<std::string, std::string> Resolve(std::vector<std::string> enumerable, std::string ambiguousIfNull)
        {
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
        }
    }
}

        
