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
		_bestMatchingPeptides = std::vector<(int, PeptideWithSetModifications)*>();
		ScanIndex = scanIndex;
		FullFilePath = scan->getFullFilePath();
		ScanNumber = scan->getOneBasedScanNumber();
		PrecursorScanNumber = scan->getOneBasedPrecursorScanNumber();
		ScanRetentionTime = scan->getRetentionTime();
		ScanExperimentalPeaks = scan->getNumPeaks();
		TotalIonCurrent = scan->getTotalIonCurrent();
		ScanPrecursorCharge = scan->getPrecursorCharge();
		ScanPrecursorMonoisotopicPeakMz = scan->getPrecursorMonoisotopicPeakMz();
		ScanPrecursorMass = scan->getPrecursorMass();
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

	std::wstring PeptideSpectralMatch::getFullSequence() const
	{
		return privateFullSequence;
	}

	void PeptideSpectralMatch::setFullSequence(const std::wstring &value)
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

	std::wstring PeptideSpectralMatch::getBaseSequence() const
	{
		return privateBaseSequence;
	}

	void PeptideSpectralMatch::setBaseSequence(const std::wstring &value)
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

	std::wstring PeptideSpectralMatch::getProteinAccession() const
	{
		return privateProteinAccession;
	}

	void PeptideSpectralMatch::setProteinAccession(const std::wstring &value)
	{
		privateProteinAccession = value;
	}

	std::wstring PeptideSpectralMatch::getOrganism() const
	{
		return privateOrganism;
	}

	void PeptideSpectralMatch::setOrganism(const std::wstring &value)
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

	std::unordered_map<std::wstring, int> PeptideSpectralMatch::getModsIdentified() const
	{
		return privateModsIdentified;
	}

	void PeptideSpectralMatch::setModsIdentified(const std::unordered_map<std::wstring, int> &value)
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

	std::wstring PeptideSpectralMatch::getFullFilePath() const
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

	private *IEnumerable < PeptideSpectralMatch::(int Notch, PeptideWithSetModifications *Peptide)
	{
		get;

	get
	{
			return _bestMatchingPeptides.OrderBy([&] (std::any p)
			{
				p::Item2->FullSequence;
			}).ThenBy([&] (std::any p)
			{
				p::Item2->Protein.Accession;
			}).ThenBy([&] (std::any p)
			{
				p::Item2->OneBasedStartResidueInProtein;
			});
	}

	std::vector<double> PeptideSpectralMatch::getFeatures() const
	{
		return std::vector<std::any> {BankersRounding::round(getScore()), getScore() - BankersRounding::round(getScore())};
	}

	std::wstring PeptideSpectralMatch::GetTabSeparatedHeader()
	{
		return std::wstring::Join(L"\t", DataDictionary(nullptr, nullptr).Keys);
	}

	void PeptideSpectralMatch::AddOrReplace(PeptideWithSetModifications *pwsm, double newScore, int notch, bool reportAllAmbiguity, std::vector<MatchedFragmentIon*> &matchedFragmentIons)
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

	std::wstring PeptideSpectralMatch::ToString()
	{
		return ToString(std::unordered_map<std::wstring, int>());
	}

	std::wstring PeptideSpectralMatch::ToString(IReadOnlyDictionary<std::wstring, int> *ModstoWritePruned)
	{
		return std::wstring::Join(L"\t", DataDictionary(this, ModstoWritePruned).Values);
	}

	std::unordered_map<std::wstring, std::wstring> PeptideSpectralMatch::DataDictionary(PeptideSpectralMatch *psm, IReadOnlyDictionary<std::wstring, int> *ModsToWritePruned)
	{
		std::unordered_map<std::wstring, std::wstring> s;
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

	void PeptideSpectralMatch::SetFdrValues(double cumulativeTarget, double cumulativeDecoy, double qValue, double cumulativeTargetNotch, double cumulativeDecoyNotch, double qValueNotch, double maximumLikelihood, double eValue, double eScore, bool calculateEValue)
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
		setIsDecoy(_bestMatchingPeptides.Any([&] (std::any p)
		{
			p::Pwsm::Protein::IsDecoy;
		}));
		setIsContaminant(_bestMatchingPeptides.Any([&] (std::any p)
		{
			p::Pwsm::Protein::IsContaminant;
		}));

		setFullSequence(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::FullSequence;
		}))->ResolvedValue);
		setBaseSequence(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::BaseSequence;
		}))->ResolvedValue);
		setPeptideLength(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm->Length;
		}))->ResolvedValue);
		setOneBasedStartResidueInProtein(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::OneBasedStartResidueInProtein;
		}))->ResolvedValue);
		setOneBasedEndResidueInProtein(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::OneBasedEndResidueInProtein;
		}))->ResolvedValue);
		setProteinLength(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::Protein->Length;
		}))->ResolvedValue);
		setPeptideMonisotopicMass(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::MonoisotopicMass;
		}))->ResolvedValue);
		setProteinAccession(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::Protein::Accession;
		}))->ResolvedValue);
		setOrganism(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::Protein::Organism;
		}))->ResolvedValue);
		setModsIdentified(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::AllModsOneIsNterminus;
		}))->ResolvedValue);
		setModsChemicalFormula(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Pwsm::AllModsOneIsNterminus->Select([&] (std::any c)
			{
				(c->Value);
			});
		}))->ResolvedValue);
		setNotch(Resolve(_bestMatchingPeptides.Select([&] (std::any b)
		{
			b::Notch;
		}))->ResolvedValue);

		// if the PSM matches a target and a decoy and they are the SAME SEQUENCE, remove the decoy
		if (getIsDecoy())
		{
			bool removedPeptides = false;
			auto hits = _bestMatchingPeptides.GroupBy([&] (std::any p)
			{
				p::Pwsm::FullSequence;
			});

			for (auto hit : hits)
			{
				if (hit->Any([&] (std::any p)
				{
					p::Pwsm::Protein::IsDecoy;
				}) && hit->Any([&] (std::any p)
				{
					!p::Pwsm::Protein::IsDecoy;
				})){_bestMatchingPeptides.RemoveAll([&] (std::any p)
				{
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

	public void <missing_class_definition>::TrimProteinMatches(std::unordered_set<Protein*> &parsimoniousProteins)
	{
		if (IsDecoy)
		{
			if (_bestMatchingPeptides::Any([&] (std::any p)
			{
				return std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(), p::Pwsm::Protein) != parsimoniousProteins.end() && p::Pwsm::Protein::IsDecoy;
			}))
			{
				_bestMatchingPeptides::RemoveAll([&] (std::any p)
				{
					!std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(), p::Item2->Protein) != parsimoniousProteins.end();
				});
			}
			// else do nothing
		}
		else
		{
			_bestMatchingPeptides::RemoveAll([&] (std::any p)
			{
				!std::find(parsimoniousProteins.begin(), parsimoniousProteins.end(), p::Item2->Protein) != parsimoniousProteins.end();
			});
		}

		ResolveAllAmbiguities();
	}

	public void <missing_class_definition>::AddProteinMatch((int, PeptideWithSetModifications)  *peptideWithNotch)
	{
		_bestMatchingPeptides->Add(peptideWithNotch);
		ResolveAllAmbiguities();
	}

	private void <missing_class_definition>::AddBasicMatchData(std::unordered_map<std::wstring, std::wstring> &s, PeptideSpectralMatch *psm)
	{
		s[L"File Name"] = psm == nullptr ? L" " : Path::GetFileNameWithoutExtension(psm->getFullFilePath());
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Scan Number"] = psm == nullptr ? L" " : psm->getScanNumber().ToString(CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Scan Retention Time"] = psm == nullptr ? L" " : psm->getScanRetentionTime().ToString(L"F5", CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Num Experimental Peaks"] = psm == nullptr ? L" " : psm->getScanExperimentalPeaks().ToString(L"F5", CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Total Ion Current"] = psm == nullptr ? L" " : psm->getTotalIonCurrent().ToString(L"F5", CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Precursor Scan Number"] = psm == nullptr ? L" " : psm->getPrecursorScanNumber().HasValue ? psm->getPrecursorScanNumber().Value.ToString(CultureInfo::InvariantCulture) : L"unknown";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Precursor Charge"] = psm == nullptr ? L" " : psm->getScanPrecursorCharge().ToString(L"F5", CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Precursor MZ"] = psm == nullptr ? L" " : psm->getScanPrecursorMonoisotopicPeakMz().ToString(L"F5", CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Precursor Mass"] = psm == nullptr ? L" " : psm->getScanPrecursorMass().ToString(L"F5", CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Score"] = psm == nullptr ? L" " : psm->getScore().ToString(L"F3", CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Delta Score"] = psm == nullptr ? L" " : psm->getDeltaScore().ToString(L"F3", CultureInfo::InvariantCulture);
		s[L"Notch"] = psm == nullptr ? L" " : Resolve(psm->BestMatchingPeptides->Select([&] (std::any p)
		{
			p::Notch;
		}))->ResolvedString;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		s[L"Different Peak Matches"] = psm == nullptr ? L" " : psm->getNumDifferentMatchingPeptides().ToString(L"F5", CultureInfo::InvariantCulture);
	}

	private void <missing_class_definition>::AddPeptideSequenceData(std::unordered_map<std::wstring, std::wstring> &s, PeptideSpectralMatch *psm, IReadOnlyDictionary<std::wstring, int> *ModsToWritePruned)
	{
		bool pepWithModsIsNull = psm == nullptr || psm->BestMatchingPeptides == nullptr || !psm->BestMatchingPeptides.Any();

		std::vector<PeptideWithSetModifications*> pepsWithMods = pepWithModsIsNull ? nullptr : psm->BestMatchingPeptides->Select([&] (std::any p)
		{
			p::Peptide;
		}).ToList();

		s[L"Base Sequence"] = pepWithModsIsNull ? L" " : Resolve(pepWithModsIsNull ? nullptr : pepsWithMods.Select([&] (std::any b)
		{
			b::BaseSequence;
		}))->ResolvedString;
		s[L"Full Sequence"] = pepWithModsIsNull ? L" " : Resolve(pepWithModsIsNull ? nullptr : pepsWithMods.Select([&] (std::any b)
		{
			b::FullSequence;
		}))->ResolvedString;
		s[L"Essential Sequence"] = pepWithModsIsNull ? L" " : Resolve(pepWithModsIsNull ? nullptr : pepsWithMods.Select([&] (std::any b)
		{
			b::EssentialSequence(ModsToWritePruned);
		}))->ResolvedString;
		s[L"Mods"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::AllModsOneIsNterminus;
		}))->ResolvedString;
		s[L"Mods Chemical Formulas"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any p)
		{
			p::AllModsOneIsNterminus->Select([&] (std::any v)
			{
				v->Value;
			});
		}))->ResolvedString;
		s[L"Mods Combined Chemical Formula"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::AllModsOneIsNterminus->Select([&] (std::any c)
			{
				(dynamic_cast<Modification*>(c->Value));
			});
		}))->ResolvedString;
		s[L"Num Variable Mods"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::NumVariableMods;
		}))->Item1;
		s[L"Missed Cleavages"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::MissedCleavages.ToString(CultureInfo::InvariantCulture);
		}))->ResolvedString;
		s[L"Peptide Monoisotopic Mass"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::MonoisotopicMass;
		}))->ResolvedString;
		s[L"Mass Diff (Da)"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			return psm->getScanPrecursorMass() - b::MonoisotopicMass;
		}))->ResolvedString;
		s[L"Mass Diff (ppm)"] = pepWithModsIsNull ? L" " : ResolveF2(pepsWithMods.Select([&] (std::any b)
		{
			((psm->getScanPrecursorMass() - b::MonoisotopicMass) / b::MonoisotopicMass * 1e6);
		})).ResolvedString;
		s[L"Protein Accession"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::Protein::Accession;
		}), psm->getFullSequence())->ResolvedString;
		s[L"Protein Name"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::Protein->FullName;
		}), psm->getFullSequence())->ResolvedString;
		s[L"Gene Name"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			std::wstring::Join(L", ", b::Protein::GeneNames->Select([&] (std::any d)
			{
				StringHelper::formatSimple(L"{0}:{1}", d::Item1, d::Item2);
			}));
		}), psm->getFullSequence())->ResolvedString;
		s[L"Intersecting Sequence Variations"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			std::wstring::Join(L", ", b::Protein::AppliedSequenceVariations::Where([&] (std::any av)
			{
				IntersectsWithVariation(b, av, false);
			})->Select([&] (std::any av)
			{
				SequenceVariantString(b, av);
			}));
		}))->ResolvedString;
		s[L"Identified Sequence Variations"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			std::wstring::Join(L", ", b::Protein::AppliedSequenceVariations::Where([&] (std::any av)
			{
				IntersectsWithVariation(b, av, true);
			})->Select([&] (std::any av)
			{
				SequenceVariantString(b, av);
			}));
		}))->ResolvedString;
		s[L"Splice Sites"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			std::wstring::Join(L", ", b::Protein::SpliceSites::Where([&] (std::any d)
			{
				Includes(b, d);
			})->Select([&] (std::any d)
			{
				StringHelper::formatSimple(L"{0}-{1}", d::OneBasedBeginPosition.ToString(), d::OneBasedEndPosition.ToString());
			}));
		}))->ResolvedString;
		s[L"Organism Name"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::Protein::Organism;
		}))->Item1;
		s[L"Contaminant"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::Protein::IsContaminant ? L"Y" : L"N";
		}))->Item1;
		s[L"Decoy"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::Protein::IsDecoy ? L"Y" : L"N";
		}))->Item1;
		s[L"Peptide Description"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::PeptideDescription;
		}))->Item1;
		s[L"Start and End Residues In Protein"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			(StringHelper::formatSimple(L"[{0} to {1}]", b::OneBasedStartResidueInProtein.ToString(CultureInfo::InvariantCulture), b::OneBasedEndResidueInProtein.ToString(CultureInfo::InvariantCulture)));
		}), psm->getFullSequence())->ResolvedString;
		s[L"Previous Amino Acid"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::PreviousAminoAcid.ToString();
		}))->ResolvedString;
		s[L"Next Amino Acid"] = pepWithModsIsNull ? L" " : Resolve(pepsWithMods.Select([&] (std::any b)
		{
			b::NextAminoAcid.ToString();
		}))->ResolvedString;

		std::wstring allScores = L" ";
		std::wstring theoreticalsSearched = L" ";
		if (!pepWithModsIsNull && psm->getFdrInfo() != nullptr && psm->getFdrInfo()->getCalculateEValue())
		{
			allScores = std::wstring::Join(L";", psm->getAllScores().Select([&] (std::any p)
			{
				p.ToString(L"F2", CultureInfo::InvariantCulture);
			}));
			theoreticalsSearched = std::to_wstring(psm->getAllScores().size());
		}

		s[L"All Scores"] = allScores;
		s[L"Theoreticals Searched"] = theoreticalsSearched;
		s[L"Decoy/Contaminant/Target"] = pepWithModsIsNull ? L" " : psm->getIsDecoy() ? L"D" : psm->getIsContaminant() ? L"C" : L"T";
	}

	private bool <missing_class_definition>::Includes(PeptideWithSetModifications *pep, SpliceSite *site)
	{
		return pep->OneBasedStartResidueInProtein <= site->OneBasedBeginPosition && pep->OneBasedEndResidueInProtein >= site->OneBasedEndPosition;
	}

	private bool <missing_class_definition>::IntersectsWithVariation(PeptideWithSetModifications *pep, SequenceVariation *appliedVariation, bool checkUnique)
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
			std::wstring originalAtIntersect = appliedVariation->OriginalSequence->substr(intersectOneBasedStart - appliedVariation->OneBasedBeginPosition, intersectSize);
			std::wstring variantAtIntersect = appliedVariation->VariantSequence->substr(intersectOneBasedStart - appliedVariation->OneBasedBeginPosition, intersectSize);
			return originalAtIntersect != variantAtIntersect;
		}
	}

	private std::wstring <missing_class_definition>::SequenceVariantString(PeptideWithSetModifications *p, SequenceVariation *applied)
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
		PeptideWithSetModifications *variantWithAnyMods = new PeptideWithSetModifications(p->Protein, p->DigestionParams, applied->OneBasedBeginPosition, applied->OneBasedEndPosition, p->CleavageSpecificityForFdrCategory, p->PeptideDescription, p->MissedCleavages, modsOnVariantOneIsNTerm, p->NumFixedMods);

		delete variantWithAnyMods;
		return StringHelper::formatSimple(L"{0}{1}{2}", applied->OriginalSequence, applied->OneBasedBeginPosition, variantWithAnyMods->FullSequence);
	}

	public void <missing_class_definition>::AddMatchedIonsData(std::unordered_map<std::wstring, std::wstring> &s, PeptideSpectralMatch *psm)
	{
		bool nullPsm = (psm == nullptr);

		StringBuilder *seriesStringBuilder = new StringBuilder();
		StringBuilder *mzStringBuilder = new StringBuilder();
		StringBuilder *fragmentDaErrorStringBuilder = new StringBuilder();
		StringBuilder *fragmentPpmErrorStringBuilder = new StringBuilder();
		StringBuilder *fragmentIntensityStringBuilder = new StringBuilder();
		std::vector<StringBuilder*> stringBuilders = {seriesStringBuilder, mzStringBuilder, fragmentDaErrorStringBuilder, fragmentPpmErrorStringBuilder, fragmentIntensityStringBuilder};

		if (!nullPsm)
		{
			auto matchedIons = psm->MatchedFragmentIons;
			if (matchedIons == nullptr)
			{
				matchedIons = psm->getPeptidesToMatchingFragments().First()->Value;
			}

			// using ", " instead of "," improves human readability
			const std::wstring delimiter = L", ";

			auto matchedIonsGroupedByProductType = matchedIons->GroupBy([&] (std::any i)
			{
				i::NeutralTheoreticalProduct::ProductType;
			}).OrderBy([&] (std::any i)
			{
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
					p->Append(L"[");
				});

				for (int i = 0; i < products.size(); i++)
				{
					MatchedFragmentIon *ion = products[i];
					std::wstring ionLabel;

					double massError = ion->Mz.ToMass(ion->Charge) - ion->NeutralTheoreticalProduct.NeutralMass;
					double ppmMassError = massError / ion->NeutralTheoreticalProduct.NeutralMass * 1e6;

					if (ion->NeutralTheoreticalProduct->NeutralLoss == 0)
					{
						// no neutral loss
						ionLabel = ion->NeutralTheoreticalProduct.ProductType + L"" + ion->NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + L"+" + ion->Charge;
					}
					else
					{
						// ion label with neutral loss
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
						ionLabel = L"(" + ion->NeutralTheoreticalProduct.ProductType + L"" + ion->NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + L"-" + ion->NeutralTheoreticalProduct.NeutralLoss.ToString(L"F2") + L")" + L"+" + ion->Charge;
					}

					// append ion label
					seriesStringBuilder->append(ionLabel);

					// append experimental m/z
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					mzStringBuilder->append(ionLabel + L":" + ion->Mz.ToString(L"F5"));

					// append absolute mass error
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					fragmentDaErrorStringBuilder->append(ionLabel + L":" + massError.ToString(L"F5"));

					// append ppm mass error
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					fragmentPpmErrorStringBuilder->append(ionLabel + L":" + ppmMassError.ToString(L"F2"));

					// append fragment ion intensity
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					fragmentIntensityStringBuilder->append(ionLabel + L":" + ion->Intensity.ToString(L"F0"));

					// append delimiter ", "
					if (i < products.size() - 1)
					{
						std::for_each(stringBuilders.begin(), stringBuilders.end(), [&] (std::any p)
						{
							p->Append(delimiter);
						});
					}
				}

				// append product type delimiter
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				>
				{
					p->Append("];
				}"));
			}}s[L"Matched Ion Series"] = nullPsm ? L" " : seriesStringBuilder->ToString()->TrimEnd(L';');
//C# TO C++ CONVERTER TODO TASK: The following lambda expression could not be converted:
				stringBuilders.ForEach(p => TangibleLambdaToken86}s[L"Matched Ion Series"];
		s[L"Matched Ion Mass-To-Charge Ratios"] = nullPsm ? L" " : mzStringBuilder->toString()->TrimEnd(L';');
		s[L"Matched Ion Mass Diff (Da)"] = nullPsm ? L" " : fragmentDaErrorStringBuilder->toString()->TrimEnd(L';');
		s[L"Matched Ion Mass Diff (Ppm)"] = nullPsm ? L" " : fragmentPpmErrorStringBuilder->toString()->TrimEnd(L';');
		s[L"Matched Ion Intensities"] = nullPsm ? L" " : fragmentIntensityStringBuilder->toString()->TrimEnd(L';');

		// number of matched ions
		s[L"Matched Ion Counts"] = nullPsm ? L" " : std::to_wstring(psm->MatchedFragmentIons->Count);
		}

//C# TO C++ CONVERTER TODO TASK: Local functions are not converted by C# to C++ Converter:
//	private static void AddMatchScoreData(Dictionary<string, string> s, PeptideSpectralMatch peptide)
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
	}

		delete fragmentIntensityStringBuilder;
		delete fragmentPpmErrorStringBuilder;
		delete fragmentDaErrorStringBuilder;
		delete mzStringBuilder;
		delete seriesStringBuilder;
	}
