#include "ProteinScoringAndFdrEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../ProteinParsimony/ProteinGroup.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "ProteinScoringAndFdrResults.h"

using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

	ProteinScoringAndFdrEngine::ProteinScoringAndFdrEngine(std::vector<ProteinGroup*> &proteinGroups, std::vector<PeptideSpectralMatch*> &newPsms, bool noOneHitWonders, bool treatModPeptidesAsDifferentPeptides, bool mergeIndistinguishableProteinGroups, CommonParameters *commonParameters, std::vector<std::wstring> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds), NewPsms(newPsms), NoOneHitWonders(noOneHitWonders), TreatModPeptidesAsDifferentPeptides(treatModPeptidesAsDifferentPeptides), MergeIndistinguishableProteinGroups(mergeIndistinguishableProteinGroups), ProteinGroups(proteinGroups)
	{
	}

	MetaMorpheusEngineResults *ProteinScoringAndFdrEngine::RunSpecific()
	{
		ProteinScoringAndFdrResults *myAnalysisResults = new ProteinScoringAndFdrResults(this);

		ScoreProteinGroups(ProteinGroups, NewPsms);
		myAnalysisResults->SortedAndScoredProteinGroups = DoProteinFdr(ProteinGroups);

//C# TO C++ CONVERTER TODO TASK: A 'delete myAnalysisResults' statement was not added since myAnalysisResults was used in a 'return' or 'throw' statement.
		return myAnalysisResults;
	}

	std::wstring ProteinScoringAndFdrEngine::StripDecoyIdentifier(const std::wstring &proteinGroupName)
	{
		return proteinGroupName.find(L"DECOY_") != std::wstring::npos ? StringHelper::replace(proteinGroupName, L"DECOY_", L"") : proteinGroupName;
	}

	void ProteinScoringAndFdrEngine::ScoreProteinGroups(std::vector<ProteinGroup*> &proteinGroups, std::vector<PeptideSpectralMatch*> &psmList)
	{
		// add each protein groups PSMs
		auto peptideToPsmMatching = std::unordered_map<PeptideWithSetModifications*, std::unordered_set<PeptideSpectralMatch*>>();
		for (auto psm : psmList)
		{
			if (psm->getFdrInfo()->getQValueNotch() <= 0.01 && psm->getFdrInfo()->getQValue() <= 0.01)
			{
				if ((TreatModPeptidesAsDifferentPeptides && psm->getFullSequence() != L"") || (!TreatModPeptidesAsDifferentPeptides && psm->getBaseSequence() != L""))
				{
					for (auto pepWithSetMods : psm->BestMatchingPeptides->Select([&] (std::any p)
					{
						p::Peptide;
					}))
					{
						HashSet<PeptideSpectralMatch*> psmsForThisPeptide;
						std::unordered_map<PeptideWithSetModifications*, std::unordered_set<PeptideSpectralMatch*>>::const_iterator peptideToPsmMatching_iterator = peptideToPsmMatching.find(pepWithSetMods);
						if (peptideToPsmMatching_iterator == peptideToPsmMatching.end())
						{
							psmsForThisPeptide = peptideToPsmMatching_iterator->second;
							peptideToPsmMatching.emplace(pepWithSetMods, std::unordered_set<PeptideSpectralMatch*> {psm});
						}
						else
						{
							psmsForThisPeptide = peptideToPsmMatching_iterator->second;
							psmsForThisPeptide->Add(psm);
						}
					}
				}
			}
		}

		for (auto proteinGroup : proteinGroups)
		{
			std::vector<PeptideWithSetModifications*> pepsToRemove;
			for (auto peptide : proteinGroup->getAllPeptides())
			{
				// build PSM list for scoring
				HashSet<PeptideSpectralMatch*> psms;
				std::unordered_map<PeptideWithSetModifications*, std::unordered_set<PeptideSpectralMatch*>>::const_iterator peptideToPsmMatching_iterator = peptideToPsmMatching.find(peptide);
				if (peptideToPsmMatching_iterator != peptideToPsmMatching.end())
				{
					psms = peptideToPsmMatching_iterator->second;
					proteinGroup->getAllPsmsBelowOnePercentFDR().UnionWith(psms);
				}
				else
				{
					psms = peptideToPsmMatching_iterator->second;
					pepsToRemove.push_back(peptide);
				}
			}

			proteinGroup->getAllPeptides().ExceptWith(pepsToRemove);
			proteinGroup->getUniquePeptides().ExceptWith(pepsToRemove);
		}

		// score the group
		for (auto proteinGroup : proteinGroups)
		{
			proteinGroup->Score();
		}

		if (MergeIndistinguishableProteinGroups)
		{
			// merge protein groups that are indistinguishable after scoring
			auto pg = proteinGroups.OrderByDescending([&] (std::any p)
			{
				p::ProteinGroupScore;
			}).ToList();
			for (int i = 0; i < (pg.size() - 1); i++)
			{
				if (pg[i].ProteinGroupScore == pg[i + 1].ProteinGroupScore && pg[i].ProteinGroupScore != 0)
				{
					auto pgsWithThisScore = pg.Where([&] (std::any p)
					{
						return p->ProteinGroupScore == pg[i].ProteinGroupScore;
					}).ToList();

					// check to make sure they have the same peptides, then merge them
					for (auto p : pgsWithThisScore)
					{
						auto seqs1 = std::unordered_set<std::wstring>(p.AllPeptides->Select([&] (std::any x)
						{
							return x::FullSequence + x::DigestionParams::Protease;
						}));
						auto seqs2 = std::unordered_set<std::wstring>(pg[i].AllPeptides->Select([&] (std::any x)
						{
							return x::FullSequence + x::DigestionParams::Protease;
						}));

						if (p != pg[i] && seqs1.SetEquals(seqs2))
						{
							pg[i].MergeProteinGroupWith(p);
						}
					}
				}
			}
		}

		// remove empty protein groups (peptides were too poor quality or group was merged)
		proteinGroups.RemoveAll([&] (std::any p)
		{
			return p->ProteinGroupScore == 0;
		});

		// calculate sequence coverage
		for (auto proteinGroup : proteinGroups)
		{
			proteinGroup->CalculateSequenceCoverage();
		}
	}

	std::vector<ProteinGroup*> ProteinScoringAndFdrEngine::DoProteinFdr(std::vector<ProteinGroup*> &proteinGroups)
	{
		if (NoOneHitWonders)
		{
			if (TreatModPeptidesAsDifferentPeptides)
			{
				proteinGroups = proteinGroups.Where([&] (std::any p)
				{
					return p::IsDecoy || (std::unordered_set<std::wstring>(p::AllPeptides->Select([&] (std::any x)
					{
						x::FullSequence;
					})))->size() > 1;
				}).ToList();
			}
			else
			{
				proteinGroups = proteinGroups.Where([&] (std::any p)
				{
					return p::IsDecoy || (std::unordered_set<std::wstring>(p::AllPeptides->Select([&] (std::any x)
					{
						x::BaseSequence;
					})))->size() > 1;
				}).ToList();
			}
		}

		// pair decoys and targets by accession
		// then use the best peptide notch-QValue as the score for the protein group
		std::unordered_map<std::wstring, std::vector<ProteinGroup*>> accessionToProteinGroup;
		for (auto pg : proteinGroups)
		{
			for (auto protein : pg->getProteins())
			{
				std::wstring stippedAccession = StripDecoyIdentifier(protein->Accession);

				List<ProteinGroup*> groups;
				std::unordered_map<std::wstring, std::vector<ProteinGroup*>>::const_iterator accessionToProteinGroup_iterator = accessionToProteinGroup.find(stippedAccession);
				if (accessionToProteinGroup_iterator != accessionToProteinGroup.end())
				{
					groups = accessionToProteinGroup_iterator->second;
					groups->Add(pg);
				}
				else
				{
					groups = accessionToProteinGroup_iterator->second;
					accessionToProteinGroup.emplace(stippedAccession, std::vector<ProteinGroup*> {pg});
				}
			}

			pg->setBestPeptideScore(pg->getAllPsmsBelowOnePercentFDR().Max([&] (std::any psm)
			{
				psm::Score;
			}));
			pg->setBestPeptideQValue(pg->getAllPsmsBelowOnePercentFDR().Min([&] (std::any psm)
			{
				psm::FdrInfo::QValueNotch;
			}));
		}

		// pick the best notch-QValue for each paired accession
		for (auto accession : accessionToProteinGroup)
		{
			if (accession.Value->Count > 1)
			{
				auto pgList = accession.Value->OrderBy([&] (std::any p)
				{
					p::BestPeptideQValue;
				}).ThenByDescending([&] (std::any p)
				{
					p::BestPeptideScore;
				}).ToList();
				auto pgToUse = pgList.front(); // pick lowest notch QValue and remove the rest
				pgList.Remove(pgToUse);
				proteinGroups = proteinGroups.Except(pgList).ToList();
			}
		}

		// order protein groups by notch-QValue
		auto sortedProteinGroups = proteinGroups.OrderBy([&] (std::any b)
		{
			b::BestPeptideQValue;
		}).ThenByDescending([&] (std::any p)
		{
			p::BestPeptideScore;
		}).ToList();

		// calculate protein QValues
		int cumulativeTarget = 0;
		int cumulativeDecoy = 0;

		for (auto proteinGroup : sortedProteinGroups)
		{
			if (proteinGroup.IsDecoy)
			{
				cumulativeDecoy++;
			}
			else
			{
				cumulativeTarget++;
			}

			proteinGroup.CumulativeTarget = cumulativeTarget;
			proteinGroup.CumulativeDecoy = cumulativeDecoy;
			proteinGroup.QValue = static_cast<double>(cumulativeDecoy) / cumulativeTarget;
		}

		return sortedProteinGroups;
	}
}
