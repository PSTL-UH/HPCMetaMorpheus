#include "ModificationAnalysisEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "ModificationAnalysisResults.h"

using namespace Chemistry;
namespace EngineLayer
{
	namespace ModificationAnalysis
	{

		ModificationAnalysisEngine::ModificationAnalysisEngine(std::vector<PeptideSpectralMatch*> &newPsms, CommonParameters *commonParameters, std::vector<std::string> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds), NewPsms(newPsms)
		{
		}

		MetaMorpheusEngineResults *ModificationAnalysisEngine::RunSpecific()
		{
			Status("Running modification analysis...");

			ModificationAnalysisResults *myAnalysisResults = new ModificationAnalysisResults(this);

			auto confidentTargetPsms = NewPsms.Where([&] (std::any b)
			{
			delete myAnalysisResults;
				return b::FdrInfo::QValue <= 0.01 && !b::IsDecoy;
			}).ToList();

			// For the database ones, only need un-ambiguous protein and location in protein
			auto forObserved = confidentTargetPsms.Where([&] (std::any b)
			{
			delete myAnalysisResults;
				return b::ProteinAccession != nullptr && b::OneBasedEndResidueInProtein != nullptr && b::OneBasedStartResidueInProtein != nullptr;
			});

			// For the unambiguously localized ones, need FullSequence and un-ambiguous protein and location in protein
			auto forUnambiguouslyLocalized = confidentTargetPsms.Where([&] (std::any b)
			{
			delete myAnalysisResults;
				return b::FullSequence != nullptr && b::ProteinAccession != nullptr && b::OneBasedEndResidueInProtein != nullptr && b::OneBasedStartResidueInProtein != nullptr;
			});

			//**DEBUG
			std::vector<PeptideSpectralMatch*> toby;
			for (auto psm : confidentTargetPsms)
			{
				if (psm->getFullSequence() != "" && psm->getProteinAccession() != "" && psm->getOneBasedEndResidueInProtein() != nullptr && psm->getOneBasedStartResidueInProtein() != nullptr)
				{
					toby.push_back(psm);
				}
			}

			int i = forUnambiguouslyLocalized->Count() + toby.size()();
			//**END DEBUG

			// For the localized but ambiguous ones, need FullSequence
			auto forAmbiguousButLocalized = confidentTargetPsms.Where([&] (std::any b)
			{
			delete myAnalysisResults;
				return b::FullSequence != nullptr && !(b::ProteinAccession != nullptr && b::OneBasedEndResidueInProtein != nullptr && b::OneBasedStartResidueInProtein != nullptr);
			}).GroupBy([&] (std::any b)
			{
				b::FullSequence;
			});

			// For unlocalized but identified modifications, skip ones with full sequences!
			auto forUnlocalized = confidentTargetPsms.Where([&] (std::any b)
			{
			delete myAnalysisResults;
				return b::BaseSequence != nullptr && b->FullSequence == nullptr && b::ModsIdentified != nullptr;
			}).GroupBy([&] (std::any b)
			{
				(b::BaseSequence, std::string::Join(" ", b::ModsIdentified->Values->OrderBy([&] (std::any c)
				{
			delete myAnalysisResults;
					return c;
				})));
			});

			// For chemical formulas of modifications, skip ones with full sequences and identified mods!
			auto forChemicalFormulas = confidentTargetPsms.Where([&] (std::any b)
			{
			delete myAnalysisResults;
				return b::BaseSequence != nullptr && b->FullSequence == nullptr && b->ModsIdentified == nullptr && b::ModsChemicalFormula != nullptr;
			}).GroupBy([&] (std::any b)
			{
				(b::BaseSequence, b::ModsChemicalFormula);
			});

			// We do not want to double-count modifications. Hence the HashSet!!!
			std::unordered_set<(std::string, std::string, int)*> modsOnProteins;
			for (auto psm : forObserved)
			{
				auto singlePeptide = psm->BestMatchingPeptides.First().Peptide;
				for (auto modInProtein : singlePeptide->Protein.OneBasedPossibleLocalizedModifications.Where([&] (std::any b)
				{
				delete myAnalysisResults;
					return b::Key >= singlePeptide->OneBasedStartResidueInProtein && b::Key <= singlePeptide->OneBasedEndResidueInProtein;
				}))
				{

					for (auto huh : modInProtein->Value)
					{
						modsOnProteins.insert((singlePeptide->Protein.Accession, huh->IdWithMotif, modInProtein::Key));
					}
				}
			}

			// We do not want to double-count modifications. Hence the HashSet!!!
			std::unordered_set<(std::string, std::string, int)*> modsSeenAndLocalized;
			for (auto psm : forUnambiguouslyLocalized)
			{
				auto singlePeptide = psm->BestMatchingPeptides.First().Peptide;
				for (auto nice : singlePeptide->AllModsOneIsNterminus)
				{
					int locInProtein;
					if (nice->Key == 1)
					{
						locInProtein = singlePeptide->OneBasedStartResidueInProtein;
					}
					else if (nice->Key == singlePeptide->Length + 2)
					{
						locInProtein = singlePeptide->OneBasedEndResidueInProtein;
					}
					else
					{
						locInProtein = singlePeptide->OneBasedStartResidueInProtein + nice->Key - 2;
					}
					modsSeenAndLocalized.insert((singlePeptide->Protein.Accession, nice->Value->IdWithMotif, locInProtein));
				}
			}

			// Might have some double counting...
			std::unordered_map<std::string, int> ambiguousButLocalizedModsSeen;
			for (auto representativePsm : forAmbiguousButLocalized->Select([&] (std::any b)
			{
				b::First();
			}))
			{
				for (auto modCountKvp : representativePsm::ModsIdentified)
				{
					if (ambiguousButLocalizedModsSeen.find(modCountKvp->Key) != ambiguousButLocalizedModsSeen.end())
					{
						ambiguousButLocalizedModsSeen[modCountKvp->Key] += modCountKvp->Value;
					}
					else
					{
						ambiguousButLocalizedModsSeen.emplace(modCountKvp->Key, modCountKvp->Value);
					}
				}
			}

			// Might have some double counting...
			std::unordered_map<std::string, int> unlocalizedMods;
			for (auto representativePsm : forUnlocalized->Select([&] (std::any b)
			{
				b::First();
			}))
			{
				for (auto modCountKvp : representativePsm::ModsIdentified)
				{
					if (unlocalizedMods.find(modCountKvp->Key) != unlocalizedMods.end())
					{
						unlocalizedMods[modCountKvp->Key] += modCountKvp->Value;
					}
					else
					{
						unlocalizedMods.emplace(modCountKvp->Key, modCountKvp->Value);
					}
				}
			}

			// Might have some double counting...
			std::unordered_map<ChemicalFormula*, int> unlocalizedFormulas;
			for (auto representativePsm : forChemicalFormulas->Select([&] (std::any b)
			{
				b::First();
			}))
			{
				if (unlocalizedFormulas.find(representativePsm::ModsChemicalFormula) != unlocalizedFormulas.end())
				{
					unlocalizedFormulas[representativePsm::ModsChemicalFormula] += 1;
				}
				else
				{
					unlocalizedFormulas.emplace(representativePsm::ModsChemicalFormula, 1);
				}
			}

			myAnalysisResults->setCountOfEachModSeenOnProteins(modsOnProteins.GroupBy([&] (std::any b)
			{
				b::Item2;
			}).ToDictionary([&] (std::any b)
			{
				b::Key;
			}, [&] (std::any b)
			{
				b->Count();
			}));
			myAnalysisResults->setCountOfModsSeenAndLocalized(modsSeenAndLocalized.GroupBy([&] (std::any b)
			{
				b::Item2;
			}).ToDictionary([&] (std::any b)
			{
				b::Key;
			}, [&] (std::any b)
			{
				b->Count();
			}));
			myAnalysisResults->setCountOfAmbiguousButLocalizedModsSeen(ambiguousButLocalizedModsSeen);
			myAnalysisResults->setCountOfUnlocalizedMods(unlocalizedMods);
			myAnalysisResults->setCountOfUnlocalizedFormulas(unlocalizedFormulas);

//C# TO C++ CONVERTER TODO TASK: A 'delete myAnalysisResults' statement was not added since myAnalysisResults was used in a 'return' or 'throw' statement.
			return myAnalysisResults;
		}
	}
}
