#include "CrosslinkedPeptides.h"
#include "Crosslinker.h"

using namespace MassSpectrometry;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
	namespace CrosslinkSearch
	{

		std::vector<std::tuple<int, std::vector<Product*>>> CrosslinkedPeptide::XlGetTheoreticalFragments(DissociationType *dissociationType, Crosslinker *crosslinker, std::vector<int> &possibleCrosslinkerPositions, double otherPeptideMass, PeptideWithSetModifications *peptide)
		{
			std::vector<double> massesToLocalize;
			if (crosslinker->getCleavable())
			{
				massesToLocalize.push_back(crosslinker->getCleaveMassShort());
				massesToLocalize.push_back(crosslinker->getCleaveMassLong());
			}
			else
			{
				massesToLocalize.push_back(crosslinker->getTotalMass() + otherPeptideMass);
			}

			for (auto crosslinkerPosition : possibleCrosslinkerPositions)
			{
				std::vector<Product*> theoreticalProducts;
				std::unordered_set<double> masses;

				for (auto massToLocalize : massesToLocalize)
				{
					std::unordered_map<int, Modification*> testMods =
					{
						{crosslinkerPosition + 1, new Modification(_monoisotopicMass: massToLocalize)}
					};
					auto testPeptide = new PeptideWithSetModifications(peptide->Protein, peptide->DigestionParams, peptide->OneBasedStartResidueInProtein, peptide->OneBasedEndResidueInProtein, peptide->CleavageSpecificityForFdrCategory, peptide->PeptideDescription, peptide->MissedCleavages, testMods, peptide->NumFixedMods);

					// add fragmentation ions for this crosslinker position guess
					for (auto fragment : testPeptide->Fragment(dissociationType, FragmentationTerminus::Both))
					{
						if (!std::find(masses.begin(), masses.end(), fragment->NeutralMass) != masses.end())
						{
							theoreticalProducts.push_back(fragment);
							masses.insert(fragment->NeutralMass);
						}
					}

					// add signature ions
					if (crosslinker->getCleavable())
					{
						Product tempVar(ProductType::M, new NeutralTerminusFragment(FragmentationTerminus::None, peptide->MonoisotopicMass + massToLocalize, peptide->Length, peptide->Length), 0);
						theoreticalProducts.push_back(&tempVar);
					}

					delete testPeptide;
				}

//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
				yield return std::tuple<int, std::vector<Product*>>(crosslinkerPosition, theoreticalProducts);
			}
		}

		std::unordered_map<std::tuple<int, int>, std::vector<Product*>> CrosslinkedPeptide::XlLoopGetTheoreticalFragments(DissociationType *dissociationType, Modification *loopMass, std::vector<int> &modPos, PeptideWithSetModifications *peptide)
		{
			std::unordered_map<std::tuple<int, int>, std::vector<Product*>> AllTheoreticalFragmentsLists;
			auto originalFragments = peptide->Fragment(dissociationType, FragmentationTerminus::Both).ToList();

			for (auto position1 : modPos)
			{
				for (auto position2 : modPos)
				{
					if (position2 <= position1)
					{
						continue;
					}

					// add N and C terminal fragments that do not contain the loop
					std::tuple<int, int> loopPositions = std::tuple<int, int>(position1, position2);
					std::vector<Product*> loopFragments = originalFragments.Where([&] (std::any p)
					{
						return p::TerminusFragment->Terminus == FragmentationTerminus::N && p::TerminusFragment::AminoAcidPosition < position1 || p::TerminusFragment->Terminus == FragmentationTerminus::C && p::TerminusFragment::AminoAcidPosition > position2;
					}).ToList();

					// add N-terminal fragments containing the loop
					std::unordered_map<int, Modification*> modDict;
					if (peptide->AllModsOneIsNterminus.Any())
					{
						double combinedModMass = loopMass->MonoisotopicMass->Value + peptide->AllModsOneIsNterminus.Where([&] (std::any v)
						{
							return v::Key <= position2 + 1;
						}).Sum([&] (std::any p)
						{
							p->Value->MonoisotopicMass->Value;
						});
						Modification *combined = new Modification(_monoisotopicMass: combinedModMass);
						modDict.emplace(position1 + 1, combined);

						for (auto mod : peptide->AllModsOneIsNterminus.Where([&] (std::any m)
						{
//C# TO C++ CONVERTER TODO TASK: A 'delete combined' statement was not added since combined was passed to a method or constructor. Handle memory management manually.
							return m::Key > position2 + 1;
						}))
						{
							modDict.emplace(mod::Key, mod->Value);
						}

//C# TO C++ CONVERTER TODO TASK: A 'delete combined' statement was not added since combined was passed to a method or constructor. Handle memory management manually.
					}
					else
					{
						modDict.emplace(position1 + 1, loopMass);
					}
					PeptideWithSetModifications *peptideWithLoop = new PeptideWithSetModifications(peptide->Protein, peptide->DigestionParams, peptide->OneBasedStartResidueInProtein, peptide->OneBasedEndResidueInProtein, peptide->CleavageSpecificityForFdrCategory, peptide->PeptideDescription, peptide->MissedCleavages, modDict, peptide->NumFixedMods);
					loopFragments.AddRange(peptideWithLoop->Fragment(dissociationType, FragmentationTerminus::Both).Where([&] (std::any p)
					{
					delete peptideWithLoop;
						return p::TerminusFragment->Terminus == FragmentationTerminus::N && p::TerminusFragment::AminoAcidPosition >= position2;
					}));

					// add C-terminal fragments containing the loop
					modDict.clear();
					if (peptide->AllModsOneIsNterminus.Any())
					{
						double combinedModMass = loopMass->MonoisotopicMass->Value + peptide->AllModsOneIsNterminus.Where([&] (std::any v)
						{
						delete peptideWithLoop;
							return v::Key >= position1 + 1;
						}).Sum([&] (std::any p)
						{
							p->Value->MonoisotopicMass->Value;
						});
						Modification *combined = new Modification(_monoisotopicMass: combinedModMass);
						modDict.emplace(position2 + 1, combined);

						for (auto mod : peptide->AllModsOneIsNterminus.Where([&] (std::any m)
						{
//C# TO C++ CONVERTER TODO TASK: A 'delete combined' statement was not added since combined was passed to a method or constructor. Handle memory management manually.
						delete peptideWithLoop;
							return m::Key < position1 + 1;
						}))
						{
							modDict.emplace(mod::Key, mod->Value);
						}

//C# TO C++ CONVERTER TODO TASK: A 'delete combined' statement was not added since combined was passed to a method or constructor. Handle memory management manually.
					}
					else
					{
						modDict.emplace(position2 + 1, loopMass);
					}
					peptideWithLoop = new PeptideWithSetModifications(peptide->Protein, peptide->DigestionParams, peptide->OneBasedStartResidueInProtein, peptide->OneBasedEndResidueInProtein, peptide->CleavageSpecificityForFdrCategory, peptide->PeptideDescription, peptide->MissedCleavages, modDict, peptide->NumFixedMods);
					loopFragments.AddRange(peptideWithLoop->Fragment(dissociationType, FragmentationTerminus::Both).Where([&] (std::any p)
					{
					delete peptideWithLoop;
						return p::TerminusFragment->Terminus == FragmentationTerminus::C && p::TerminusFragment::AminoAcidPosition <= position1;
					}));

					AllTheoreticalFragmentsLists.emplace(loopPositions, loopFragments);

					delete peptideWithLoop;
				}
			}

			return AllTheoreticalFragmentsLists;
		}
	}
}
