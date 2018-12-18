#include "GptmdEngine.h"
#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"
#include "GptmdResults.h"

using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
	namespace Gptmd
	{

		GptmdEngine::GptmdEngine(std::vector<PeptideSpectralMatch*> &allIdentifications, std::vector<Modification*> &gptmdModifications, std::vector<std::tuple<double, double>> &combos, std::unordered_map<std::wstring, Tolerance*> &filePathToPrecursorMassTolerance, CommonParameters *commonParameters, std::vector<std::wstring> &nestedIds) : MetaMorpheusEngine(commonParameters, nestedIds), AllIdentifications(allIdentifications), Combos(combos), GptmdModifications(gptmdModifications), FilePathToPrecursorMassTolerance(filePathToPrecursorMassTolerance)
		{
		}

		bool GptmdEngine::ModFits(Modification *attemptToLocalize, Protein *protein, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex)
		{
			//the peptideOneBasedIndex and proteinOneBasedIndex are for the position of the modification on the sequence

			auto motif = attemptToLocalize->Target;
			// First find the capital letter...
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			auto hehe = motif->ToString()->find(motif->ToString()->First([&] (std::any b)
			{
				std::isupper(b);
			}));

			auto proteinToMotifOffset = proteinOneBasedIndex - hehe - 1;
			auto indexUp = 0;
			// Look up starting at and including the capital letter
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			while (indexUp < motif->ToString()->Length)
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				if (indexUp + proteinToMotifOffset < 0 || indexUp + proteinToMotifOffset >= protein->Length || (!std::toupper(motif->ToString()[indexUp])->Equals(L'X') && !std::toupper(motif->ToString()[indexUp])->Equals(protein->BaseSequence[indexUp + proteinToMotifOffset])))
				{
					return false;
				}
				indexUp++;
			}
			if (attemptToLocalize->LocationRestriction == L"Anywhere.")
			{
				return true;
			}
			if (attemptToLocalize->LocationRestriction == L"N-terminal." && (proteinOneBasedIndex <= 2))
			{
				return true;
			}
			if (attemptToLocalize->LocationRestriction == L"Peptide N-terminal." && peptideOneBasedIndex == 1)
			{
				return true;
			}
			if (attemptToLocalize->LocationRestriction == L"Peptide C-terminal." && peptideOneBasedIndex == peptideLength)
			{
				return true;
			}
			if (attemptToLocalize->LocationRestriction == L"C-terminal." && proteinOneBasedIndex == protein->Length)
			{
				return true;
			}
			return false;
		}

		MetaMorpheusEngineResults *GptmdEngine::RunSpecific()
		{
			auto modDict = std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>();

			int modsAdded = 0;
			//foreach peptide in each psm and for each modification that matches the notch,
			//add that modification to every allowed residue
			for (auto psm : AllIdentifications.Where([&] (std::any b)
			{
				return b::FdrInfo::QValueNotch <= 0.05 && !b::IsDecoy;
			}))
			{
				// get file-specific precursor tolerance
				Tolerance *precursorMassTolerance = FilePathToPrecursorMassTolerance[psm::FullFilePath];

				// get mods to annotate database with
				for (auto pepWithSetMods : psm::BestMatchingPeptides->Select([&] (std::any v)
				{
					v::Peptide;
				}))
				{
					for (auto mod : GetPossibleMods(psm::ScanPrecursorMass, GptmdModifications, Combos, precursorMassTolerance, pepWithSetMods))
					{
						auto isVariantProtein = pepWithSetMods::Protein != pepWithSetMods::Protein::NonVariantProtein;

						for (int i = 0; i < pepWithSetMods->Length; i++)
						{
							int indexInProtein = pepWithSetMods::OneBasedStartResidueInProtein + i;

							if (ModFits(mod, pepWithSetMods::Protein, i + 1, pepWithSetMods->Length, indexInProtein))
							{
								// if not a variant protein, index to base protein sequence
								if (!isVariantProtein)
								{
									AddIndexedMod(modDict, pepWithSetMods::Protein::Accession, std::tuple<int, Modification*>(indexInProtein, mod));
									modsAdded++;
								}

								// if a variant protein, index to variant protein if on variant, or to the original protein if not
								else
								{
									bool foundSite = false;
									int offset = 0;
									for (auto variant : pepWithSetMods::Protein::AppliedSequenceVariations::OrderBy([&] (std::any v)
									{
										v::OneBasedBeginPosition;
									}))
									{
										bool modIsBeforeVariant = indexInProtein < variant::OneBasedBeginPosition + offset;
										bool modIsOnVariant = variant::OneBasedBeginPosition + offset <= indexInProtein && indexInProtein <= variant::OneBasedEndPosition + offset;

										// if a variant protein and the mod is on the variant, index to the variant protein sequence
										if (modIsOnVariant)
										{
											AddIndexedMod(modDict, pepWithSetMods::Protein::Accession, std::tuple<int, Modification*>(indexInProtein, mod));
											modsAdded++;
											foundSite = true;
											break;
										}

										// otherwise back calculate the index to the original protein sequence
										if (modIsBeforeVariant)
										{
											AddIndexedMod(modDict, pepWithSetMods::Protein::NonVariantProtein::Accession, std::tuple<int, Modification*>(indexInProtein - offset, mod));
											modsAdded++;
											foundSite = true;
											break;
										}

										offset += variant::VariantSequence->Length - variant::OriginalSequence->Length;
									}
									if (!foundSite)
									{
										AddIndexedMod(modDict, pepWithSetMods::Protein::NonVariantProtein::Accession, std::tuple<int, Modification*>(indexInProtein - offset, mod));
										modsAdded++;
									}
								}
							}
						}
					}
				}
			}

			return new GptmdResults(this, modDict, modsAdded);
		}

		void GptmdEngine::AddIndexedMod(std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> &modDict, const std::wstring &proteinAccession, std::tuple<int, Modification*> &indexedMod)
		{
			TValue hash;
			std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>>::const_iterator modDict_iterator = modDict.find(proteinAccession);
			if (modDict_iterator != modDict.end())
			{
				hash = modDict_iterator->second;
				hash->Add(indexedMod);
			}
			else
			{
				hash = modDict_iterator->second;
				modDict[proteinAccession] = {indexedMod};
			}
		}

		std::vector<Modification*> GptmdEngine::GetPossibleMods(double totalMassToGetTo, std::vector<Modification*> &allMods, std::vector<std::tuple<double, double>> &combos, Tolerance *precursorTolerance, PeptideWithSetModifications *peptideWithSetModifications)
		{
			for (auto Mod : allMods.Where([&] (std::any b)
			{
				return b->ValidModification == true;
			}))
			{
				if (precursorTolerance->Within(totalMassToGetTo, peptideWithSetModifications->MonoisotopicMass + static_cast<double>(Mod::MonoisotopicMass)))
				{
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
					yield return Mod;
				}
				for (auto modOnPsm : peptideWithSetModifications->AllModsOneIsNterminus->Values->Where([&] (std::any b)
				{
					return b->ValidModification == true;
				}))
				{
					if (modOnPsm::Target->Equals(Mod::Target))
					{
						if (precursorTolerance->Within(totalMassToGetTo, peptideWithSetModifications->MonoisotopicMass + static_cast<double>(Mod::MonoisotopicMass) - static_cast<double>(modOnPsm::MonoisotopicMass)))
						{

							//TODO: not necessarily here. I think we're creating ambiguity. If we're going to add a gptmd mod to a peptide that already has that mod, then we need info to suggest that it is at a postion other than that in the database. could be presence of frag for unmodified or presence of frag with modified at alternative location.
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
							yield return Mod;
						}
					}
				}
			}

			for (auto combo : combos)
			{
				auto m1 = combo.Item1;
				auto m2 = combo.Item2;
				auto combined = m1 + m2;
				if (precursorTolerance->Within(totalMassToGetTo, peptideWithSetModifications->MonoisotopicMass + combined))
				{
					for (auto mod : GetPossibleMods(totalMassToGetTo - m1, allMods, combos, precursorTolerance, peptideWithSetModifications))
					{
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
						yield return mod;
					}
					for (auto mod : GetPossibleMods(totalMassToGetTo - m2, allMods, combos, precursorTolerance, peptideWithSetModifications))
					{
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
						yield return mod;
					}
				}
			}
		}
	}
}
