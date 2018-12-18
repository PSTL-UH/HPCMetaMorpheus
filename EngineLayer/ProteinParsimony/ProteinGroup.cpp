#include "ProteinGroup.h"
#include "../PeptideSpectralMatch.h"
#include "../GlobalVariables.h"

using namespace FlashLFQ;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

	ProteinGroup::ProteinGroup(std::unordered_set<Protein*> &proteins, std::unordered_set<PeptideWithSetModifications*> &peptides, std::unordered_set<PeptideWithSetModifications*> &uniquePeptides)
	{
		setProteins(proteins);
		ListOfProteinsOrderedByAccession = getProteins().OrderBy([&] (std::any p)
		{
			p::Accession;
		}).ToList();
		setProteinGroupName(std::wstring::Join(L"|", ListOfProteinsOrderedByAccession.Select([&] (std::any p)
		{
			p::Accession;
		})));
		setAllPeptides(peptides);
		setUniquePeptides(uniquePeptides);
		setAllPsmsBelowOnePercentFDR(std::unordered_set<PeptideSpectralMatch*>());
		setSequenceCoveragePercent(std::vector<double>());
		setSequenceCoverageDisplayList(std::vector<std::wstring>());
		setSequenceCoverageDisplayListWithMods(std::vector<std::wstring>());
		setProteinGroupScore(0);
		setBestPeptideScore(0);
		setQValue(0);
		IsDecoy = false;
		IsContaminant = false;
		setModsInfo(std::vector<std::wstring>());

		// if any of the proteins in the protein group are decoys, the protein group is a decoy
		for (auto protein : proteins)
		{
			if (protein->IsDecoy)
			{
				IsDecoy = true;
				break;
			}
			if (protein->IsContaminant)
			{
				IsContaminant = true;
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

	std::wstring ProteinGroup::getProteinGroupName() const
	{
		return privateProteinGroupName;
	}

	void ProteinGroup::setProteinGroupName(const std::wstring &value)
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

	std::vector<std::wstring> ProteinGroup::getSequenceCoverageDisplayList() const
	{
		return privateSequenceCoverageDisplayList;
	}

	void ProteinGroup::setSequenceCoverageDisplayList(const std::vector<std::wstring> &value)
	{
		privateSequenceCoverageDisplayList = value;
	}

	std::vector<std::wstring> ProteinGroup::getSequenceCoverageDisplayListWithMods() const
	{
		return privateSequenceCoverageDisplayListWithMods;
	}

	void ProteinGroup::setSequenceCoverageDisplayListWithMods(const std::vector<std::wstring> &value)
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

	std::vector<std::wstring> ProteinGroup::getModsInfo() const
	{
		return privateModsInfo;
	}

	void ProteinGroup::setModsInfo(const std::vector<std::wstring> &value)
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

	std::wstring ProteinGroup::GetTabSeparatedHeader()
	{
		auto sb = new StringBuilder();
		sb->append(L"Protein Accession" + StringHelper::toString(L'\t'));
		sb->append(L"Gene" + StringHelper::toString(L'\t'));
		sb->append(L"Organism" + StringHelper::toString(L'\t'));
		sb->append(L"Protein Full Name" + StringHelper::toString(L'\t'));
		sb->append(L"Protein Unmodified Mass" + StringHelper::toString(L'\t'));
		sb->append(L"Number of Proteins in Group" + StringHelper::toString(L'\t'));
		sb->append(L"Unique Peptides" + StringHelper::toString(L'\t'));
		sb->append(L"Shared Peptides" + StringHelper::toString(L'\t'));
		sb->append(L"Number of Peptides" + StringHelper::toString(L'\t'));
		sb->append(L"Number of Unique Peptides" + StringHelper::toString(L'\t'));
		sb->append(L"Sequence Coverage %" + StringHelper::toString(L'\t'));
		sb->append(L"Sequence Coverage" + StringHelper::toString(L'\t'));
		sb->append(L"Sequence Coverage with Mods" + StringHelper::toString(L'\t'));
		sb->append(std::wstring(L"Modification Info List") + L"\t");
		if (getFilesForQuantification().size() > 0)
		{
			for (int i = 0; i < getFilesForQuantification().size(); i++)
			{
				sb->append(L"Intensity_" + getFilesForQuantification()[i]->FilenameWithoutExtension + L'\t');
			}
		}
		sb->append(L"Number of PSMs" + StringHelper::toString(L'\t'));
		sb->append(L"Protein Decoy/Contaminant/Target" + StringHelper::toString(L'\t'));
		sb->append(L"Protein Cumulative Target" + StringHelper::toString(L'\t'));
		sb->append(L"Protein Cumulative Decoy" + StringHelper::toString(L'\t'));
		sb->append(L"Protein QValue" + StringHelper::toString(L'\t'));
		sb->append(L"Best Peptide Score" + StringHelper::toString(L'\t'));
		sb->append(L"Best Peptide Notch QValue");

		delete sb;
		return sb->toString();
	}

	std::wstring ProteinGroup::ToString()
	{
		auto sb = new StringBuilder();

		// list of protein accession numbers
		sb->append(getProteinGroupName());
		sb->append(L"\t");

		// genes
		sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", ListOfProteinsOrderedByAccession.Select([&] (std::any p)
		{
			p::GeneNames->Select([&] (std::any x)
			{
				x::Item2;
			}).FirstOrDefault();
		}))));
		sb->append(L"\t");

		// organisms
		sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", ListOfProteinsOrderedByAccession.Select([&] (std::any p)
		{
			p::Organism;
		}).Distinct())));
		sb->append(L"\t");

		// list of protein names
		sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", ListOfProteinsOrderedByAccession.Select([&] (std::any p)
		{
			p->FullName;
		}).Distinct())));
		sb->append(L"\t");

		// list of masses
		auto sequences = ListOfProteinsOrderedByAccession.Select([&] (std::any p)
		{
			p::BaseSequence;
		}).Distinct();
		std::vector<double> masses;
		for (auto sequence : sequences)
		{
			try
			{
				Proteomics::AminoAcidPolymer::Peptide tempVar(sequence);
				masses.push_back((&tempVar)->MonoisotopicMass);
			}
			catch (const std::runtime_error &e1)
			{
				masses.push_back(NAN);
			}
		}
		sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", masses)));
		sb->append(L"\t");

		// number of proteins in group
		sb->append(L"" + std::to_wstring(getProteins().size()));
		sb->append(L"\t");

		// list of unique peptides
		if (!getDisplayModsOnPeptides())
		{
			sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", getUniquePeptides().Select([&] (std::any p)
			{
				p::BaseSequence;
			}).Distinct())));
		}
		else
		{
			sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", getUniquePeptides().Select([&] (std::any p)
			{
				p::FullSequence;
			}).Distinct())));
		}
		sb->append(L"\t");

		// list of shared peptides
		auto SharedPeptides = getAllPeptides().Except(getUniquePeptides());
		if (!getDisplayModsOnPeptides())
		{
			sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", SharedPeptides->Select([&] (std::any p)
			{
				p::BaseSequence;
			}).Distinct())));
		}
		else
		{
			sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", SharedPeptides->Select([&] (std::any p)
			{
				p::FullSequence;
			}).Distinct())));
		}
		sb->append(L"\t");

		// number of peptides
		if (!getDisplayModsOnPeptides())
		{
			sb->append(L"" + getAllPeptides().Select([&] (std::any p)
			{
				p::BaseSequence;
			}).Distinct()->Count());
		}
		else
		{
			sb->append(L"" + getAllPeptides().Select([&] (std::any p)
			{
				p::FullSequence;
			}).Distinct()->Count());
		}
		sb->append(L"\t");

		// number of unique peptides
		if (!getDisplayModsOnPeptides())
		{
			sb->append(L"" + getUniquePeptides().Select([&] (std::any p)
			{
				p::BaseSequence;
			}).Distinct()->Count());
		}
		else
		{
			sb->append(L"" + getUniquePeptides().Select([&] (std::any p)
			{
				p::FullSequence;
			}).Distinct()->Count());
		}
		sb->append(L"\t");

		// sequence coverage percent
		sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", getSequenceCoveragePercent().Select([&] (std::any p)
		{
			std::wstring::Format(std::wstring(L"{0:0}") + L"%", (p * 100));
		}))));
		sb->append(L"\t");

		// sequence coverage
		sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", getSequenceCoverageDisplayList())));
		sb->append(L"\t");

		// sequence coverage with mods
		sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", getSequenceCoverageDisplayListWithMods())));
		sb->append(L"\t");

		//Detailed mods information list
		sb->append(GlobalVariables::CheckLengthOfOutput(std::wstring::Join(L"|", getModsInfo())));
		sb->append(L"\t");

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
					sb->append(L"");
				}
				sb->append(L"\t");
			}
		}

		// number of PSMs for listed peptides
		sb->append(L"" + std::to_wstring(getAllPsmsBelowOnePercentFDR().size()));
		sb->append(L"\t");

		// isDecoy
		if (getIsDecoy())
		{
			sb->append(L"D");
		}
		else if (getIsContaminant())
		{
			sb->append(L"C");
		}
		else
		{
			sb->append(L"T");
		}
		sb->append(L"\t");

		// cumulative target
		sb->append(getCumulativeTarget());
		sb->append(L"\t");

		// cumulative decoy
		sb->append(getCumulativeDecoy());
		sb->append(L"\t");

		// q value
		sb->append(getQValue());
		sb->append(L"\t");

		// best peptide score
		sb->append(getBestPeptideScore());
		sb->append(L"\t");

		// best peptide q value
		sb->append(getBestPeptideQValue());
		sb->append(L"\t");

		delete sb;
		return sb->toString();
	}

	void ProteinGroup::Score()
	{
		// sum the scores of the best PSM per base sequence
		setProteinGroupScore(getAllPsmsBelowOnePercentFDR().GroupBy([&] (std::any p)
		{
			p::BaseSequence;
		})->Select([&] (std::any p)
		{
			p->Select([&] (std::any x)
			{
				x::Score;
			}).Max();
		}).Sum());
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
			if (psm->BaseSequence != nullptr)
			{
				auto peptides = psm->BestMatchingPeptides->Select([&] (std::any p)
				{
					p::Peptide;
				});
				for (auto peptide : peptides)
				{
					// might be unambiguous but also shared; make sure this protein group contains this peptide+protein combo
					if (std::find(Proteins.begin(), Proteins.end(), peptide->Protein) != Proteins.end())
					{
						proteinsWithUnambigSeqPsms[peptide->Protein].push_back(peptide);

						// null FullSequence means that mods were not successfully localized; do not display them on the sequence coverage mods info
						if (psm->FullSequence != nullptr)
						{
							proteinsWithPsmsWithLocalizedMods[peptide->Protein].push_back(peptide);
						}
					}
				}
			}
		}

		for (auto protein : ListOfProteinsOrderedByAccession)
		{
			bool errorResult = false;
			auto sequenceCoverageDisplay = protein->BaseSequence->ToLower(CultureInfo::InvariantCulture);
			std::unordered_set<int> coveredOneBasedResidues;

			// get residue numbers of each peptide in the protein and identify them as observed if the sequence is unambiguous
			for (auto peptide : proteinsWithUnambigSeqPsms[protein])
			{
				std::wstring sequenceExtractedFromProtein = L"";
				for (int i = peptide->OneBasedStartResidueInProtein; i <= peptide->OneBasedEndResidueInProtein; i++)
				{
					// check for bugs in sequence coverage; make sure we have the right amino acids!
					sequenceExtractedFromProtein += sequenceCoverageDisplay[i - 1];
					coveredOneBasedResidues.insert(i);
				}

				if (StringHelper::toUpper(sequenceExtractedFromProtein) != peptide->BaseSequence)
				{
					errorResult = true;
				}
			}

			// calculate sequence coverage percent
			double seqCoveragePercent = static_cast<double>(coveredOneBasedResidues.size()) / protein->Length;
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
			auto coverageArray = sequenceCoverageDisplay->ToCharArray();
			for (auto obsResidueLocation : coveredOneBasedResidues)
			{
				coverageArray[obsResidueLocation - 1] = std::toupper(coverageArray[obsResidueLocation - 1]);
			}
			sequenceCoverageDisplay = std::wstring(coverageArray);

			// check to see if there was an errored result; if not, add the coverage display
			if (!errorResult)
			{
				getSequenceCoverageDisplayList().push_back(sequenceCoverageDisplay);
			}
			else
			{
				getSequenceCoverageDisplayList().push_back(L"Error calculating sequence coverage");
			}

			// put mods in the sequence coverage display
			if (!errorResult)
			{
				// get mods to display in sequence (only unambiguously identified mods)
				auto modsOnThisProtein = std::unordered_set<KeyValuePair<int, Modification*>*>();
				for (auto pep : proteinsWithPsmsWithLocalizedMods[protein])
				{
					for (auto mod : pep->AllModsOneIsNterminus)
					{
						if (!mod->Value->ModificationType->Contains(L"PeptideTermMod") && !mod->Value->ModificationType->Contains(L"Common Variable") && !mod->Value->ModificationType->Contains(L"Common Fixed"))
						{
							modsOnThisProtein.insert(KeyValuePair<int, Modification*>(pep->OneBasedStartResidueInProtein + mod->Key - 2, mod->Value));
						}
					}
				}

				auto temp1 = modsOnThisProtein.OrderBy([&] (std::any p)
				{
					p::Key;
				}).ToList();

				for (auto mod : temp1)
				{
					if (mod.Value->LocationRestriction->Equals(L"N-terminal."))
					{
						sequenceCoverageDisplay = sequenceCoverageDisplay->Insert(0, L"[" + mod.Value->IdWithMotif + L"]-");
					}
					else if (mod.Value->LocationRestriction->Equals(L"Anywhere."))
					{
						int modStringIndex = sequenceCoverageDisplay->Length - (protein->Length - mod.Key);
						sequenceCoverageDisplay = sequenceCoverageDisplay->Insert(modStringIndex, L"[" + mod.Value->IdWithMotif + L"]");
					}
					else if (mod.Value->LocationRestriction->Equals(L"C-terminal."))
					{
						sequenceCoverageDisplay = sequenceCoverageDisplay->Insert(sequenceCoverageDisplay->Length, L"-[" + mod.Value->IdWithMotif + L"]");
					}
				}

				getSequenceCoverageDisplayListWithMods().push_back(sequenceCoverageDisplay);

				if (modsOnThisProtein.Any())
				{
					// calculate spectral count percentage of modified observation
					std::wstring tempModStrings = L""; //The whole string
					std::vector<int> tempPepModTotals; //The List of (For one mod, The Modified Pep Num)
					std::vector<int> tempPepTotals; //The List of (For one mod, The total Pep Num)
					std::vector<std::wstring> tempPepModValues; //The List of (For one mod, the Modified Name)
					std::vector<int> tempModIndex; //The Index of the modified position.

					for (auto pep : proteinsWithPsmsWithLocalizedMods[protein])
					{
						for (auto mod : pep->AllModsOneIsNterminus)
						{
							int tempPepNumTotal = 0; //For one mod, The total Pep Num
							if (!mod->Value->ModificationType->Contains(L"Common Variable") && !mod->Value->ModificationType->Contains(L"Common Fixed") && !mod->Value->LocationRestriction->Equals(ModLocationOnPeptideOrProtein::PepC) && !mod->Value->LocationRestriction->Equals(ModLocationOnPeptideOrProtein::NPep))
							{
								int tempIndexInProtein;
								if (mod->Value->LocationRestriction->Equals(L"N-terminal."))
								{
									tempIndexInProtein = 1;
								}
								else if (mod->Value->LocationRestriction->Equals(L"Anywhere."))
								{
									tempIndexInProtein = pep->OneBasedStartResidueInProtein + mod->Key - 2;
								}
								else if (mod->Value->LocationRestriction->Equals(L"C-terminal."))
								{
									tempIndexInProtein = protein->Length;
								}
								else
								{
									// In case it's a peptide terminal mod, skip!
									// we don't want this annotated in the protein's modifications
									continue;
								}

								if (std::find(tempModIndex.begin(), tempModIndex.end(), tempIndexInProtein) != tempModIndex.end() && tempPepModValues[tempModIndex.find(tempIndexInProtein)] == mod->Value->IdWithMotif)
								{
									tempPepModTotals[tempModIndex.find(tempIndexInProtein)] += 1;
								}
								else
								{
									tempModIndex.push_back(tempIndexInProtein);
									for (auto pept : proteinsWithPsmsWithLocalizedMods[protein])
									{
										if (tempIndexInProtein >= pept->OneBasedStartResidueInProtein - (tempIndexInProtein == 1 ? 1 : 0) && tempIndexInProtein <= pept->OneBasedEndResidueInProtein)
										{
											tempPepNumTotal += 1;
										}
									}
									tempPepTotals.push_back(tempPepNumTotal);
									tempPepModValues.push_back(mod->Value->IdWithMotif);
									tempPepModTotals.push_back(1);
								}
							}
						}
					}
					for (int i = 0; i < tempPepModTotals.size(); i++)
					{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
						std::wstring tempString = (L"#aa" + std::to_wstring(tempModIndex[i]) + L"[" + tempPepModValues[i].ToString() + L",info:occupancy=" + (static_cast<double>(tempPepModTotals[i]) / static_cast<double>(tempPepTotals[i])).ToString(L"F2") + L"(" + std::to_wstring(tempPepModTotals[i]) + L"/" + std::to_wstring(tempPepTotals[i]) + L")" + L"];");
						tempModStrings += tempString;
					}

					if (!tempModStrings.empty())
					{
						getModsInfo().push_back(tempModStrings);
					}
				}
			}

//C# TO C++ CONVERTER TODO TASK: A 'delete sequenceCoverageDisplay' statement was not added since sequenceCoverageDisplay was passed to a method or constructor. Handle memory management manually.
		}
	}

	void ProteinGroup::MergeProteinGroupWith(ProteinGroup *other)
	{
		this->getProteins().UnionWith(other->getProteins());
		this->getAllPeptides().UnionWith(other->getAllPeptides());
		this->getUniquePeptides().UnionWith(other->getUniquePeptides());
		this->getAllPsmsBelowOnePercentFDR().UnionWith(other->getAllPsmsBelowOnePercentFDR());
		other->setProteinGroupScore(0);

		ListOfProteinsOrderedByAccession = getProteins().OrderBy([&] (std::any p)
		{
			p::Accession;
		}).ToList();

		setProteinGroupName(std::wstring::Join(L"|", ListOfProteinsOrderedByAccession.Select([&] (std::any p)
		{
			p::Accession;
		})));
	}

	ProteinGroup *ProteinGroup::ConstructSubsetProteinGroup(const std::wstring &fullFilePath)
	{
		auto allPsmsForThisFile = std::unordered_set<PeptideSpectralMatch*>(this->getAllPsmsBelowOnePercentFDR().Where([&] (std::any p)
		{
			p::FullFilePath->Equals(fullFilePath);
		}));
		auto allPeptidesForThisFile = std::unordered_set<PeptideWithSetModifications*>(allPsmsForThisFile.SelectMany([&] (std::any p)
		{
			p::BestMatchingPeptides->Select([&] (std::any v)
			{
				v::Peptide;
			});
		}));
		auto allUniquePeptidesForThisFile = std::unordered_set<PeptideWithSetModifications*>(this->getUniquePeptides().Intersect(allPeptidesForThisFile));

		ProteinGroup *subsetPg = new ProteinGroup(this->getProteins(), allPeptidesForThisFile, allUniquePeptidesForThisFile);
		subsetPg->setAllPsmsBelowOnePercentFDR(allPsmsForThisFile);
		subsetPg->setDisplayModsOnPeptides(this->getDisplayModsOnPeptides());

		SpectraFileInfo *spectraFileInfo = nullptr;
		if (getFilesForQuantification().size() > 0)
		{
			spectraFileInfo = getFilesForQuantification().Where([&] (std::any p)
			{
			delete subsetPg;
				return p->FullFilePathWithExtension == fullFilePath;
			}).First();
			subsetPg->setFilesForQuantification(std::vector<SpectraFileInfo*> {spectraFileInfo});
		}

		if (getIntensitiesByFile().empty())
		{
			subsetPg->getIntensitiesByFile().clear();
		}
		else
		{
			subsetPg->setIntensitiesByFile(std::unordered_map<SpectraFileInfo*, double>
			{
				{spectraFileInfo, getIntensitiesByFile()[spectraFileInfo]}
			});
		}

//C# TO C++ CONVERTER TODO TASK: A 'delete subsetPg' statement was not added since subsetPg was used in a 'return' or 'throw' statement.
		return subsetPg;
	}
}
