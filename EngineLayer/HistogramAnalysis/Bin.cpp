#include "Bin.h"
#include "../PeptideSpectralMatch.h"
#include "../GlobalVariables.h"

using namespace Proteomics;
using namespace Proteomics::AminoAcidPolymer;
namespace EngineLayer
{
	namespace HistogramAnalysis
	{

		Bin::Bin(double massShift)
		{
			this->MassShift = massShift;
			UniquePSMs = std::unordered_map<std::string, std::tuple<std::string, std::string, PeptideSpectralMatch*>>();
		}

		int Bin::getPepNlocCount() const
		{
			return privatePepNlocCount;
		}

		void Bin::setPepNlocCount(int value)
		{
			privatePepNlocCount = value;
		}

		int Bin::getPepClocCount() const
		{
			return privatePepClocCount;
		}

		void Bin::setPepClocCount(int value)
		{
			privatePepClocCount = value;
		}

		int Bin::getProtNlocCount() const
		{
			return privateProtNlocCount;
		}

		void Bin::setProtNlocCount(int value)
		{
			privateProtNlocCount = value;
		}

		int Bin::getProtClocCount() const
		{
			return privateProtClocCount;
		}

		void Bin::setProtClocCount(int value)
		{
			privateProtClocCount = value;
		}

		std::string Bin::getCombos() const
		{
			return privateCombos;
		}

		void Bin::setCombos(const std::string &value)
		{
			privateCombos = value;
		}

		std::string Bin::getUnimodDiffs() const
		{
			return privateUnimodDiffs;
		}

		void Bin::setUnimodDiffs(const std::string &value)
		{
			privateUnimodDiffs = value;
		}

		std::string Bin::getUniprotID() const
		{
			return privateUniprotID;
		}

		void Bin::setUniprotID(const std::string &value)
		{
			privateUniprotID = value;
		}

		std::string Bin::getUnimodFormulas() const
		{
			return privateUnimodFormulas;
		}

		void Bin::setUnimodFormulas(const std::string &value)
		{
			privateUnimodFormulas = value;
		}

		std::string Bin::getUnimodId() const
		{
			return privateUnimodId;
		}

		void Bin::setUnimodId(const std::string &value)
		{
			privateUnimodId = value;
		}

		double Bin::getMassShift() const
		{
			return privateMassShift;
		}

		int Bin::getCount() const
		{
			return UniquePSMs.size();
		}

		int Bin::getCountDecoy() const
		{
			return UniquePSMs.Values->Count([&] (std::any b)
			{
				b::Item3->IsDecoy;
			});
		}

		int Bin::getCountTarget() const
		{
			return getCount() - getCountDecoy();
		}

		int Bin::getLocalizeableTarget() const
		{
			return UniquePSMs.Values->Where([&] (std::any b)
			{
				return b::Item3->LocalizedScores != nullptr;
			})->Count([&] (std::any b)
			{
				return !b::Item3->IsDecoy && b::Item3->LocalizedScores.Max() >= b::Item3->Score + 1;
			});
		}

		std::string Bin::getMine() const
		{
			return privateMine;
		}

		void Bin::setMine(const std::string &value)
		{
			privateMine = value;
		}

		std::unordered_map<char, int> Bin::getAAsInCommon() const
		{
			return privateAAsInCommon;
		}

		void Bin::setAAsInCommon(const std::unordered_map<char, int> &value)
		{
			privateAAsInCommon = value;
		}

		int Bin::getOverlapping() const
		{
			return privateOverlapping;
		}

		void Bin::setOverlapping(int value)
		{
			privateOverlapping = value;
		}

		double Bin::getFracWithSingle() const
		{
			return privateFracWithSingle;
		}

		void Bin::setFracWithSingle(double value)
		{
			privateFracWithSingle = value;
		}

		double Bin::getMedianLength() const
		{
			return privateMedianLength;
		}

		void Bin::setMedianLength(double value)
		{
			privateMedianLength = value;
		}

		void Bin::IdentifyResidues()
		{
			ResidueCount = std::unordered_map<char, int>();
			for (auto hehe : UniquePSMs.Values->Where([&] (std::any b)
			{
				return b::Item3->LocalizedScores != nullptr;
			}))
			{
				double bestScore = hehe::Item3->LocalizedScores.Max();
				if (bestScore >= hehe::Item3->Score + 1 && !hehe::Item3->IsDecoy)
				{
					for (int i = 0; i < hehe::Item1->Count(); i++)
					{
						if (bestScore - hehe::Item3->LocalizedScores[i] < 0.5)
						{
							if (ResidueCount.find(hehe::Item1[i]) != ResidueCount.end())
							{
								ResidueCount[hehe::Item1[i]]++;
							}
							else
							{
								ResidueCount.emplace(hehe::Item1[i], 1);
							}
						}
					}
					if (hehe::Item3->LocalizedScores.Max() - hehe::Item3->LocalizedScores[0] < 0.5)
					{
						setPepNlocCount(getPepNlocCount() + 1);
						if (hehe::Item3->OneBasedStartResidueInProtein.HasValue && hehe::Item3->OneBasedStartResidueInProtein->Value <= 2)
						{
							setProtNlocCount(getProtNlocCount() + 1);
						}
					}
					if (hehe::Item3->LocalizedScores.Max() - hehe::Item3->LocalizedScores.Last() < 0.5)
					{
						setPepClocCount(getPepClocCount() + 1);
						if (hehe::Item3->OneBasedEndResidueInProtein.HasValue && hehe::Item3->ProteinLength.HasValue && hehe::Item3->OneBasedEndResidueInProtein->Value == hehe::Item3->ProteinLength->Value)
						{
							setProtClocCount(getProtClocCount() + 1);
						}
					}
				}
			}
		}

		void Bin::IdentifyCombos(double v, std::unordered_set<std::tuple<double, double, double>> &ok)
		{
			auto okk = std::unordered_set<std::string>();
			for (auto hm : ok)
			{
				if (std::abs(hm->Item1 + hm->Item2 - getMassShift()) <= v && getCountTarget() < hm->Item3)
				{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					okk.insert("Combo " + std::min(hm->Item1, hm->Item2).ToString("F3", CultureInfo::InvariantCulture) + " and " + std::max(hm->Item1, hm->Item2).ToString("F3", CultureInfo::InvariantCulture));
				}
			}
			setCombos(std::string::Join("|", okk));
		}

		double Bin::ComputeZ(double v)
		{
			return std::sqrt(getCount()) * (static_cast<double>(getCountDecoy()) / getCount() - v) / (v * (1 - v));
		}

		void Bin::IdentifyUniprotBins(double v)
		{
			auto modIdOptions = std::unordered_set<std::string>();
			for (auto mod : GlobalVariables::getUniprotDeseralized())
			{
				if (mod->MonoisotopicMass.HasValue && std::abs(mod->MonoisotopicMass->Value - getMassShift()) <= v)
				{
					modIdOptions.insert(mod->IdWithMotif);
				}
			}
			setUniprotID(std::string::Join("|", modIdOptions));
		}

		void Bin::IdentifyAA(double v)
		{
			auto ok = std::unordered_set<std::string>();
			for (char c = 'A'; c <= 'Z'; c++)
			{
				Residue residue;
				if (Residue::TryGetResidue(c, residue))
				{
					if (std::abs(residue::MonoisotopicMass - getMassShift()) <= v)
					{
						ok.insert("Add " + residue->Name);
					}
					if (std::abs(residue::MonoisotopicMass + getMassShift()) <= v)
					{
						ok.insert("Remove " + residue->Name);
					}
					for (char cc = 'A'; cc <= 'Z'; cc++)
					{
						Residue residueCC;
						if (Residue::TryGetResidue(cc, residueCC))
						{
							if (std::abs(residueCC::MonoisotopicMass + residue::MonoisotopicMass - getMassShift()) <= v)
							{
								ok.insert("Add (" + residue->Name + "+" + residueCC->Name + ")");
							}
							if (std::abs(residueCC::MonoisotopicMass + residue::MonoisotopicMass + getMassShift()) <= v)
							{
								ok.insert("Remove (" + residue->Name + "+" + residueCC->Name + ")");
							}
						}
					}
				}
			}
			AA = std::string::Join("|", ok);
		}

		void Bin::IdentifyUnimodBins(double v)
		{
			auto ok = std::unordered_set<std::string>();
			auto okformula = std::unordered_set<std::string>();
			auto okDiff = std::unordered_set<double>();
			for (auto hm : GlobalVariables::getUnimodDeserialized())
			{
				auto theMod = dynamic_cast<Modification*>(hm);
				if (std::abs(theMod->MonoisotopicMass->Value - getMassShift()) <= v)
				{
					ok.insert(hm->IdWithMotif);
					okformula.insert(theMod->ChemicalFormula.Formula);
					okDiff.insert(theMod->MonoisotopicMass->Value - getMassShift());
				}
			}
			setUnimodId(std::string::Join("|", ok));
			setUnimodFormulas(std::string::Join("|", okformula));
			setUnimodDiffs(std::string::Join("|", okDiff));
		}

		void Bin::Add(PeptideSpectralMatch *ok)
		{
			if (ok->getFullSequence() != "")
			{
				if (UniquePSMs.find(ok->getFullSequence()) != UniquePSMs.end())
				{
					auto current = UniquePSMs[ok->getFullSequence()];
					if (current.Item3->Score < ok->getScore())
					{
						UniquePSMs[ok->getFullSequence()] = std::tuple<std::string, std::string, PeptideSpectralMatch*>(ok->getBaseSequence(), ok->getFullSequence(), ok);
					}
				}
				else
				{
					UniquePSMs.emplace(ok->getFullSequence(), std::tuple<std::string, std::string, PeptideSpectralMatch*>(ok->getBaseSequence(), ok->getFullSequence(), ok));
				}
			}
		}
	}
}
