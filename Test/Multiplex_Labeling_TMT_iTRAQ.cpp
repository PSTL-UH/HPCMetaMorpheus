#include "Multiplex_Labeling_TMT_iTRAQ.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/CommonParameters.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace Test
{

	void Multiplex_Labeling_TMT_iTRAQ::TestChemicalFormulaWithIsotopesTMT(const std::wstring &formula, double mass)
	{
		ChemicalFormula *cf = ChemicalFormula::ParseFormula(formula);
		Assert::AreEqual(mass, ClassExtensions::RoundedDouble(cf->MonoisotopicMass));
	}

	void Multiplex_Labeling_TMT_iTRAQ::TestPeptideLabelledWithTMT(const std::wstring &peptide, double totalMass)
	{
		std::vector<Modification*> gptmdModifications;
		gptmdModifications.insert(gptmdModifications.end(), GlobalVariables::getAllModsKnown().begin(), GlobalVariables::getAllModsKnown().end());
		std::vector<Modification*> tmt10Mods = gptmdModifications.Where([&] (std::any m)
		{
			return m->ModificationType == L"Multiplex Label" && m::IdWithMotif->Contains(L"TMT10");
		}).ToList();

		Protein *P = new Protein(peptide, L"", L"", nullptr, nullptr, nullptr, nullptr, nullptr, false, false, nullptr, nullptr, nullptr, nullptr);
		DigestionParams tempVar(minPeptideLength: 1);
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar);
		auto p = P->Digest(CommonParameters->getDigestionParams(), tmt10Mods, std::vector<Modification*>()).First();
		auto f = p->Fragment(DissociationType::HCD, FragmentationTerminus::Both);

		std::vector<double> productMasses = f->Select([&] (std::any m)
		{
			m::NeutralMass::ToMz(1);
		}).ToList();
		productMasses.Distinct();
		std::sort(productMasses.begin(), productMasses.end());

		Assert::AreEqual(totalMass, ClassExtensions::RoundedDouble(p->MonoisotopicMass.ToMz(1), 4));

		delete CommonParameters;
		delete P;
	}

	void Multiplex_Labeling_TMT_iTRAQ::TestPeptideLabelledWith_iTRAQ_4plex(const std::wstring &peptide, double totalMass)
	{
		std::vector<Modification*> gptmdModifications;
		gptmdModifications.insert(gptmdModifications.end(), GlobalVariables::getAllModsKnown().begin(), GlobalVariables::getAllModsKnown().end());
		std::vector<Modification*> itraq4plex = gptmdModifications.Where([&] (std::any m)
		{
			return m->ModificationType == L"Multiplex Label" && m::IdWithMotif->Contains(L"iTRAQ-4plex");
		}).ToList();

		Protein *P = new Protein(peptide, L"", L"", nullptr, nullptr, nullptr, nullptr, nullptr, false, false, nullptr, nullptr, nullptr, nullptr);
		DigestionParams tempVar(minPeptideLength: 1);
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar);
		auto p = P->Digest(CommonParameters->getDigestionParams(), itraq4plex, std::vector<Modification*>()).First();
		auto f = p->Fragment(DissociationType::HCD, FragmentationTerminus::Both);

		std::vector<double> productMasses = f->Select([&] (std::any m)
		{
			m::NeutralMass::ToMz(1);
		}).ToList();
		productMasses.Distinct();
		std::sort(productMasses.begin(), productMasses.end());

		Assert::AreEqual(totalMass, ClassExtensions::RoundedDouble(p->MonoisotopicMass.ToMz(1), 4));

		delete CommonParameters;
		delete P;
	}

	void Multiplex_Labeling_TMT_iTRAQ::TestPeptideLabelledWith_DiLeu_4plex(const std::wstring &peptide, double totalMass)
	{
		std::vector<Modification*> gptmdModifications;
		gptmdModifications.insert(gptmdModifications.end(), GlobalVariables::getAllModsKnown().begin(), GlobalVariables::getAllModsKnown().end());
		std::vector<Modification*> itraq4plex = gptmdModifications.Where([&] (std::any m)
		{
			return m->ModificationType == L"Multiplex Label" && m::IdWithMotif->Contains(L"DiLeu-4plex");
		}).ToList();

		Protein *P = new Protein(peptide, L"", L"", nullptr, nullptr, nullptr, nullptr, nullptr, false, false, nullptr, nullptr, nullptr, nullptr);
		DigestionParams tempVar(minPeptideLength: 1);
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar);
		auto p = P->Digest(CommonParameters->getDigestionParams(), itraq4plex, std::vector<Modification*>()).First();
		auto f = p->Fragment(DissociationType::HCD, FragmentationTerminus::Both);

		std::vector<double> productMasses = f->Select([&] (std::any m)
		{
			m::NeutralMass::ToMz(1);
		}).ToList();
		productMasses.Distinct();
		std::sort(productMasses.begin(), productMasses.end());

		Assert::AreEqual(totalMass, ClassExtensions::RoundedDouble(p->MonoisotopicMass.ToMz(1), 4));

		delete CommonParameters;
		delete P;
	}

	void Multiplex_Labeling_TMT_iTRAQ::TestPeptideLabelledWith_iTRAQ_8plex(const std::wstring &peptide, double totalMass)
	{
		std::vector<Modification*> gptmdModifications;
		gptmdModifications.insert(gptmdModifications.end(), GlobalVariables::getAllModsKnown().begin(), GlobalVariables::getAllModsKnown().end());
		std::vector<Modification*> itraq8plex = gptmdModifications.Where([&] (std::any m)
		{
			return m->ModificationType == L"Multiplex Label" && m::IdWithMotif->Contains(L"iTRAQ-8plex");
		}).ToList();

		Protein *P = new Protein(peptide, L"", L"", nullptr, nullptr, nullptr, nullptr, nullptr, false, false, nullptr, nullptr, nullptr, nullptr);
		DigestionParams tempVar(minPeptideLength: 1);
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar);
		auto p = P->Digest(CommonParameters->getDigestionParams(), itraq8plex, std::vector<Modification*>()).First();
		auto f = p->Fragment(DissociationType::HCD, FragmentationTerminus::Both);

		std::vector<double> productMasses = f->Select([&] (std::any m)
		{
			m::NeutralMass::ToMz(1);
		}).ToList();
		productMasses.Distinct();
		std::sort(productMasses.begin(), productMasses.end());

		Assert::AreEqual(totalMass, ClassExtensions::RoundedDouble(p->MonoisotopicMass.ToMz(1), 4));

		delete CommonParameters;
		delete P;
	}

	void Multiplex_Labeling_TMT_iTRAQ::TestChemicalFormulaWithIsotopes_iTRAQ(const std::wstring &formula, double mass, bool mz)
	{
		ChemicalFormula *cf = ChemicalFormula::ParseFormula(formula);
		if (mz)
		{
			Assert::AreEqual(mass, ClassExtensions::RoundedDouble(cf->MonoisotopicMass.ToMz(1)));
		}
		else
		{
			Assert::AreEqual(mass, ClassExtensions::RoundedDouble(cf->MonoisotopicMass));
		}
	}

	void Multiplex_Labeling_TMT_iTRAQ::TestChemicalFormulaWithIsotopes_DiLeu4plex(const std::wstring &formula, double mass, bool mz)
	{
		ChemicalFormula *cf = ChemicalFormula::ParseFormula(formula);
		if (mz)
		{
			Assert::AreEqual(mass, ClassExtensions::RoundedDouble(cf->MonoisotopicMass.ToMz(1)));
		}
		else
		{
			Assert::AreEqual(mass, ClassExtensions::RoundedDouble(cf->MonoisotopicMass));
		}
	}

	void Multiplex_Labeling_TMT_iTRAQ::TestChemicalFormulaWithIsotopes_DiLeu12plex(const std::wstring &formula, double mass, bool mz)
	{
		ChemicalFormula *cf = ChemicalFormula::ParseFormula(formula);
		if (mz)
		{
			Assert::AreEqual(mass, ClassExtensions::RoundedDouble(cf->MonoisotopicMass.ToMz(1),5));
		}
		else
		{
			Assert::AreEqual(mass, ClassExtensions::RoundedDouble(cf->MonoisotopicMass));
		}
	}
}
