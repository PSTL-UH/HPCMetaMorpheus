#include "SeqCoverageTest.h"
#include "../EngineLayer/IScan.h"
#include "ModificationAnalysisTest.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/ProteinParsimony/ProteinParsimonyEngine.h"
#include "../EngineLayer/ProteinParsimony/ProteinParsimonyResults.h"
#include "../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.h"

using namespace EngineLayer;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace MassSpectrometry;

namespace Test
{

	void SeqCoverageTest::TryFailSequenceCoverage()
	{
		auto prot1 = new Protein(L"MMKMMK", L"prot1");

		ModificationMotif motifM;
		ModificationMotif::TryGetMotif(L"M", motifM);
		Modification *mod1 = new Modification(_originalId: L"mod1", _modificationType: L"mt", _target: motifM, _locationRestriction: L"N-terminal.", _monoisotopicMass: 10);
		Modification *mod2 = new Modification(_originalId: L"mod2", _modificationType: L"mt", _target: motifM, _locationRestriction: L"Peptide N-terminal.", _monoisotopicMass: 10);
		Modification *mod3 = new Modification(_originalId: L"mod3", _modificationType: L"mt", _target: motifM, _locationRestriction: L"Anywhere.", _monoisotopicMass: 10);
		ModificationMotif motifK;
		ModificationMotif::TryGetMotif(L"K", motifK);
		Modification *mod4 = new Modification(_originalId: L"mod4", _modificationType: L"mt", _target: motifK, _locationRestriction: L"Peptide C-terminal.", _monoisotopicMass: 10);
		Modification *mod5 = new Modification(_originalId: L"mod5", _modificationType: L"mt", _target: motifK, _locationRestriction: L"C-terminal.", _monoisotopicMass: 10);

		std::unordered_map<int, Modification*> modsFor1 =
		{
			{1, mod1},
			{3, mod3},
			{5, mod4}
		};
		std::unordered_map<int, Modification*> modsFor2 =
		{
			{1, mod2},
			{5, mod5}
		};
		std::unordered_map<int, Modification*> modsFor3 =
		{
			{1, mod1},
			{5, mod3},
			{8, mod5}
		};

		DigestionParams *digestionParams = new DigestionParams();
		auto pwsm1 = new PeptideWithSetModifications(prot1, digestionParams, 1, 3, CleavageSpecificity::Unknown, L"", 0, modsFor1, 0);
		auto pwsm2 = new PeptideWithSetModifications(prot1, digestionParams, 4, 6, CleavageSpecificity::Unknown, L"", 0, modsFor2, 0);
		auto pwsm3 = new PeptideWithSetModifications(prot1, digestionParams, 1, 6, CleavageSpecificity::Unknown, L"", 0, modsFor3, 0);

		std::unordered_set<PeptideWithSetModifications*> peptides = {pwsm1, pwsm2, pwsm3};

		IScan *scan = new ThisTestScan();
		auto psm1 = new PeptideSpectralMatch(pwsm1, 0, 1, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		psm1->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);

		auto psm2 = new PeptideSpectralMatch(pwsm2, 0, 1, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		psm2->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);

		auto psm3 = new PeptideSpectralMatch(pwsm3, 0, 1, 0, scan, digestionParams, std::vector<MatchedFragmentIon*>());
		psm3->SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);


		std::vector<PeptideSpectralMatch*> newPsms = {psm1, psm2, psm3};

		std::for_each(newPsms.begin(), newPsms.end(), [&] (std::any p)
		{
			p::ResolveAllAmbiguities();
		});

		CommonParameters tempVar();
		ProteinParsimonyEngine *ppe = new ProteinParsimonyEngine(newPsms, true, &tempVar, std::vector<std::wstring>());
		ProteinParsimonyResults *fjkd = static_cast<ProteinParsimonyResults*>(ppe->Run());

		CommonParameters tempVar2();
		ProteinScoringAndFdrEngine *psafe = new ProteinScoringAndFdrEngine(fjkd->getProteinGroups(), newPsms, true, true, true, &tempVar2, std::vector<std::wstring>());

		psafe->Run();

		fjkd->getProteinGroups().front().CalculateSequenceCoverage();

		auto firstSequenceCoverageDisplayList = fjkd->getProteinGroups().front().SequenceCoverageDisplayList.First();
		Assert::AreEqual(L"MMKMMK", firstSequenceCoverageDisplayList);
		auto firstSequenceCoverageDisplayListWithMods = fjkd->getProteinGroups().front().SequenceCoverageDisplayListWithMods.First();
		Assert::AreEqual(L"[mod1 on M]-MM[mod3 on M]KM[mod3 on M]MK-[mod5 on K]", firstSequenceCoverageDisplayListWithMods);

		auto firstModInfo = fjkd->getProteinGroups().front().ModsInfo.First();
		Assert::IsTrue(firstModInfo->Contains(LR"(#aa1[mod1 on M,info:occupancy=1.00(2/2)])"));
		Assert::IsTrue(firstModInfo->Contains(LR"(#aa2[mod3 on M,info:occupancy=0.50(1/2)])"));
		Assert::IsFalse(firstModInfo->Contains(LR"(#aa3)"));
		Assert::IsTrue(firstModInfo->Contains(LR"(#aa4[mod3 on M,info:occupancy=0.50(1/2)])"));
		Assert::IsFalse(firstModInfo->Contains(LR"(#aa5)"));
		Assert::IsTrue(firstModInfo->Contains(LR"(#aa6[mod5 on K,info:occupancy=1.00(2/2)])"));

		delete psafe;
		delete ppe;
		delete psm3;
		delete psm2;
		delete psm1;
//C# TO C++ CONVERTER TODO TASK: A 'delete scan' statement was not added since scan was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pwsm3' statement was not added since pwsm3 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pwsm2' statement was not added since pwsm2 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete pwsm1' statement was not added since pwsm1 was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete digestionParams' statement was not added since digestionParams was passed to a method or constructor. Handle memory management manually.
		delete mod5;
		delete mod4;
		delete mod3;
		delete mod2;
		delete mod1;
//C# TO C++ CONVERTER TODO TASK: A 'delete prot1' statement was not added since prot1 was passed to a method or constructor. Handle memory management manually.
	}
}
