#include "IndexEngineTest.h"
#include "../EngineLayer/CommonParameters.h"

using namespace EngineLayer;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace NUnit::Framework;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace Test
{

	void IndexEngineTest::TestIndexEngine()
	{
		auto proteinList = std::vector<Protein*> {new Protein(L"MNNNKQQQ", nullptr)};
		auto variableModifications = std::vector<Modification*>();
		auto fixedModifications = std::vector<Modification*>();
		auto localizeableModifications = std::vector<Modification*>();

		std::unordered_map<Modification*, unsigned short> modsDictionary;
		for (auto mod : fixedModifications)
		{
			modsDictionary.emplace(mod, 0);
		}
		int i = 1;
		for (auto mod : variableModifications)
		{
			modsDictionary.emplace(mod, static_cast<unsigned short>(i));
			i++;
		}
		for (auto mod : localizeableModifications)
		{
			modsDictionary.emplace(mod, static_cast<unsigned short>(i));
			i++;
		}

		std::vector<DigestionMotif*> motifs = {new DigestionMotif(L"K", nullptr, 1, nullptr)};
		Protease *p = new Protease(L"Custom Protease2",CleavageSpecificity::Full, nullptr, nullptr, motifs);
		ProteaseDictionary::Dictionary->Add(p->Name, p);
		DigestionParams tempVar(protease: p->Name, minPeptideLength: 1);
		CommonParameters *CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: &tempVar);

		auto engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters, 30000, false, std::vector<FileInfo*>(), std::vector<std::wstring>());

		auto results = static_cast<IndexingResults*>(engine->Run());

		Assert::AreEqual(5, results->getPeptideIndex().size());

		auto digestedList = proteinList[0]->Digest(CommonParameters->getDigestionParams(), std::vector<Modification*>(), variableModifications).ToList();

		Assert::AreEqual(5, digestedList.size());
		for (auto fdfd : digestedList)
		{
			Assert->Contains(fdfd, results->getPeptideIndex());
		}

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete p' statement was not added since p was passed to a method or constructor. Handle memory management manually.
	}

	void IndexEngineTest::TestIndexEngineWithWeirdSeq()
	{
		auto proteinList = std::vector<Protein*> {new Protein(L"MQXQ", nullptr)};
		auto variableModifications = std::vector<Modification*>();
		auto fixedModifications = std::vector<Modification*>();
		auto localizeableModifications = std::vector<Modification*>();

		std::unordered_map<Modification*, unsigned short> modsDictionary;
		for (auto mod : fixedModifications)
		{
			modsDictionary.emplace(mod, 0);
		}
		int i = 1;
		for (auto mod : variableModifications)
		{
			modsDictionary.emplace(mod, static_cast<unsigned short>(i));
			i++;
		}
		for (auto mod : localizeableModifications)
		{
			modsDictionary.emplace(mod, static_cast<unsigned short>(i));
			i++;
		}

		std::vector<DigestionMotif*> motifs = {new DigestionMotif(L"K", nullptr, 1, nullptr)};
		Protease *protease = new Protease(L"Custom Protease", CleavageSpecificity::Full, nullptr, nullptr, motifs);
		ProteaseDictionary::Dictionary->Add(protease->Name, protease);
		DigestionParams tempVar(protease: protease->Name, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain);
		CommonParameters *CommonParameters = new CommonParameters(digestionParams: &tempVar, scoreCutoff: 1);

		auto engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType::Reverse, CommonParameters, 30000, false, std::vector<FileInfo*>(), std::vector<std::wstring>());

		auto results = static_cast<IndexingResults*>(engine->Run());

		Assert::AreEqual(1, results->getPeptideIndex().size());

		Assert::IsNaN(results->getPeptideIndex()[0]->MonoisotopicMass);
		Assert::AreEqual(30000000 + 1, results->getFragmentIndex().size());

		delete engine;
//C# TO C++ CONVERTER TODO TASK: A 'delete CommonParameters' statement was not added since CommonParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete protease' statement was not added since protease was passed to a method or constructor. Handle memory management manually.
	}
}
