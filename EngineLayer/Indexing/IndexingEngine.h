﻿#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>
#include <cmath>
#include <mutex>
#include "exceptionhelper.h"
#include "stringhelper.h"
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }

using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace EngineLayer
{
	namespace Indexing
	{
		class IndexingEngine : public MetaMorpheusEngine
		{
		private:
			static constexpr int FragmentBinsPerDalton = 1000;
			const std::vector<Protein*> ProteinList;
			const std::vector<Modification*> FixedModifications;
			const std::vector<Modification*> VariableModifications;
			const int CurrentPartition;
			DecoyType *const DecoyType;
			const double MaxFragmentSize;
		public:
			const bool GeneratePrecursorIndex;
			const std::vector<FileInfo*> ProteinDatabases;

			virtual ~IndexingEngine()
			{
				delete DecoyType;
			}

			IndexingEngine(std::vector<Protein*> &proteinList, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, int currentPartition, DecoyType *decoyType, CommonParameters *commonParams, double maxFragmentSize, bool generatePrecursorIndex, std::vector<FileInfo*> &proteinDatabases, std::vector<std::wstring> &nestedIds);

			std::wstring ToString() override;

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;
		};
	}
}