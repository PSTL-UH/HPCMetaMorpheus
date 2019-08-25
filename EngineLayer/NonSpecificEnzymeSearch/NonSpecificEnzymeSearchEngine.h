#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class Ms2ScanWithSpecificMass; }
namespace EngineLayer { class MassDiffAcceptor; }
namespace EngineLayer { class MetaMorpheusEngineResults; }

using namespace Chemistry;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::ModernSearch;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace NonSpecificEnzymeSearch
	{
		class NonSpecificEnzymeSearchEngine : public ModernSearchEngine
		{
		private:
			static const double WaterMonoisotopicMass;

			std::vector<std::vector<int>> const PrecursorIndex;
			const int MinimumPeptideLength;
			std::vector<std::vector<PeptideSpectralMatch*>> GlobalCategorySpecificPsms;
			CommonParameters *ModifiedParametersNoComp;
			std::vector<ProductType*> ProductTypesToSearch;

		public:
			virtual ~NonSpecificEnzymeSearchEngine()
			{
				delete ModifiedParametersNoComp;
			}

			NonSpecificEnzymeSearchEngine(std::vector<std::vector<PeptideSpectralMatch*>> &globalPsms, std::vector<Ms2ScanWithSpecificMass*> &listOfSortedms2Scans, std::vector<PeptideWithSetModifications*> &peptideIndex, std::vector<std::vector<int>&> &fragmentIndex, std::vector<std::vector<int>&> &precursorIndex, int currentPartition, CommonParameters *CommonParameters, MassDiffAcceptor *massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled, std::vector<std::string> &nestedIds);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;

		private:
			std::tuple<int, PeptideWithSetModifications*> Accepts(std::vector<Product*> &fragments, double scanPrecursorMass, PeptideWithSetModifications *peptide, FragmentationTerminus *fragmentationTerminus, MassDiffAcceptor *searchMode);

		public:
			static std::vector<PeptideSpectralMatch*> ResolveFdrCategorySpecificPsms(std::vector<std::vector<PeptideSpectralMatch*>&> &AllPsms, int numNotches, const std::string &taskId, CommonParameters *commonParameters);
		};
	}
}
