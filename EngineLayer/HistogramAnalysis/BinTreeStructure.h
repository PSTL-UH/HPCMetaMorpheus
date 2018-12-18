#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <limits>
#include <tuple>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { namespace HistogramAnalysis { class Bin; } }
namespace EngineLayer { class PeptideSpectralMatch; }

using namespace MassSpectrometry;
using namespace MathNet::Numerics::Statistics;
using namespace Proteomics::Fragmentation;

namespace EngineLayer
{
	namespace HistogramAnalysis
	{
		class BinTreeStructure
		{
		private:
			std::vector<Bin*> privateFinalBins;

			static constexpr int MinAdditionalPsmsInBin = 1;

			public:
				std::vector<Bin*> getFinalBins() const;
				void setFinalBins(const std::vector<Bin*> &value);

			void GenerateBins(std::vector<PeptideSpectralMatch*> &targetAndDecoyMatches, double dc);

		private:
			static double GetSigma(double thisMassShift, int thisP, int i, std::vector<double> &listOfMassShifts, std::vector<int> &p);

			void IdentifyFracWithSingle();

			void IdentifyMedianLength();

			void IdentifyAAsInCommon();

			void IdentifyMods();

			void IdentifyResidues();

			void IdentifyUnimodBins(double v);

			void IdentifyUniprotBins(double v);

			void IdentifyCombos(double v);

			void IdentifyAA(double v);

			void IdentifyMine(double v);
		};
	}
}
