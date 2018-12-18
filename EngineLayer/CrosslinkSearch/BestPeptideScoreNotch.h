#pragma once

#include <vector>

using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace CrosslinkSearch
	{
		class BestPeptideScoreNotch
		{
		private:
			PeptideWithSetModifications *privateBestPeptide;
			double privateBestScore = 0;
			int privateBestNotch = 0;
			std::vector<int> privateTopPosition;

		public:
			BestPeptideScoreNotch(PeptideWithSetModifications *bestPeptide, double bestScore, int bestNotch);

				PeptideWithSetModifications *getBestPeptide() const;
				void setBestPeptide(PeptideWithSetModifications *value);
				double getBestScore() const;
				void setBestScore(double value);
				int getBestNotch() const;
				void setBestNotch(int value);
				std::vector<int> getTopPosition() const;
				void setTopPosition(const std::vector<int> &value);
		};
	}
}
