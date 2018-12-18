#pragma once

namespace EngineLayer
{
	namespace HistogramAnalysis
	{
		class OkBin
		{
		public:
			double MassShift = 0;
			double Sigma = 0;
			int P = 0;

			OkBin(double massShift, double sigma, int p);
		};
	}
}
