#include "OkBin.h"

namespace EngineLayer
{
	namespace HistogramAnalysis
	{

		OkBin::OkBin(double massShift, double sigma, int p)
		{
			MassShift = massShift;
			Sigma = sigma;
			P = p;
		}
	}
}
