#pragma once

#include <string>

namespace EngineLayer
{
	namespace HistogramAnalysis
	{
		class MyInfo
		{
		private:
			double privateMassShift = 0;

		public:
			std::string infostring;

			MyInfo(double MassShift, const std::string &infostring);

				double getMassShift() const;
				void setMassShift(double value);
		};
	}
}
