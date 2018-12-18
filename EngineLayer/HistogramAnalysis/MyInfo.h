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
			std::wstring infostring;

			MyInfo(double MassShift, const std::wstring &infostring);

				double getMassShift() const;
				void setMassShift(double value);
		};
	}
}
