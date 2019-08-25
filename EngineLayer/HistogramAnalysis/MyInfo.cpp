#include "MyInfo.h"

namespace EngineLayer
{
	namespace HistogramAnalysis
	{

		MyInfo::MyInfo(double MassShift, const std::string &infostring)
		{
			this->setMassShift(MassShift);
			this->infostring = infostring;
		}

		double MyInfo::getMassShift() const
		{
			return privateMassShift;
		}

		void MyInfo::setMassShift(double value)
		{
			privateMassShift = value;
		}
	}
}
