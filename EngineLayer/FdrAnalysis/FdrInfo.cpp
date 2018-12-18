#include "FdrInfo.h"

namespace EngineLayer
{
	namespace FdrAnalysis
	{

		double FdrInfo::getCumulativeTarget() const
		{
			return privateCumulativeTarget;
		}

		void FdrInfo::setCumulativeTarget(double value)
		{
			privateCumulativeTarget = value;
		}

		double FdrInfo::getCumulativeDecoy() const
		{
			return privateCumulativeDecoy;
		}

		void FdrInfo::setCumulativeDecoy(double value)
		{
			privateCumulativeDecoy = value;
		}

		double FdrInfo::getCumulativeTargetNotch() const
		{
			return privateCumulativeTargetNotch;
		}

		void FdrInfo::setCumulativeTargetNotch(double value)
		{
			privateCumulativeTargetNotch = value;
		}

		double FdrInfo::getCumulativeDecoyNotch() const
		{
			return privateCumulativeDecoyNotch;
		}

		void FdrInfo::setCumulativeDecoyNotch(double value)
		{
			privateCumulativeDecoyNotch = value;
		}

		double FdrInfo::getQValue() const
		{
			return privateQValue;
		}

		void FdrInfo::setQValue(double value)
		{
			privateQValue = value;
		}

		double FdrInfo::getQValueNotch() const
		{
			return privateQValueNotch;
		}

		void FdrInfo::setQValueNotch(double value)
		{
			privateQValueNotch = value;
		}

		bool FdrInfo::getCalculateEValue() const
		{
			return privateCalculateEValue;
		}

		void FdrInfo::setCalculateEValue(bool value)
		{
			privateCalculateEValue = value;
		}

		double FdrInfo::getMaximumLikelihood() const
		{
			return privateMaximumLikelihood;
		}

		void FdrInfo::setMaximumLikelihood(double value)
		{
			privateMaximumLikelihood = value;
		}

		double FdrInfo::getEValue() const
		{
			return privateEValue;
		}

		void FdrInfo::setEValue(double value)
		{
			privateEValue = value;
		}

		double FdrInfo::getEScore() const
		{
			return privateEScore;
		}

		void FdrInfo::setEScore(double value)
		{
			privateEScore = value;
		}
	}
}
