#pragma once

namespace EngineLayer
{
	namespace FdrAnalysis
	{
		class FdrInfo
		{
		private:
			double privateCumulativeTarget = 0;
			double privateCumulativeDecoy = 0;
			double privateCumulativeTargetNotch = 0;
			double privateCumulativeDecoyNotch = 0;
			double privateQValue = 0;
			double privateQValueNotch = 0;
			bool privateCalculateEValue = false;
			double privateMaximumLikelihood = 0;
			double privateEValue = 0;
			double privateEScore = 0;

			public:
				double getCumulativeTarget() const;
				void setCumulativeTarget(double value);
				double getCumulativeDecoy() const;
				void setCumulativeDecoy(double value);
				double getCumulativeTargetNotch() const;
				void setCumulativeTargetNotch(double value);
				double getCumulativeDecoyNotch() const;
				void setCumulativeDecoyNotch(double value);
				double getQValue() const;
				void setQValue(double value);
				double getQValueNotch() const;
				void setQValueNotch(double value);
				bool getCalculateEValue() const;
				void setCalculateEValue(bool value);
				double getMaximumLikelihood() const;
				void setMaximumLikelihood(double value);
				double getEValue() const;
				void setEValue(double value);
				double getEScore() const;
				void setEScore(double value);
		};
	}
}
