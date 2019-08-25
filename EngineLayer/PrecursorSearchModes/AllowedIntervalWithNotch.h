#pragma once

using namespace MzLibUtil;

namespace EngineLayer
{
	class AllowedIntervalWithNotch
	{
	private:
		int privateNotch = 0;

	public:
		DoubleRange *AllowedInterval;

		virtual ~AllowedIntervalWithNotch()
		{
			delete AllowedInterval;
		}

		AllowedIntervalWithNotch(DoubleRange *doubleRange, int j);

		int getNotch() const;
	};
}
