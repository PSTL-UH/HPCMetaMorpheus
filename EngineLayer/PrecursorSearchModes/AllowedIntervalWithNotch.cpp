#include "AllowedIntervalWithNotch.h"

using namespace MzLibUtil;

namespace EngineLayer
{

	AllowedIntervalWithNotch::AllowedIntervalWithNotch(DoubleRange *doubleRange, int j)
	{
		AllowedInterval = doubleRange;
		Notch = j;
	}

	int AllowedIntervalWithNotch::getNotch() const
	{
		return privateNotch;
	}
}
