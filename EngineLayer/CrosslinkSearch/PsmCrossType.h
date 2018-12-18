#pragma once

namespace EngineLayer
{
	namespace CrosslinkSearch
	{
		enum class PsmCrossType
		{
			Single,
			Cross,
			DeadEnd,
			Loop,
			Inter,
			Intra,
			DeadEndH2O,
			DeadEndNH2,
			DeadEndTris
		};
	}
}
