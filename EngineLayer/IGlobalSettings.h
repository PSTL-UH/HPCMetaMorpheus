#pragma once

namespace EngineLayer
{
	class IGlobalSettings
	{
	public:
		virtual bool getWriteExcelCompatibleTSVs() const = 0;
	};
}
