#include "GlobalSettings.h"

namespace EngineLayer
{

	bool GlobalSettings::getWriteExcelCompatibleTSVs() const
	{
		return privateWriteExcelCompatibleTSVs;
	}

	void GlobalSettings::setWriteExcelCompatibleTSVs(bool value)
	{
		privateWriteExcelCompatibleTSVs = value;
	}
}
