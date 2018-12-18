#pragma once

#include "IGlobalSettings.h"

namespace EngineLayer
{
	class GlobalSettings : public IGlobalSettings
	{
	private:
		bool privateWriteExcelCompatibleTSVs = false;

	public:
		bool getWriteExcelCompatibleTSVs() const override;
		void setWriteExcelCompatibleTSVs(bool value) override;
	};
}
