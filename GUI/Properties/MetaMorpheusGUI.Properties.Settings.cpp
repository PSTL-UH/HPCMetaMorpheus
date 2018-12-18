#include "MetaMorpheusGUI.Properties.Settings.h"

namespace MetaMorpheusGUI
{
	namespace Properties
	{

Settings *Settings::defaultInstance = (static_cast<Settings*>(System::Configuration::ApplicationSettingsBase::Synchronized(new Settings())));

		Settings *Settings::getDefault()
		{
			return defaultInstance;
		}
	}
}
