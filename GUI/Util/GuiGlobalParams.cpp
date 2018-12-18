#include "GuiGlobalParams.h"

namespace MetaMorpheusGUI
{

	bool GuiGlobalParams::getAskAboutUpdating() const
	{
		return privateAskAboutUpdating;
	}

	void GuiGlobalParams::setAskAboutUpdating(bool value)
	{
		privateAskAboutUpdating = value;
	}

	bool GuiGlobalParams::getDisableCloseWindow() const
	{
		return privateDisableCloseWindow;
	}

	void GuiGlobalParams::setDisableCloseWindow(bool value)
	{
		privateDisableCloseWindow = value;
	}
}
