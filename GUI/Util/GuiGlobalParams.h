#pragma once

namespace MetaMorpheusGUI
{
	class GuiGlobalParams
	{
	private:
		bool privateAskAboutUpdating = true;
		bool privateDisableCloseWindow = false;

	public:
		bool getAskAboutUpdating() const;
		void setAskAboutUpdating(bool value);
		bool getDisableCloseWindow() const;
		void setDisableCloseWindow(bool value);
	};
}
