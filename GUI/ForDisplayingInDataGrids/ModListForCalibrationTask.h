#pragma once

#include <string>

namespace MetaMorpheusGUI
{
	class ModListForCalibrationTask
	{
	private:
		bool privateFixed = false;
		bool privateVariable = false;
		bool privateLocalize = false;
		std::wstring privateFileName;

//		#region Public Constructors

	public:
		ModListForCalibrationTask(const std::wstring &filePath);

//		#endregion Public Constructors

//		#region Public Properties

		bool getFixed() const;
		void setFixed(bool value);
		bool getVariable() const;
		void setVariable(bool value);
		bool getLocalize() const;
		void setLocalize(bool value);

		std::wstring getFileName() const;

//		#endregion Public Properties
	};
}
