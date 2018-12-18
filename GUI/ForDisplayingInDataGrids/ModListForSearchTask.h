#pragma once

#include <string>

namespace MetaMorpheusGUI
{
	class ModListForSearchTask
	{
	private:
		bool privateFixed = false;
		bool privateVariable = false;
		bool privateLocalize = false;
		bool privateAlwaysKeep = false;
		std::wstring privateFileName;

//		#region Public Constructors

	public:
		ModListForSearchTask(const std::wstring &filePath);

//		#endregion Public Constructors

//		#region Public Properties

		bool getFixed() const;
		void setFixed(bool value);
		bool getVariable() const;
		void setVariable(bool value);
		bool getLocalize() const;
		void setLocalize(bool value);
		bool getAlwaysKeep() const;
		void setAlwaysKeep(bool value);

		std::wstring getFileName() const;

//		#endregion Public Properties
	};
}
