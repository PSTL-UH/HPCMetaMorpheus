#pragma once

#include <string>

//Object that is placed in the DataGrid for ModSelection when writing a pruned DB
namespace MetaMorpheusGUI
{
	class ModTypeForGrid
	{
	private:
		std::wstring privateModName;
		bool privateItem2 = false;
		bool privateItem3 = false;
		bool privateItem4 = false;
		bool privateItem5 = false;

	public:
		ModTypeForGrid(const std::wstring &modName);

		//types
		std::wstring getModName() const;
		void setModName(const std::wstring &value);

		bool getItem2() const;
		void setItem2(bool value);

		bool getItem3() const;
		void setItem3(bool value);

		bool getItem4() const;
		void setItem4(bool value);

		bool getItem5() const;
		void setItem5(bool value);
	};
}
