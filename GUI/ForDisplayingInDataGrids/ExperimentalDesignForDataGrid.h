#pragma once

#include <string>


namespace MetaMorpheusGUI
{
	class ExperimentalDesignForDataGrid
	{
	private:
		std::wstring privateFileName;
		std::wstring privateCondition;
		std::wstring privateBiorep;
		std::wstring privateFraction;
		std::wstring privateTechrep;

	public:
		ExperimentalDesignForDataGrid(const std::wstring &filename);

		std::wstring getFileName() const;
		void setFileName(const std::wstring &value);
		std::wstring getCondition() const;
		void setCondition(const std::wstring &value);
		std::wstring getBiorep() const;
		void setBiorep(const std::wstring &value);
		std::wstring getFraction() const;
		void setFraction(const std::wstring &value);
		std::wstring getTechrep() const;
		void setTechrep(const std::wstring &value);

		void SetQconditionText(const std::wstring &text);
		void SetQbioRepText(const std::wstring &text);
		void SetQfractionText(const std::wstring &text);
		void SetQtechRepText(const std::wstring &text);
	};
}
