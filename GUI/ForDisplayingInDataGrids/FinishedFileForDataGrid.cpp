#include "FinishedFileForDataGrid.h"

namespace MetaMorpheusGUI
{

	FinishedFileForDataGrid::FinishedFileForDataGrid(const std::wstring &filePath)
	{
		setFilePath(filePath);
	}

	std::wstring FinishedFileForDataGrid::getFilePath() const
	{
		return privateFilePath;
	}

	void FinishedFileForDataGrid::setFilePath(const std::wstring &value)
	{
		privateFilePath = value;
	}
}
