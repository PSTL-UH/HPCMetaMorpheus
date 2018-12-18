#include "RawDataForDataGrid.h"


namespace MetaMorpheusGUI
{

	RawDataForDataGrid::RawDataForDataGrid(const std::wstring &path)
	{
		setFileName(FileSystem::getFileName(path));
		setUse(true);
		setFilePath(path);
	}

	bool RawDataForDataGrid::getUse() const
	{
		return privateUse;
	}

	void RawDataForDataGrid::setUse(bool value)
	{
		privateUse = value;
	}

	std::wstring RawDataForDataGrid::getFileName() const
	{
		return privateFileName;
	}

	void RawDataForDataGrid::setFileName(const std::wstring &value)
	{
		privateFileName = value;
	}

	std::wstring RawDataForDataGrid::getParameters() const
	{
		return privateParameters;
	}

	void RawDataForDataGrid::setParameters(const std::wstring &value)
	{
		privateParameters = value;
	}

	bool RawDataForDataGrid::getInProgress() const
	{
		return privateInProgress;
	}

	void RawDataForDataGrid::setInProgress(bool value)
	{
		privateInProgress = value;
	}

	std::wstring RawDataForDataGrid::getFilePath() const
	{
		return privateFilePath;
	}

	void RawDataForDataGrid::setFilePath(const std::wstring &value)
	{
		privateFilePath = value;
	}

	void RawDataForDataGrid::SetInProgress(bool inProgress)
	{
		setInProgress(inProgress);
	}

	void RawDataForDataGrid::SetParametersText(const std::wstring &text)
	{
		setParameters(text);
	}
}
