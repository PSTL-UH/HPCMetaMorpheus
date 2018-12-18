#include "ProteinDbForDataGrid.h"
#include "../../TaskLayer/DbForTask.h"

using namespace TaskLayer;

namespace MetaMorpheusGUI
{

	ProteinDbForDataGrid::ProteinDbForDataGrid(const std::wstring &FilePath)
	{
		setUse(true);
		this->setFilePath(FilePath);
//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
		if (FilePath.ToUpperInvariant().find((std::wstring(L"contaminant")).ToUpperInvariant()) != std::wstring::npos || FilePath.ToUpperInvariant().find(L"CRAP") != std::wstring::npos)
		{
			setContaminant(true);
		}
		setFileName(FileSystem::getFileName(FilePath));
	}

	ProteinDbForDataGrid::ProteinDbForDataGrid(DbForTask *uu)
	{
		setUse(true);
		setContaminant(uu->getIsContaminant());
		setFilePath(uu->getFilePath());
		setFileName(uu->getFileName());
	}

	bool ProteinDbForDataGrid::getUse() const
	{
		return privateUse;
	}

	void ProteinDbForDataGrid::setUse(bool value)
	{
		privateUse = value;
	}

	bool ProteinDbForDataGrid::getContaminant() const
	{
		return privateContaminant;
	}

	void ProteinDbForDataGrid::setContaminant(bool value)
	{
		privateContaminant = value;
	}

	std::wstring ProteinDbForDataGrid::getFileName() const
	{
		return privateFileName;
	}

	void ProteinDbForDataGrid::setFileName(const std::wstring &value)
	{
		privateFileName = value;
	}

	std::wstring ProteinDbForDataGrid::getFilePath() const
	{
		return privateFilePath;
	}

	void ProteinDbForDataGrid::setFilePath(const std::wstring &value)
	{
		privateFilePath = value;
	}

	bool ProteinDbForDataGrid::getInProgress() const
	{
		return privateInProgress;
	}

	void ProteinDbForDataGrid::setInProgress(bool value)
	{
		privateInProgress = value;
	}

	void ProteinDbForDataGrid::SetInProgress(bool inProgress)
	{
		setInProgress(inProgress);
	}
}
