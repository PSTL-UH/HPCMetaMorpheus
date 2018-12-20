#include "DbForTask.h"

namespace TaskLayer
{

	DbForTask::DbForTask(const std::wstring &filePath, bool isContaminant)
	{
		FilePath = filePath;
		IsContaminant = isContaminant;
		FileName = FileSystem::getFileName(filePath);
	}

	std::wstring DbForTask::getFilePath() const
	{
		return privateFilePath;
	}

	bool DbForTask::getIsContaminant() const
	{
		return privateIsContaminant;
	}

	std::wstring DbForTask::getFileName() const
	{
		return privateFileName;
	}
}
