#include "DbForTask.h"

namespace TaskLayer
{

	DbForTask::DbForTask(const std::string &filePath, bool isContaminant)
	{
            privateFilePath = filePath;
            privateIsContaminant = isContaminant;
            privateFileName = FileSystem::getFileName(filePath);
	}

	std::string DbForTask::getFilePath() const
	{
		return privateFilePath;
	}

	bool DbForTask::getIsContaminant() const
	{
		return privateIsContaminant;
	}

	std::string DbForTask::getFileName() const
	{
		return privateFileName;
	}
}
