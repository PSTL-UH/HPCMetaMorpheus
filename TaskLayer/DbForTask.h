#pragma once

#include <string>
#include "tangible_filesystem.h"

namespace TaskLayer
{
	class DbForTask
	{
	private:
		std::string privateFilePath;
		bool privateIsContaminant = false;
		std::string privateFileName;

	public:
		DbForTask(const std::string &filePath, bool isContaminant);

		std::string getFilePath() const;
		bool getIsContaminant() const;
		std::string getFileName() const;
	};
}
