#pragma once

#include <string>
#include "tangible_filesystem.h"

namespace TaskLayer
{
	class DbForTask
	{
	private:
		std::wstring privateFilePath;
		bool privateIsContaminant = false;
		std::wstring privateFileName;

	public:
		DbForTask(const std::wstring &filePath, bool isContaminant);

		std::wstring getFilePath() const;
		bool getIsContaminant() const;
		std::wstring getFileName() const;
	};
}
