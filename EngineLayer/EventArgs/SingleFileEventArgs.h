#pragma once

#include "MyRecursiveEventArgs.h"
#include <string>
#include <vector>


namespace EngineLayer
{
	class SingleFileEventArgs : public MyRecursiveEventArgs
	{
	private:
		std::wstring privateWrittenFile;

	public:
		SingleFileEventArgs(const std::wstring &writtenFile, std::vector<std::wstring> &nestedIds);

		std::wstring getWrittenFile() const;
		void setWrittenFile(const std::wstring &value);
	};
}
