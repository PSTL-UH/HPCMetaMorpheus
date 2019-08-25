#pragma once

#include "MyRecursiveEventArgs.h"
#include <string>
#include <vector>


namespace EngineLayer
{
	class SingleFileEventArgs : public MyRecursiveEventArgs
	{
	private:
		std::string privateWrittenFile;

	public:
		SingleFileEventArgs(const std::string &writtenFile, std::vector<std::string> &nestedIds);

		std::string getWrittenFile() const;
		void setWrittenFile(const std::string &value);
	};
}
