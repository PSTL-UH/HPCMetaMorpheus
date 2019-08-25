#include "SingleFileEventArgs.h"


namespace EngineLayer
{

	SingleFileEventArgs::SingleFileEventArgs(const std::string &writtenFile, std::vector<std::string> &nestedIds) : MyRecursiveEventArgs(nestedIds)
	{
		setWrittenFile(writtenFile);
	}

	std::string SingleFileEventArgs::getWrittenFile() const
	{
		return privateWrittenFile;
	}

	void SingleFileEventArgs::setWrittenFile(const std::string &value)
	{
		privateWrittenFile = value;
	}
}
