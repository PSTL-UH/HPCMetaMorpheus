#include "SingleFileEventArgs.h"


namespace EngineLayer
{

	SingleFileEventArgs::SingleFileEventArgs(const std::wstring &writtenFile, std::vector<std::wstring> &nestedIds) : MyRecursiveEventArgs(nestedIds)
	{
		setWrittenFile(writtenFile);
	}

	std::wstring SingleFileEventArgs::getWrittenFile() const
	{
		return privateWrittenFile;
	}

	void SingleFileEventArgs::setWrittenFile(const std::wstring &value)
	{
		privateWrittenFile = value;
	}
}
