#include "OutputFileForTreeView.h"

namespace MetaMorpheusGUI
{

	OutputFileForTreeView::OutputFileForTreeView(const std::wstring &fullPath, const std::wstring &displayName) : ForTreeView(displayName, fullPath)
	{
		FullPath = fullPath;
	}

	std::wstring OutputFileForTreeView::getFullPath() const
	{
		return privateFullPath;
	}
}
