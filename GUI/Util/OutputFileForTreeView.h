#pragma once

#include "ForTreeView.h"
#include <string>

namespace MetaMorpheusGUI
{
	class OutputFileForTreeView : public ForTreeView
	{
	private:
		std::wstring privateFullPath;

	public:
		OutputFileForTreeView(const std::wstring &fullPath, const std::wstring &displayName);

		std::wstring getFullPath() const;
	};
}
