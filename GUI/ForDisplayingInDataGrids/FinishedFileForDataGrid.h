#pragma once

#include <string>

namespace MetaMorpheusGUI
{
	class FinishedFileForDataGrid
	{
	private:
		std::wstring privateFilePath;

//		#region Public Constructors

	public:
		FinishedFileForDataGrid(const std::wstring &filePath);

//		#endregion Public Constructors

//		#region Public Properties

		std::wstring getFilePath() const;
		void setFilePath(const std::wstring &value);

//		#endregion Public Properties
	};
}
