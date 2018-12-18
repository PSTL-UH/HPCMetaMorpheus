#pragma once

#include <string>
#include "tangible_filesystem.h"


namespace MetaMorpheusGUI
{
	class RawDataForDataGrid
	{
	private:
		bool privateUse = false;
		std::wstring privateFileName;
		std::wstring privateParameters;
		bool privateInProgress = false;
		std::wstring privateFilePath;

//		#region Public Constructors

	public:
		RawDataForDataGrid(const std::wstring &path);

//		#endregion Public Constructors

//		#region Public Properties

		bool getUse() const;
		void setUse(bool value);
		std::wstring getFileName() const;
		void setFileName(const std::wstring &value);
		std::wstring getParameters() const;
		void setParameters(const std::wstring &value);
		bool getInProgress() const;
		void setInProgress(bool value);
		std::wstring getFilePath() const;
		void setFilePath(const std::wstring &value);



//		#endregion Public Properties

//		#region Public Methods

		/// <summary>
		/// Method to mark as in progress. Need the property setter to be private so user could not check off in GUI
		/// </summary>
		/// <param name="inProgress"></param>
		void SetInProgress(bool inProgress);

		void SetParametersText(const std::wstring &text);


//		#endregion Public Methods
	};
}
