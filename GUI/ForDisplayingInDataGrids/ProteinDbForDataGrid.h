#pragma once

#include <string>
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class DbForTask; }

using namespace TaskLayer;

namespace MetaMorpheusGUI
{
	class ProteinDbForDataGrid
	{
	private:
		bool privateUse = false;
		bool privateContaminant = false;
		std::wstring privateFileName;
		std::wstring privateFilePath;
		bool privateInProgress = false;

//		#region Public Constructors

	public:
		ProteinDbForDataGrid(const std::wstring &FilePath);

		ProteinDbForDataGrid(DbForTask *uu);

//		#endregion Public Constructors

//		#region Public Properties

		bool getUse() const;
		void setUse(bool value);
		bool getContaminant() const;
		void setContaminant(bool value);
		std::wstring getFileName() const;
		void setFileName(const std::wstring &value);
		std::wstring getFilePath() const;
		void setFilePath(const std::wstring &value);
		bool getInProgress() const;
		void setInProgress(bool value);

//		#endregion Public Properties

//		#region Public Methods

		/// <summary>
		/// Method to mark as in progress. Need the property setter to be private so user could not check off in GUI
		/// </summary>
		/// <param name="inProgress"></param>
		void SetInProgress(bool inProgress);

//		#endregion Public Methods
	};
}
