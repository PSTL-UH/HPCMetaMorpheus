#pragma once

#include <string>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class MetaMorpheusTask; }

using namespace TaskLayer;

namespace MetaMorpheusGUI
{
	class PreRunTask
	{
	private:
		std::wstring privateDisplayName;

//		#region Public Fields

	public:
		MetaMorpheusTask *const metaMorpheusTask;

//		#endregion Public Fields

//		#region Public Constructors

		virtual ~PreRunTask()
		{
			delete metaMorpheusTask;
		}

		PreRunTask(MetaMorpheusTask *theTask);

//		#endregion Public Constructors

//		#region Public Properties

		std::wstring getDisplayName() const;
		void setDisplayName(const std::wstring &value);

//		#endregion Public Properties

//		#region Public Methods

		PreRunTask *Clone();

//		#endregion Public Methods
	};
}
