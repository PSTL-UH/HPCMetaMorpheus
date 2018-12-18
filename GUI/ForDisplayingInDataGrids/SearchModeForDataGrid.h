#pragma once

#include <string>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MassDiffAcceptor; }

using namespace EngineLayer;

namespace MetaMorpheusGUI
{
	class SearchModeForDataGrid
	{
	private:
		bool privateUse = false;

//		#region Public Fields

	public:
		MassDiffAcceptor *const searchMode;

//		#endregion Public Fields

//		#region Public Constructors

		virtual ~SearchModeForDataGrid()
		{
			delete searchMode;
		}

		SearchModeForDataGrid(MassDiffAcceptor *searchMode);

//		#endregion Public Constructors

//		#region Public Properties

		bool getUse() const;
		void setUse(bool value);

		std::wstring getName() const;

//		#endregion Public Properties
	};
}
