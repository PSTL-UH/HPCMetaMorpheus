#include "SearchModeForDataGrid.h"
#include "../../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"

using namespace EngineLayer;

namespace MetaMorpheusGUI
{

	SearchModeForDataGrid::SearchModeForDataGrid(MassDiffAcceptor *searchMode) : searchMode(searchMode)
	{
	}

	bool SearchModeForDataGrid::getUse() const
	{
		return privateUse;
	}

	void SearchModeForDataGrid::setUse(bool value)
	{
		privateUse = value;
	}

	std::wstring SearchModeForDataGrid::getName() const
	{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		return searchMode->ToString();
	}
}
