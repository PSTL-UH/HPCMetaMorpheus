#pragma once

#include "ForTreeView.h"
#include <string>

namespace MetaMorpheusGUI
{
	class CollectionForTreeView : public ForTreeView
	{
	public:
		CollectionForTreeView(const std::wstring &displayName, const std::wstring &id);
	};
}
