#pragma once

#include <string>
#include <vector>


namespace EngineLayer
{
	class StringListEventArgs : public EventArgs
	{
	private:
		std::vector<std::wstring> privateStringList;

	public:
		StringListEventArgs(std::vector<std::wstring> &stringList);

		std::vector<std::wstring> getStringList() const;
	};
}
