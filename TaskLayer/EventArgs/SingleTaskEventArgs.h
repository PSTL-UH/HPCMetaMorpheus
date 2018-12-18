#pragma once

#include <string>


namespace TaskLayer
{
	class SingleTaskEventArgs : public EventArgs
	{
	private:
		std::wstring privateDisplayName;

	public:
		SingleTaskEventArgs(const std::wstring &displayName);

		std::wstring getDisplayName() const;
		void setDisplayName(const std::wstring &value);
	};
}
