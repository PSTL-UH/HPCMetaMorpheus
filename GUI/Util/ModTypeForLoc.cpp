#include "ModTypeForLoc.h"


namespace MetaMorpheusGUI
{

	ModTypeForLoc::ModTypeForLoc(const std::wstring &displayName)
	{
		DisplayName = displayName;
	}

	std::optional<bool> ModTypeForLoc::getUse() const
	{
		return _isChecked;
	}

	void ModTypeForLoc::setUse(const std::optional<bool> &value)
	{
		SetUseStatus(value);
	}

	std::wstring ModTypeForLoc::getDisplayName() const
	{
		return privateDisplayName;
	}

	void ModTypeForLoc::RaisePropertyChanged(const std::wstring &name)
	{
		PropertyChangedEventArgs tempVar(name);
		PropertyChanged +== nullptr ? nullptr : PropertyChanged::Invoke(this, &tempVar);
	}

	void ModTypeForLoc::SetUseStatus(std::optional<bool> &value)
	{
		_isChecked = value;

		RaisePropertyChanged(L"Use");
	}
}
