#pragma once

#include <string>
#include <optional>
#include "tangible_event.h"


namespace MetaMorpheusGUI
{
	class ModTypeForLoc : public INotifyPropertyChanged
	{
	private:
		std::wstring privateDisplayName;

		std::optional<bool> _isChecked;

	public:
		ModTypeForLoc(const std::wstring &displayName);

		TangibleEvent<PropertyChangedEventHandler> *PropertyChanged = new TangibleEvent<PropertyChangedEventHandler>();

		std::optional<bool> getUse() const;
		void setUse(const std::optional<bool> &value);

		std::wstring getDisplayName() const;

	protected:
		void RaisePropertyChanged(const std::wstring &name);

	private:
		void SetUseStatus(std::optional<bool> &value);
	};
}
