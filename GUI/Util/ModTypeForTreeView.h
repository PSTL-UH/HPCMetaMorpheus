#pragma once

#include <string>
#include <optional>
#include "tangible_event.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace MetaMorpheusGUI { class ModForTreeView; }


namespace MetaMorpheusGUI
{
	class ModTypeForTreeView : public INotifyPropertyChanged
	{
	private:
		bool privateExpanded = false;
		std::wstring privateDisplayName;
		ObservableCollection<ModForTreeView*> *privateChildren;
		Brush *privateBackground;

		std::optional<bool> _isChecked;

	public:
		ModTypeForTreeView(const std::wstring &displayName, bool bad);

		TangibleEvent<PropertyChangedEventHandler> *PropertyChanged = new TangibleEvent<PropertyChangedEventHandler>();

		std::optional<bool> getUse() const;
		void setUse(const std::optional<bool> &value);

		bool getExpanded() const;
		void setExpanded(bool value);

		std::wstring getDisplayName() const;

		ObservableCollection<ModForTreeView*> *getChildren() const;

		Brush *getBackground() const;

		void VerifyCheckState();

	protected:
		void RaisePropertyChanged(const std::wstring &name);

	private:
		void SetUseStatus(std::optional<bool> &value);
	};
}
