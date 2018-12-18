#pragma once

#include <string>
#include "stringhelper.h"
#include "tangible_event.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace MetaMorpheusGUI { class ModTypeForTreeView; }


namespace MetaMorpheusGUI
{
	class ModForTreeView : public INotifyPropertyChanged
	{
	private:
		ModTypeForTreeView *privateParent;
		std::wstring privateToolTipStuff;
		std::wstring privateModName;
		std::wstring privateDisplayName;
		Brush *privateBackground;

		bool _isChecked = false;

	public:
		ModForTreeView(const std::wstring &toolTip, bool use, const std::wstring &modName, bool bad, ModTypeForTreeView *parent);

		TangibleEvent<PropertyChangedEventHandler> *PropertyChanged = new TangibleEvent<PropertyChangedEventHandler>();

		ModTypeForTreeView *getParent() const;
		std::wstring getToolTipStuff() const;

		bool getUse() const;
		void setUse(bool value);

		std::wstring getModName() const;
		std::wstring getDisplayName() const;
		Brush *getBackground() const;

		void SetUseStatus(bool value);

	protected:
		void RaisePropertyChanged(const std::wstring &name);
	};
}
