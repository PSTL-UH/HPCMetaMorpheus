#pragma once

#include <string>
#include "tangible_event.h"


namespace MetaMorpheusGUI
{
	class ForTreeView : public INotifyPropertyChanged
	{
	private:
		ObservableCollection<ForTreeView*> *privateChildren;
		std::wstring privateDisplayName;
		std::wstring privateId;

		std::wstring _Status;
		int _Progress = 0;
		bool _IsIndeterminate = false;

	public:
		ForTreeView(const std::wstring &displayName, const std::wstring &id);

		TangibleEvent<PropertyChangedEventHandler> *PropertyChanged = new TangibleEvent<PropertyChangedEventHandler>();

		ObservableCollection<ForTreeView*> *getChildren() const;
		void setChildren(ObservableCollection<ForTreeView*> *value);

		std::wstring getStatus() const;
		void setStatus(const std::wstring &value);

		int getProgress() const;
		void setProgress(int value);

		std::wstring getDisplayName() const;
		std::wstring getId() const;

		bool getIsIndeterminate() const;
		void setIsIndeterminate(bool value);

	protected:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: protected void OnPropertyChanged([CallerMemberName] string propertyName = null)
		void OnPropertyChanged(const std::wstring &propertyName = L"");
	};
}
