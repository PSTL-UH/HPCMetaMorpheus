#include "ForTreeView.h"


namespace MetaMorpheusGUI
{

	ForTreeView::ForTreeView(const std::wstring &displayName, const std::wstring &id)
	{
		DisplayName = displayName;
		Id = id;
		ObservableCollection<ForTreeView*> tempVar();
		setChildren(&tempVar);
	}

	ObservableCollection<ForTreeView*> *ForTreeView::getChildren() const
	{
		return privateChildren;
	}

	void ForTreeView::setChildren(ObservableCollection<ForTreeView*> *value)
	{
		privateChildren = value;
	}

	std::wstring ForTreeView::getStatus() const
	{
		return _Status;
	}

	void ForTreeView::setStatus(const std::wstring &value)
	{
		_Status = value;
		OnPropertyChanged();
	}

	int ForTreeView::getProgress() const
	{
		return _Progress;
	}

	void ForTreeView::setProgress(int value)
	{
		_Progress = value;
		OnPropertyChanged();
	}

	std::wstring ForTreeView::getDisplayName() const
	{
		return privateDisplayName;
	}

	std::wstring ForTreeView::getId() const
	{
		return privateId;
	}

	bool ForTreeView::getIsIndeterminate() const
	{
		return _IsIndeterminate;
	}

	void ForTreeView::setIsIndeterminate(bool value)
	{
		_IsIndeterminate = value;
		OnPropertyChanged();
	}

//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: protected void OnPropertyChanged([CallerMemberName] string propertyName = null)
	void ForTreeView::OnPropertyChanged(const std::wstring &propertyName)
	{
		PropertyChangedEventArgs tempVar(propertyName);
		PropertyChanged +== nullptr ? nullptr : PropertyChanged::Invoke(this, &tempVar);
	}
}
