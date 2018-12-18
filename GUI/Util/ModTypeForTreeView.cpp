#include "ModTypeForTreeView.h"
#include "ModForTreeView.h"


namespace MetaMorpheusGUI
{

	ModTypeForTreeView::ModTypeForTreeView(const std::wstring &displayName, bool bad)
	{
		Children = new ObservableCollection<ModForTreeView*>();
		setExpanded(false);
		DisplayName = displayName;
		if (bad)
		{
			Background = new SolidColorBrush(Colors::Red);
		}
		else
		{
			Background = new SolidColorBrush(Colors::Transparent);
		}
	}

	std::optional<bool> ModTypeForTreeView::getUse() const
	{
		return _isChecked;
	}

	void ModTypeForTreeView::setUse(const std::optional<bool> &value)
	{
		SetUseStatus(value);
	}

	bool ModTypeForTreeView::getExpanded() const
	{
		return privateExpanded;
	}

	void ModTypeForTreeView::setExpanded(bool value)
	{
		privateExpanded = value;
	}

	std::wstring ModTypeForTreeView::getDisplayName() const
	{
		return privateDisplayName;
	}

	ObservableCollection<ModForTreeView*> *ModTypeForTreeView::getChildren() const
	{
		return privateChildren;
	}

	Brush *ModTypeForTreeView::getBackground() const
	{
		return privateBackground;
	}

	void ModTypeForTreeView::VerifyCheckState()
	{
		std::optional<bool> state = std::nullopt;
		for (int i = 0; i < getChildren()->Count; ++i)
		{
			bool current = getChildren()[i]->getUse();
			if (i == 0)
			{
				state = std::make_optional(current);
			}
			else if (state != current)
			{
				state = std::nullopt;
				break;
			}
		}
		SetUseStatus(state);
	}

	void ModTypeForTreeView::RaisePropertyChanged(const std::wstring &name)
	{
		PropertyChangedEventArgs tempVar(name);
		PropertyChanged +== nullptr ? nullptr : PropertyChanged::Invoke(this, &tempVar);
	}

	void ModTypeForTreeView::SetUseStatus(std::optional<bool> &value)
	{
		if (value == getUse())
		{
			return;
		}

		_isChecked = value;

		if (value)
		{
			for (auto child : getChildren())
			{
				child->SetUseStatus(value.value());
			}
		}

		RaisePropertyChanged(L"Use");
	}
}
