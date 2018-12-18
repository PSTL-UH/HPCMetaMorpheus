#include "ModForTreeView.h"
#include "ModTypeForTreeView.h"


namespace MetaMorpheusGUI
{

	ModForTreeView::ModForTreeView(const std::wstring &toolTip, bool use, const std::wstring &modName, bool bad, ModTypeForTreeView *parent)
	{
		ToolTipStuff = toolTip;
		Parent = parent;
		setUse(use);
		ModName = modName;

		DisplayName = modName;

		if (StringHelper::toLower(toolTip).find(L"terminal") != std::wstring::npos)
		{
			auto split = StringHelper::split(toolTip, L'\n');

			std::wstring location = split.First([&] (std::any p)
			{
				p->StartsWith(L"PP   ");
			});
			location = StringHelper::trim(location.substr(5, location.length() - 5));

//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//			switch (location)
//ORIGINAL LINE: case "N-terminal.":
			if (location == L"N-terminal.")
			{
					DisplayName += L" (Prot N-Term)";
			}
//ORIGINAL LINE: case "C-terminal.":
			else if (location == L"C-terminal.")
			{
					DisplayName += L" (Prot C-Term)";
			}
//ORIGINAL LINE: case "Peptide N-terminal.":
			else if (location == L"Peptide N-terminal.")
			{
					DisplayName += L" (Pep N-Term)";
			}
//ORIGINAL LINE: case "Peptide C-terminal.":
			else if (location == L"Peptide C-terminal.")
			{
					DisplayName += L" (Pep C-Term)";
			}
		}

		if (bad)
		{
			Background = new SolidColorBrush(Colors::Red);
		}
		else
		{
			Background = new SolidColorBrush(Colors::Transparent);
		}
	}

	ModTypeForTreeView *ModForTreeView::getParent() const
	{
		return privateParent;
	}

	std::wstring ModForTreeView::getToolTipStuff() const
	{
		return privateToolTipStuff;
	}

	bool ModForTreeView::getUse() const
	{
		return _isChecked;
	}

	void ModForTreeView::setUse(bool value)
	{
		SetUseStatus(value);
	}

	std::wstring ModForTreeView::getModName() const
	{
		return privateModName;
	}

	std::wstring ModForTreeView::getDisplayName() const
	{
		return privateDisplayName;
	}

	Brush *ModForTreeView::getBackground() const
	{
		return privateBackground;
	}

	void ModForTreeView::SetUseStatus(bool value)
	{
		if (value == getUse())
		{
			return;
		}

		_isChecked = value;
		getParent()->VerifyCheckState();

		RaisePropertyChanged(L"Use");
	}

	void ModForTreeView::RaisePropertyChanged(const std::wstring &name)
	{
		PropertyChangedEventArgs tempVar(name);
		PropertyChanged +== nullptr ? nullptr : PropertyChanged::Invoke(this, &tempVar);
	}
}
