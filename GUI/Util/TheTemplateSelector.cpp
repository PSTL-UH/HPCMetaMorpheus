#include "TheTemplateSelector.h"
#include "Parameter.h"


namespace MetaMorpheusGUI
{

	DataTemplate *TheTemplateSelector::getComboBoxProtease() const
	{
		return privateComboBoxProtease;
	}

	void TheTemplateSelector::setComboBoxProtease(DataTemplate *value)
	{
		privateComboBoxProtease = value;
	}

	DataTemplate *TheTemplateSelector::getComboBoxInit() const
	{
		return privateComboBoxInit;
	}

	void TheTemplateSelector::setComboBoxInit(DataTemplate *value)
	{
		privateComboBoxInit = value;
	}

	DataTemplate *TheTemplateSelector::getComboBoxTolerance() const
	{
		return privateComboBoxTolerance;
	}

	void TheTemplateSelector::setComboBoxTolerance(DataTemplate *value)
	{
		privateComboBoxTolerance = value;
	}

	DataTemplate *TheTemplateSelector::getBool() const
	{
		return privateBool;
	}

	void TheTemplateSelector::setBool(DataTemplate *value)
	{
		privateBool = value;
	}

	DataTemplate *TheTemplateSelector::getTextBox() const
	{
		return privateTextBox;
	}

	void TheTemplateSelector::setTextBox(DataTemplate *value)
	{
		privateTextBox = value;
	}

	System::Windows::DataTemplate *TheTemplateSelector::SelectTemplate(std::any item, System::Windows::DependencyObject *container)
	{
		if (dynamic_cast<Parameter*>(item) != nullptr)
		{
			Parameter *settings = dynamic_cast<Parameter*>(item);
			if (settings->getValueType() == L"ComboBoxProtease")
			{
				return getComboBoxProtease();
			}
			else if (settings->getValueType() == L"Bool")
			{
				return getBool();
			}
			else if (settings->getValueType() == L"ComboBoxInit")
			{
				return getComboBoxInit();
			}
			else if (settings->getValueType() == L"ProductMassToleranceList")
			{
				return getComboBoxTolerance();
			}
			else if (settings->getValueType() == L"TextBox")
			{
				return getTextBox();
			}
		}
		return nullptr;
	}
}
