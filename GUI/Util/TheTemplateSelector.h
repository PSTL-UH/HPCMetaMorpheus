#pragma once

#include <any>


namespace MetaMorpheusGUI
{
	class TheTemplateSelector : public DataTemplateSelector
	{
	private:
		DataTemplate *privateComboBoxProtease;
		DataTemplate *privateComboBoxInit;
		DataTemplate *privateComboBoxTolerance;
		DataTemplate *privateBool;
		DataTemplate *privateTextBox;

	public:
		DataTemplate *getComboBoxProtease() const;
		void setComboBoxProtease(DataTemplate *value);
		DataTemplate *getComboBoxInit() const;
		void setComboBoxInit(DataTemplate *value);
		DataTemplate *getComboBoxTolerance() const;
		void setComboBoxTolerance(DataTemplate *value);
		DataTemplate *getBool() const;
		void setBool(DataTemplate *value);
		DataTemplate *getTextBox() const;
		void setTextBox(DataTemplate *value);

		System::Windows::DataTemplate *SelectTemplate(std::any item, System::Windows::DependencyObject *container) override;
	};
}
