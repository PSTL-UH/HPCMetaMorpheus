#include "MetaMorpheusGUI.CustomMsgBox.h"


namespace MetaMorpheusGUI
{

	CustomMsgBox::CustomMsgBox(const std::wstring &title, const std::wstring &caption, const std::wstring &btn1, const std::wstring &btn2, const std::wstring &btn3)
	{
		InitializeComponent();
		this->Title = title;
		this->Label->Content = caption;
		this->Button1->Content = btn1;
		this->Button2->Content = btn2;
		this->Button3->Content = btn3;

		Button1::Click += new RoutedEventHandler(Button1_Click);
		Button2::Click += new RoutedEventHandler(Button2_Click);
		Button3::Click += new RoutedEventHandler(Button3_Click);
	}

CustomMsgBox *CustomMsgBox::MsgBox;
MessageBoxResult *CustomMsgBox::result = MessageBoxResult::No;

	MessageBoxResult *CustomMsgBox::Show(const std::wstring &title, const std::wstring &caption, const std::wstring &btn1, const std::wstring &btn2, const std::wstring &btn3)
	{
		MsgBox = new CustomMsgBox(title, caption, btn1, btn2, btn3);
		MsgBox->ShowDialog();
		return result;
	}

	void CustomMsgBox::Button1_Click(std::any sender, RoutedEventArgs *e)
	{
		result = MessageBoxResult::Yes;
		MsgBox->Close();
	}

	void CustomMsgBox::Button2_Click(std::any sender, RoutedEventArgs *e)
	{
		result = MessageBoxResult::No;
		MsgBox->Close();
	}

	void CustomMsgBox::Button3_Click(std::any sender, RoutedEventArgs *e)
	{
		result = MessageBoxResult::OK;
		MsgBox->Close();
	}
}
