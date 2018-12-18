#pragma once

#include <string>
#include <any>


namespace MetaMorpheusGUI
{
	/// <summary>
	/// Interaction logic for CustomMsgBox.xaml
	/// </summary>
	class CustomMsgBox : public Window
	{
	public:
		CustomMsgBox(const std::wstring &title, const std::wstring &caption, const std::wstring &btn1, const std::wstring &btn2, const std::wstring &btn3);

	private:
		static CustomMsgBox *MsgBox;
		static MessageBoxResult *result;

	public:
		static MessageBoxResult *Show(const std::wstring &title, const std::wstring &caption, const std::wstring &btn1, const std::wstring &btn2, const std::wstring &btn3);

	private:
		void Button1_Click(std::any sender, RoutedEventArgs *e);

		void Button2_Click(std::any sender, RoutedEventArgs *e);

		void Button3_Click(std::any sender, RoutedEventArgs *e);
	};
}
