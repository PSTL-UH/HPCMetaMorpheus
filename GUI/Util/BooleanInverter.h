#pragma once

#include <any>
#include <typeinfo>


namespace MetaMorpheusGUI
{
	class BooleanInverter : public IValueConverter
	{
	public:
		std::any Convert(std::any value, std::type_info targetType, std::any parameter, System::Globalization::CultureInfo *culture);

		std::any ConvertBack(std::any value, std::type_info targetType, std::any parameter, System::Globalization::CultureInfo *culture);
	};
}
