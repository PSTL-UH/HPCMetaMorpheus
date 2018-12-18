#include "BooleanInverter.h"


namespace MetaMorpheusGUI
{

	std::any BooleanInverter::Convert(std::any value, std::type_info targetType, std::any parameter, System::Globalization::CultureInfo *culture)
	{
		return !std::any_cast<bool>(value);
	}

	std::any BooleanInverter::ConvertBack(std::any value, std::type_info targetType, std::any parameter, System::Globalization::CultureInfo *culture)
	{
		return std::any_cast<bool>(value);
	}
}
