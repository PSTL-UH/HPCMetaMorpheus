#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <cctype>
#include <stdexcept>
#include <any>
#include <optional>
#include "stringhelper.h"
#include "tangible_filesystem.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace Proteomics;

namespace MetaMorpheusGUI
{
	/// <summary>
	/// Interaction logic for AddCustomModWindow.xaml
	/// </summary>
	class CustomModButtonWindow : public Window
	{
	public:
		static std::unordered_map<std::wstring, std::wstring> locationRestrictions;

		CustomModButtonWindow();

		void SaveCustomMod_Click(std::any sender, RoutedEventArgs *e);

	private:
		void CancelCustomMod_Click(std::any sender, RoutedEventArgs *e);

		bool ErrorsDetected(const std::wstring &myModName, const std::wstring &motif, const std::wstring &mass, const std::wstring &chemFormula, const std::wstring &neutralLoss, const std::wstring &modType, const std::wstring &diagnosticIon);
	};
}
