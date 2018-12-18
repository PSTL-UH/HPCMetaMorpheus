#pragma once

#include <string>
#include <limits>
#include <any>
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace MetaMorpheusGUI { class RawDataForDataGrid; }

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Nett;
using namespace TaskLayer;
using namespace Proteomics::ProteolyticDigestion;

namespace MetaMorpheusGUI
{
	/// <summary>
	/// Interaction logic for ChangeParametersWindow.xaml
	/// </summary>
	class FileSpecificParametersWindow : public Window
	{
	private:
		ObservableCollection<RawDataForDataGrid*> *privateSelectedSpectra;


		//Window that is opened if user wishes to change file specific settings (TOML) for 
		//individual or multiple spectra files. Creates a toml file where settings can be
		//viewed, loaded, and changed from it.
	public:
		FileSpecificParametersWindow(ObservableCollection<RawDataForDataGrid*> *selectedSpectraFiles);

		ObservableCollection<RawDataForDataGrid*> *getSelectedSpectra() const;
		void setSelectedSpectra(ObservableCollection<RawDataForDataGrid*> *value);

		// write the toml settings file on clicking "save"
	private:
		void Save_Click(std::any sender, RoutedEventArgs *e);

		// exits dialog; nothing is written
		void Cancel_Click(std::any sender, RoutedEventArgs *e);

		void PopulateChoices();

		void CheckIfNumber(std::any sender, TextCompositionEventArgs *e);

		void KeyPressed(std::any sender, KeyEventArgs *e);
	};
}
