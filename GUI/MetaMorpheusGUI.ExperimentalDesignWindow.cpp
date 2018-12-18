#include "MetaMorpheusGUI.ExperimentalDesignWindow.h"
#include "ForDisplayingInDataGrids/ExperimentalDesignForDataGrid.h"
#include "ForDisplayingInDataGrids/RawDataForDataGrid.h"
#include "../EngineLayer/GlobalVariables.h"

using namespace EngineLayer;

namespace MetaMorpheusGUI
{

	ExperimentalDesignWindow::ExperimentalDesignWindow(ObservableCollection<RawDataForDataGrid*> *spectraFilesObservableCollection)
	{
		InitializeComponent();

		for (auto item : spectraFilesObservableCollection->Where([&] (std::any p)
		{
			p::Use;
		})->Select([&] (std::any p)
		{
			p::FilePath;
		}))
		{
		}
	}

	private *ExperimentalDesignWindow::if_Renamed(spectraFilesObservableCollection::Any())
	{
		outputPath = Directory::GetParent(spectraFilesObservableCollection::Where([&] (std::any p)
		{
			p::Use;
		}).First().FilePath)->FullName;
		outputPath = FileSystem::combine(outputPath, GlobalVariables::getExperimentalDesignFileName());
	}
}
