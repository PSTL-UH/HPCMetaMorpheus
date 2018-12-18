#pragma once

#include <string>
#include <vector>
#include <limits>
#include <any>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class CalibrationTask; }
namespace MetaMorpheusGUI { class ModTypeForTreeView; }
namespace MetaMorpheusGUI { class ModTypeForLoc; }

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace TaskLayer;
using namespace Proteomics::ProteolyticDigestion;
using namespace MassSpectrometry;

namespace MetaMorpheusGUI
{
	/// <summary>
	/// Interaction logic for CalibrateTaskWindow.xaml
	/// </summary>
	class CalibrateTaskWindow : public Window
	{
	private:
		CalibrationTask *privateTheTask;

		ObservableCollection<ModTypeForTreeView*> *const FixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView*>();
		ObservableCollection<ModTypeForTreeView*> *const VariableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView*>();
		ObservableCollection<ModTypeForLoc*> *const LocalizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc*>();

	public:
		virtual ~CalibrateTaskWindow()
		{
			delete FixedModTypeForTreeViewObservableCollection;
			delete VariableModTypeForTreeViewObservableCollection;
			delete LocalizeModTypeForTreeViewObservableCollection;
		}

		CalibrateTaskWindow();

		CalibrateTaskWindow(CalibrationTask *myCalibrateTask);

		CalibrationTask *getTheTask() const;
		void setTheTask(CalibrationTask *value);

	private:
		void UpdateFieldsFromTask(CalibrationTask *task);

		void PopulateChoices();

		void CancelButton_Click(std::any sender, RoutedEventArgs *e);

		void SaveButton_Click(std::any sender, RoutedEventArgs *e);

		void CheckIfNumber(std::any sender, TextCompositionEventArgs *e);

		void KeyPressed(std::any sender, KeyEventArgs *e);

		void TextChanged_Fixed(std::any sender, TextChangedEventArgs *args);

		void TextChanged_Var(std::any sender, TextChangedEventArgs *args);

		void TextChangeTimerHandler(std::any sender, EventArgs *e);
	};
}
