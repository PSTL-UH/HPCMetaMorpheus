#pragma once

#include <string>
#include <vector>
#include <limits>
#include <any>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class GptmdTask; }
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
	/// Interaction logic for GptmdTaskWindow.xaml
	/// </summary>
	class GptmdTaskWindow : public Window
	{
	private:
		GptmdTask *privateTheTask;

		ObservableCollection<ModTypeForTreeView*> *const fixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView*>();
		ObservableCollection<ModTypeForTreeView*> *const variableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView*>();
		ObservableCollection<ModTypeForLoc*> *const localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc*>();
		ObservableCollection<ModTypeForTreeView*> *const gptmdModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView*>();

	public:
		virtual ~GptmdTaskWindow()
		{
			delete fixedModTypeForTreeViewObservableCollection;
			delete variableModTypeForTreeViewObservableCollection;
			delete localizeModTypeForTreeViewObservableCollection;
			delete gptmdModTypeForTreeViewObservableCollection;
		}

		GptmdTaskWindow();

		GptmdTaskWindow(GptmdTask *myGPTMDtask);

		GptmdTask *getTheTask() const;
		void setTheTask(GptmdTask *value);

	private:
		void Row_DoubleClick(std::any sender, MouseButtonEventArgs *e);

		void UpdateFieldsFromTask(GptmdTask *task);

		void PopulateChoices();

		void CancelButton_Click(std::any sender, RoutedEventArgs *e);

		void SaveButton_Click(std::any sender, RoutedEventArgs *e);

		void CheckIfNumber(std::any sender, TextCompositionEventArgs *e);

		void KeyPressed(std::any sender, KeyEventArgs *e);

		void TextChanged_Fixed(std::any sender, TextChangedEventArgs *args);

		void TextChanged_Var(std::any sender, TextChangedEventArgs *args);

		void TextChanged_GPTMD(std::any sender, TextChangedEventArgs *args);

		void TextChangeTimerHandler(std::any sender, EventArgs *e);
	};
}
