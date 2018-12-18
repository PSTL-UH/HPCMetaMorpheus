#pragma once

#include <string>
#include "tangible_event.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace MassSpectrometry;

namespace MetaMorpheusGUI
{
	/// <summary>
	/// Interaction logic for SearchTaskWindow.xaml
	/// </summary>

	class DataContextForSearchTaskWindow : public INotifyPropertyChanged
	{
	private:
		std::wstring _ExpanderTitle;
		std::wstring _SearchModeExpanderTitle;
		std::wstring _ModExpanderTitle;
		std::wstring _AnalysisExpanderTitle;

	public:
		TangibleEvent<PropertyChangedEventHandler> *PropertyChanged = new TangibleEvent<PropertyChangedEventHandler>();

		std::wstring getExpanderTitle() const;
		void setExpanderTitle(const std::wstring &value);

		std::wstring getAnalysisExpanderTitle() const;
		void setAnalysisExpanderTitle(const std::wstring &value);

		std::wstring getModExpanderTitle() const;
		void setModExpanderTitle(const std::wstring &value);

		std::wstring getSearchModeExpanderTitle() const;
		void setSearchModeExpanderTitle(const std::wstring &value);

	protected:
		void RaisePropertyChanged(const std::wstring &name);
	};
}
