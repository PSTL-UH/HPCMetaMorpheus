#include "SearchTaskWindow.xaml.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace MassSpectrometry;

namespace MetaMorpheusGUI
{

	std::wstring DataContextForSearchTaskWindow::getExpanderTitle() const
	{
		return _ExpanderTitle;
	}

	void DataContextForSearchTaskWindow::setExpanderTitle(const std::wstring &value)
	{
		_ExpanderTitle = value;
		RaisePropertyChanged(L"ExpanderTitle");
	}

	std::wstring DataContextForSearchTaskWindow::getAnalysisExpanderTitle() const
	{
		return _AnalysisExpanderTitle;
	}

	void DataContextForSearchTaskWindow::setAnalysisExpanderTitle(const std::wstring &value)
	{
		_AnalysisExpanderTitle = value;
		RaisePropertyChanged(L"AnalysisExpanderTitle");
	}

	std::wstring DataContextForSearchTaskWindow::getModExpanderTitle() const
	{
		return _ModExpanderTitle;
	}

	void DataContextForSearchTaskWindow::setModExpanderTitle(const std::wstring &value)
	{
		_ModExpanderTitle = value;
		RaisePropertyChanged(L"ModExpanderTitle");
	}

	std::wstring DataContextForSearchTaskWindow::getSearchModeExpanderTitle() const
	{
		return _SearchModeExpanderTitle;
	}

	void DataContextForSearchTaskWindow::setSearchModeExpanderTitle(const std::wstring &value)
	{
		_SearchModeExpanderTitle = value;
		RaisePropertyChanged(L"SearchModeExpanderTitle");
	}

	void DataContextForSearchTaskWindow::RaisePropertyChanged(const std::wstring &name)
	{
		PropertyChangedEventArgs tempVar(name);
		PropertyChanged +== nullptr ? nullptr : PropertyChanged::Invoke(this, &tempVar);
	}
}
