#pragma once

#include <string>
#include <any>
#include "tangible_event.h"

using namespace EngineLayer;
using namespace Proteomics::ProteolyticDigestion;

namespace MetaMorpheusGUI
{
	class Parameter : public INotifyPropertyChanged
	{
	private:
		std::wstring privateParamName;
		std::wstring privateValueType;
		bool privateDifferent = false;
		ObservableCollection<Protease*> *privateProtList;
		ObservableCollection<std::wstring> *privateInitList;
		ObservableCollection<std::wstring> *privateProductMassToleranceList;

		std::any _value;
		bool Status = false;

	public:
		Parameter();

		Parameter(const std::wstring &name, const std::wstring &valueType);

		TangibleEvent<PropertyChangedEventHandler> *PropertyChanged = new TangibleEvent<PropertyChangedEventHandler>();

		std::wstring getParamName() const;
		void setParamName(const std::wstring &value);

		std::wstring getValueType() const;
		void setValueType(const std::wstring &value);

		std::any getValue() const;
		void setValue(const std::any &value);

		bool getDifferent() const;
		void setDifferent(bool value);

		bool getHasChanged() const;
		void setHasChanged(bool value);

		ObservableCollection<Protease*> *getProtList() const;
		void setProtList(ObservableCollection<Protease*> *value);

		ObservableCollection<std::wstring> *getInitList() const;
		void setInitList(ObservableCollection<std::wstring> *value);

		ObservableCollection<std::wstring> *getProductMassToleranceList() const;
		void setProductMassToleranceList(ObservableCollection<std::wstring> *value);

	protected:
		virtual void OnPropertyChanged(const std::wstring &propertyName);
	};
}
