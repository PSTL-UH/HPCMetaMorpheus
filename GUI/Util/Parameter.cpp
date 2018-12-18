#include "Parameter.h"

using namespace EngineLayer;
using namespace Proteomics::ProteolyticDigestion;

namespace MetaMorpheusGUI
{

	Parameter::Parameter()
	{
	}

	Parameter::Parameter(const std::wstring &name, const std::wstring &valueType)
	{
		setParamName(name);
		setValueType(valueType);

		ObservableCollection<Protease*> tempVar();
		setProtList(&tempVar);
		ObservableCollection<std::wstring> tempVar2();
		setInitList(&tempVar2);
		ObservableCollection<std::wstring> tempVar3();
		setProductMassToleranceList(&tempVar3);

		for (Protease *protease : ProteaseDictionary::Dictionary->Values)
		{
			getProtList()->Add(protease);
		}
		for (auto initiatior_methionine_behavior : Enum::GetNames(InitiatorMethionineBehavior::typeid))
		{
			getInitList()->Add(initiatior_methionine_behavior);
		}
		getProductMassToleranceList()->Add(L"Absolute");
		getProductMassToleranceList()->Add(L"ppm");
	}

	std::wstring Parameter::getParamName() const
	{
		return privateParamName;
	}

	void Parameter::setParamName(const std::wstring &value)
	{
		privateParamName = value;
	}

	std::wstring Parameter::getValueType() const
	{
		return privateValueType;
	}

	void Parameter::setValueType(const std::wstring &value)
	{
		privateValueType = value;
	}

	std::any Parameter::getValue() const
	{
		return _value;
	}

	void Parameter::setValue(const std::any &value)
	{
		_value = value;
		OnPropertyChanged(L"Value");
	}

	bool Parameter::getDifferent() const
	{
		return privateDifferent;
	}

	void Parameter::setDifferent(bool value)
	{
		privateDifferent = value;
	}

	bool Parameter::getHasChanged() const
	{
		return Status;
	}

	void Parameter::setHasChanged(bool value)
	{
		Status = value;
		if (value == Status)
		{
			return;
		}
		Status = value;
	}

	ObservableCollection<Protease*> *Parameter::getProtList() const
	{
		return privateProtList;
	}

	void Parameter::setProtList(ObservableCollection<Protease*> *value)
	{
		privateProtList = value;
	}

	ObservableCollection<std::wstring> *Parameter::getInitList() const
	{
		return privateInitList;
	}

	void Parameter::setInitList(ObservableCollection<std::wstring> *value)
	{
		privateInitList = value;
	}

	ObservableCollection<std::wstring> *Parameter::getProductMassToleranceList() const
	{
		return privateProductMassToleranceList;
	}

	void Parameter::setProductMassToleranceList(ObservableCollection<std::wstring> *value)
	{
		privateProductMassToleranceList = value;
	}

	void Parameter::OnPropertyChanged(const std::wstring &propertyName)
	{
		PropertyChangedEventHandler handler = PropertyChanged;
		if (handler != nullptr)
		{
			PropertyChangedEventArgs tempVar(propertyName);
			handler(this, &tempVar);
		}

		this->setHasChanged(true);
	}
}
