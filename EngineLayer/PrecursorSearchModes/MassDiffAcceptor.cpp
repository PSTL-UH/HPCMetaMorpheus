#include "MassDiffAcceptor.h"
#include "AllowedIntervalWithNotch.h"


namespace EngineLayer
{

	MassDiffAcceptor::MassDiffAcceptor(const std::wstring &fileNameAddition)
	{
		FileNameAddition = fileNameAddition;
		setNumNotches(1);
	}

	int MassDiffAcceptor::getNumNotches() const
	{
		return privateNumNotches;
	}

	void MassDiffAcceptor::setNumNotches(int value)
	{
		privateNumNotches = value;
	}

	std::wstring MassDiffAcceptor::getFileNameAddition() const
	{
		return privateFileNameAddition;
	}
}
