#include "ProteinParsimonyResults.h"
#include "ProteinGroup.h"
#include "ProteinParsimonyEngine.h"


namespace EngineLayer
{

	ProteinParsimonyResults::ProteinParsimonyResults(ProteinParsimonyEngine *proteinAnalysisEngine) : MetaMorpheusEngineResults(proteinAnalysisEngine)
	{
	}

	std::vector<ProteinGroup*> ProteinParsimonyResults::getProteinGroups() const
	{
		return privateProteinGroups;
	}

	void ProteinParsimonyResults::setProteinGroups(const std::vector<ProteinGroup*> &value)
	{
		privateProteinGroups = value;
	}

	std::wstring ProteinParsimonyResults::ToString()
	{
		auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		sb->appendLine(MetaMorpheusEngineResults::ToString());

		delete sb;
		return sb->toString();
	}
}
