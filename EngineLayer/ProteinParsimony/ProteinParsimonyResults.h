#pragma once

#include "../MetaMorpheusEngineResults.h"
#include <string>
#include <vector>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class ProteinGroup; }
namespace EngineLayer { class ProteinParsimonyEngine; }


namespace EngineLayer
{
	class ProteinParsimonyResults : public MetaMorpheusEngineResults
	{
	private:
		std::vector<ProteinGroup*> privateProteinGroups;

	public:
		ProteinParsimonyResults(ProteinParsimonyEngine *proteinAnalysisEngine);

		std::vector<ProteinGroup*> getProteinGroups() const;
		void setProteinGroups(const std::vector<ProteinGroup*> &value);

		std::wstring ToString() override;
	};
}
