#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "stringhelper.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class ProteinGroup; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }

using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	class ProteinScoringAndFdrEngine : public MetaMorpheusEngine
	{
	private:
		const std::vector<PeptideSpectralMatch*> NewPsms;
		const bool NoOneHitWonders;
		const bool TreatModPeptidesAsDifferentPeptides;
		const bool MergeIndistinguishableProteinGroups;
		const std::vector<ProteinGroup*> ProteinGroups;

	public:
		ProteinScoringAndFdrEngine(std::vector<ProteinGroup*> &proteinGroups, std::vector<PeptideSpectralMatch*> &newPsms, bool noOneHitWonders, bool treatModPeptidesAsDifferentPeptides, bool mergeIndistinguishableProteinGroups, CommonParameters *commonParameters, std::vector<std::wstring> &nestedIds);

	protected:
		MetaMorpheusEngineResults *RunSpecific() override;

	private:
		static std::wstring StripDecoyIdentifier(const std::wstring &proteinGroupName); //we're keeping only the better scoring protein group for each target/decoy pair. to do that we need to strip decoy from the name temporarily. this is the "top-picked" method

		void ScoreProteinGroups(std::vector<ProteinGroup*> &proteinGroups, std::vector<PeptideSpectralMatch*> &psmList);

		std::vector<ProteinGroup*> DoProteinFdr(std::vector<ProteinGroup*> &proteinGroups);
	};
}
