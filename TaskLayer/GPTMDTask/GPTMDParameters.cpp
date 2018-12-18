#include "GPTMDParameters.h"
#include "../../EngineLayer/GlobalVariables.h"


namespace EngineLayer
{

	GptmdParameters::GptmdParameters()
	{
		ListOfModsGptmd = GlobalVariables::getAllModsKnown().Where([&] (std::any b)
		{
			return b::ModificationType->Equals(L"Common Artifact") || b::ModificationType->Equals(L"Common Biological") || b::ModificationType->Equals(L"Metal") || b::ModificationType->Equals(L"Less Common");
		})->Select([&] (std::any b)
		{
			(b::ModificationType, b::IdWithMotif);
		}).ToList();
	}

	private *List < GptmdParameters::(std::wstring, std::wstring)
	{
		get;;

	get;
		set;;

	set;
	}
}
	}
