#include "GPTMDParameters.h"
#include "../../EngineLayer/GlobalVariables.h"


namespace EngineLayer
{

	GptmdParameters::GptmdParameters()
	{
		ListOfModsGptmd = GlobalVariables::getAllModsKnown().Where([&] (std::any b)
		{
			return b::ModificationType->Equals("Common Artifact") || b::ModificationType->Equals("Common Biological") || b::ModificationType->Equals("Metal") || b::ModificationType->Equals("Less Common");
		})->Select([&] (std::any b)
		{
			(b::ModificationType, b::IdWithMotif);
		}).ToList();
	}

	private *List < GptmdParameters::(std::string, std::string)
	{
		get;;

	get;
		set;;

	set;
	}
}
}
