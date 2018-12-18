#pragma once

using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	enum class FdrCategory
	{
		//Cleavage Specificity
		FullySpecific = 0,
		SemiSpecific = 1,
		NonSpecific = 2,

		//New category here
	};
}
