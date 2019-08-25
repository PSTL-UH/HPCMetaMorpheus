#pragma once

#include "FdrCategory.h"
#include "exceptionhelper.h"

using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	class FdrClassifier final
	{
	public:
		static FdrCategory GetCleavageSpecificityCategory(CleavageSpecificity *cleavageSpecificity);
	};
}
