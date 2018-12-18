#include "FdrClassifier.h"

using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

	FdrCategory FdrClassifier::GetCleavageSpecificityCategory(CleavageSpecificity *cleavageSpecificity)
	{
		if (cleavageSpecificity == CleavageSpecificity::Full)
		{
			return FdrCategory::FullySpecific;
		}
		else if (cleavageSpecificity == CleavageSpecificity::Semi)
		{
			return FdrCategory::SemiSpecific;
		}
		else if (cleavageSpecificity == CleavageSpecificity::None)
		{
			return FdrCategory::NonSpecific;
		}
		else
		{
			throw NotImplementedException("Cleavage specificity '" + cleavageSpecificity + L"' has not been immplemented for local FDR calculations.");
		}
	}
}
