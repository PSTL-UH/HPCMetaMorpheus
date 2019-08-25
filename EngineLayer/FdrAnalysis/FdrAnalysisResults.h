#pragma once

#include "../MetaMorpheusEngineResults.h"
#include <string>
#include "stringhelper.h"
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { namespace FdrAnalysis { class FdrAnalysisEngine; } }


namespace EngineLayer
{
	namespace FdrAnalysis
	{
		class FdrAnalysisResults : public MetaMorpheusEngineResults
		{
		private:
			int privatePsmsWithin1PercentFdr = 0;
			bool privateDeltaScoreImprovement = false;

		public:
			FdrAnalysisResults(FdrAnalysisEngine *s);

				int getPsmsWithin1PercentFdr() const;
				void setPsmsWithin1PercentFdr(int value);
				bool getDeltaScoreImprovement() const;
				void setDeltaScoreImprovement(bool value);

			std::string ToString();
		};
	}
}
