#include "FdrAnalysisResults.h"
#include "FdrAnalysisEngine.h"

namespace EngineLayer
{
	namespace FdrAnalysis
	{

		FdrAnalysisResults::FdrAnalysisResults(FdrAnalysisEngine *s) : MetaMorpheusEngineResults(s)
		{
			setDeltaScoreImprovement(false);
		}

		int FdrAnalysisResults::getPsmsWithin1PercentFdr() const
		{
			return privatePsmsWithin1PercentFdr;
		}

		void FdrAnalysisResults::setPsmsWithin1PercentFdr(int value)
		{
			privatePsmsWithin1PercentFdr = value;
		}

		bool FdrAnalysisResults::getDeltaScoreImprovement() const
		{
			return privateDeltaScoreImprovement;
		}

		void FdrAnalysisResults::setDeltaScoreImprovement(bool value)
		{
			privateDeltaScoreImprovement = value;
		}

		std::wstring FdrAnalysisResults::ToString()
		{
			auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->appendLine(MetaMorpheusEngineResults::ToString());
			sb->appendLine(L"PSMs within 1% fdr: " + std::to_wstring(getPsmsWithin1PercentFdr()));
			sb->appendLine(L"Delta Score Used for FDR Analysis: " + StringHelper::toString(getDeltaScoreImprovement()));

			delete sb;
			return sb->toString();
		}
	}
}
