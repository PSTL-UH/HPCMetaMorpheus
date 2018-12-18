#pragma once

#include <string>
#include <vector>

using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	class CommonParameters
	{
	private:
		std::wstring privateTaskDescriptor;
		int privateMaxThreadsToUsePerFile = 0;
		bool privateDoPrecursorDeconvolution = false;
		bool privateUseProvidedPrecursorInfo = false;
		double privateDeconvolutionIntensityRatio = 0;
		int privateDeconvolutionMaxAssumedChargeState = 0;
		Tolerance *privateDeconvolutionMassTolerance;
		int privateTotalPartitions = 0;
		Tolerance *privateProductMassTolerance;
		Tolerance *privatePrecursorMassTolerance;
		bool privateAddCompIons = false;
		double privateScoreCutoff = 0;
		DigestionParams *privateDigestionParams;
		bool privateReportAllAmbiguity = false;
		int privateTopNpeaks = 0;
		double privateMinRatio = 0;
		bool privateTrimMs1Peaks = false;
		bool privateTrimMsMsPeaks = false;
		bool privateUseDeltaScore = false;
		bool privateCalculateEValue = false;
		double privateQValueOutputFilter = 0;
		DissociationType *privateDissociationType;
		bool privateAssumeOrphanPeaksAreZ1Fragments = false;
		int privateMaxHeterozygousVariants = 0;
		int privateMinVariantDepth = 0;

		// this parameterless constructor needs to exist to read the toml.
		// if you can figure out a way to get rid of it, feel free...
	public:
		CommonParameters();

		CommonParameters(const std::wstring &taskDescriptor = L"", DissociationType *dissociationType = DissociationType->HCD, bool doPrecursorDeconvolution = true, bool useProvidedPrecursorInfo = true, double deconvolutionIntensityRatio = 3, int deconvolutionMaxAssumedChargeState = 12, bool reportAllAmbiguity = true, bool addCompIons = false, int totalPartitions = 1, double scoreCutoff = 5, int topNpeaks = 200, double minRatio = 0.01, bool trimMs1Peaks = false, bool trimMsMsPeaks = true, bool useDeltaScore = false, bool calculateEValue = false, Tolerance *productMassTolerance = nullptr, Tolerance *precursorMassTolerance = nullptr, Tolerance *deconvolutionMassTolerance = nullptr, int maxThreadsToUsePerFile = -1, DigestionParams *digestionParams = nullptr, std::vector<(std::wstring, std::wstring)*> &listOfModsVariable, std::vector<(std::wstring, std::wstring)*> &listOfModsFixed, double qValueOutputFilter = 1.0, bool assumeOrphanPeaksAreZ1Fragments = true, int maxHeterozygousVariants = 4, int minVariantDepth = 1);

		// Notes:
		// 1) Any new property must not be nullable (such as int?) or else if it is null,
		//    the null setting will not be written to a toml
		//    and the default will override (so it's okay ONLY if the default is null)
		// 2) All setters should be private unless necessary

		std::wstring getTaskDescriptor() const;
		void setTaskDescriptor(const std::wstring &value);
		int getMaxThreadsToUsePerFile() const;
		void setMaxThreadsToUsePerFile(int value);
		private *IEnumerable < (std::wstring, std::wstring);
