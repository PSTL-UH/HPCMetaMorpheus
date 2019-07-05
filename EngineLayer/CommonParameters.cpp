#include "CommonParameters.h"

using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

	CommonParameters::CommonParameters() : CommonParameters(, = DissociationType->HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, nullptr)
	{
	}

	CommonParameters::CommonParameters(const std::string &taskDescriptor, DissociationType *dissociationType, bool doPrecursorDeconvolution, bool useProvidedPrecursorInfo, double deconvolutionIntensityRatio, int deconvolutionMaxAssumedChargeState, bool reportAllAmbiguity, bool addCompIons, int totalPartitions, double scoreCutoff, int topNpeaks, double minRatio, bool trimMs1Peaks, bool trimMsMsPeaks, bool useDeltaScore, bool calculateEValue, Tolerance *productMassTolerance, Tolerance *precursorMassTolerance, Tolerance *deconvolutionMassTolerance, int maxThreadsToUsePerFile, DigestionParams *digestionParams, std::vector<(std::string, std::string)*> &listOfModsVariable, std::vector<(std::string, std::string)*> &listOfModsFixed, double qValueOutputFilter, bool assumeOrphanPeaksAreZ1Fragments, int maxHeterozygousVariants, int minVariantDepth)
	{
		setTaskDescriptor(taskDescriptor);
		setDoPrecursorDeconvolution(doPrecursorDeconvolution);
		setUseProvidedPrecursorInfo(useProvidedPrecursorInfo);
		setDeconvolutionIntensityRatio(deconvolutionIntensityRatio);
		setDeconvolutionMaxAssumedChargeState(deconvolutionMaxAssumedChargeState);
		setReportAllAmbiguity(reportAllAmbiguity);
		setAddCompIons(addCompIons);
		setTotalPartitions(totalPartitions);
		setScoreCutoff(scoreCutoff);
		setTopNpeaks(topNpeaks);
		setMinRatio(minRatio);
		setTrimMs1Peaks(trimMs1Peaks);
		setTrimMsMsPeaks(trimMsMsPeaks);
		setUseDeltaScore(useDeltaScore);
		setCalculateEValue(calculateEValue);
		setMaxThreadsToUsePerFile(maxThreadsToUsePerFile == -1 ? Environment::ProcessorCount > 1 ? Environment::ProcessorCount - 1 : 1 : maxThreadsToUsePerFile);

		PpmTolerance tempVar(20);
		setProductMassTolerance((productMassTolerance != nullptr) ? productMassTolerance : &tempVar);
		PpmTolerance tempVar2(5);
		setPrecursorMassTolerance((precursorMassTolerance != nullptr) ? precursorMassTolerance : &tempVar2);
		PpmTolerance tempVar3(4);
		setDeconvolutionMassTolerance((deconvolutionMassTolerance != nullptr) ? deconvolutionMassTolerance : &tempVar3);
		DigestionParams tempVar4();
		setDigestionParams((digestionParams != nullptr) ? digestionParams : &tempVar4);
		ListOfModsVariable = listOfModsVariable ? listOfModsVariable : {("Common Variable", "Oxidation on M")};
		ListOfModsFixed = listOfModsFixed ? listOfModsFixed : {("Common Fixed", "Carbamidomethyl on C"), ("Common Fixed", "Carbamidomethyl on U")};
		setDissociationType(dissociationType);
		setQValueOutputFilter(qValueOutputFilter);

		setAssumeOrphanPeaksAreZ1Fragments(assumeOrphanPeaksAreZ1Fragments);

		setMaxHeterozygousVariants(maxHeterozygousVariants);
		setMinVariantDepth(minVariantDepth);
	}

	std::string CommonParameters::getTaskDescriptor() const
	{
		return privateTaskDescriptor;
	}

	void CommonParameters::setTaskDescriptor(const std::string &value)
	{
		privateTaskDescriptor = value;
	}

	int CommonParameters::getMaxThreadsToUsePerFile() const
	{
		return privateMaxThreadsToUsePerFile;
	}

	void CommonParameters::setMaxThreadsToUsePerFile(int value)
	{
		privateMaxThreadsToUsePerFile = value;
	}

	private *IEnumerable < CommonParameters::(std::string, std::string)
	{
		get;;

	get;
	private:
		set;;

	set;
	}
//C# TO C++ CONVERTER TODO TASK: Local functions are not converted by C# to C++ Converter:
//	public IEnumerable<(string, string)> ListOfModsVariable
	//		{
	//			get;
	//			private set;
	//		}
	bool getDoPrecursorDeconvolution() {get;private set;}
	bool getUseProvidedPrecursorInfo() {get;private set;}
	double getDeconvolutionIntensityRatio() {get;private set;}
	int getDeconvolutionMaxAssumedChargeState() {get;private set;}
	Tolerance getDeconvolutionMassTolerance() {get;private set;}
	int getTotalPartitions() {get;private set;}
	Tolerance getProductMassTolerance() {get;set;} // public setter required for calibration task
	Tolerance getPrecursorMassTolerance() {get;set;} // public setter required for calibration task
	bool getAddCompIons() {get;private set;}
	double getScoreCutoff() {get;private set;}
	DigestionParams getDigestionParams() {get;private set;}
	bool getReportAllAmbiguity() {get;private set;}
	int getTopNpeaks() {get;private set;}
	double getMinRatio() {get;private set;}
	bool getTrimMs1Peaks() {get;private set;}
	bool getTrimMsMsPeaks() {get;private set;}
	bool getUseDeltaScore() {get;private set;}
	bool getCalculateEValue() {get;private set;}
	double getQValueOutputFilter() {get;private set;}
	DissociationType getDissociationType() {get;private set;}
	bool getAssumeOrphanPeaksAreZ1Fragments() {get;private set;}
	int getMaxHeterozygousVariants() {get;private set;}
	int getMinVariantDepth() {get;private set;}

//C# TO C++ CONVERTER TODO TASK: Local functions are not converted by C# to C++ Converter:
//	public CommonParameters Clone()
	//		{
	//			CommonParameters c = new CommonParameters();
	//			foreach (PropertyInfo @property in typeof(CommonParameters).GetProperties())
	//			{
	//				@property.SetValue(c, @property.GetValue(this));
	//			}
	//			return c;
	//		}

//C# TO C++ CONVERTER TODO TASK: Local functions are not converted by C# to C++ Converter:
//	public CommonParameters CloneWithNewTerminus(Nullable<FragmentationTerminus> terminus = nullptr, Nullable<bool> addCompIons = nullptr) //for use with speedy semi-specific searches to get both termini
	//		{
	//			if (terminus == nullptr)
	//			{
	//				terminus = DigestionParams.FragmentationTerminus;
	//			}
	//			if (addCompIons == nullptr)
	//			{
	//				addCompIons = AddCompIons;
	//			}
	//			return new CommonParameters(TaskDescriptor, DissociationType, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, ReportAllAmbiguity, addCompIons.Value, TotalPartitions, ScoreCutoff, TopNpeaks, MinRatio, TrimMs1Peaks, TrimMsMsPeaks, UseDeltaScore, CalculateEValue, ProductMassTolerance, PrecursorMassTolerance, DeconvolutionMassTolerance, MaxThreadsToUsePerFile, new DigestionParams(DigestionParams.Protease.Name, DigestionParams.MaxMissedCleavages, DigestionParams.MinPeptideLength, DigestionParams.MaxPeptideLength, DigestionParams.MaxModificationIsoforms, DigestionParams.InitiatorMethionineBehavior, DigestionParams.MaxModsForPeptide, DigestionParams.SearchModeType, terminus.Value), ListOfModsVariable, ListOfModsFixed, QValueOutputFilter, AssumeOrphanPeaksAreZ1Fragments);
	//		}
}
	}
