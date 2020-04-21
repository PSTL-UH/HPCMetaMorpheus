#include "CommonParameters.h"

using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{

    //CommonParameters::CommonParameters() : CommonParameters(, = DissociationType->HCD, = true, = true, = 3, = 12, = true, = false, = 1, = 5, = 200, = 0.01, = false, = true, = false, = false, = nullptr, = nullptr, = nullptr, = -1, nullptr)
    //{
    //}
    CommonParameters::CommonParameters( CommonParameters *cp): CommonParameters(cp->getTaskDescriptor(),
                                                                                cp->getDissociationType(),
                                                                                cp->getDoPrecursorDeconvolution(),
                                                                                cp->getUseProvidedPrecursorInfo(),
                                                                                cp->getDeconvolutionIntensityRatio(),
                                                                                cp->getDeconvolutionMaxAssumedChargeState(),
                                                                                cp->getReportAllAmbiguity(),
                                                                                cp->getAddCompIons(),
                                                                                cp->getTotalPartitions(),
                                                                                cp->getScoreCutoff(),
                                                                                cp->getTopNpeaks(),
                                                                                cp->getMinRatio(),
                                                                                cp->getTrimMs1Peaks(),
                                                                                cp->getTrimMsMsPeaks(), 
                                                                                cp->getUseDeltaScore(),
                                                                                cp->getCalculateEValue(),
                                                                                cp->getProductMassTolerance(),
                                                                                cp->getPrecursorMassTolerance(),
                                                                                cp->getDeconvolutionMassTolerance(),
                                                                                cp->getMaxThreadsToUsePerFile(),
                                                                                cp->getDigestionParams(),
                                                                                cp->getListOfModsVariable(),
                                                                                cp->getListOfModsFixed(),
                                                                                cp->getQValueOutputFilter(),
                                                                                cp->getAssumeOrphanPeaksAreZ1Fragments(),
                                                                                cp->getMaxHeterozygousVariants(),
                                                                                cp->getMinVariantDepth() )
    {
    }
    CommonParameters::CommonParameters (CommonParameters *cp, std::optional<FragmentationTerminus> terminus,
                                        std::optional<bool> addCompIons ): CommonParameters (cp)
    {
        if (terminus.has_value() ) {
            DigestionParams *dp_org = getDigestionParams();
            DigestionParams *dp = new DigestionParams(dp_org->getProtease()->getName(),
                                                      dp_org->getMaxMissedCleavages(),
                                                      dp_org->getMinPeptideLength(),
                                                      dp_org->getMaxPeptideLength(),
                                                      dp_org->getMaxModificationIsoforms(),
                                                      dp_org->getInitiatorMethionineBehavior(),
                                                      dp_org->getMaxModsForPeptide(),
                                                      dp_org->getSearchModeType(),
                                                      terminus.value() );
            setDigestionParams(dp);
        }
        if ( addCompIons.has_value() ) {
            setAddCompIons(addCompIons.value() );
        }
    }
    
	CommonParameters::CommonParameters(const std::string &taskDescriptor,
                                           DissociationType dissociationType,
                                           bool doPrecursorDeconvolution,
                                           bool useProvidedPrecursorInfo,
                                           double deconvolutionIntensityRatio,
                                           int deconvolutionMaxAssumedChargeState,
                                           bool reportAllAmbiguity,
                                           bool addCompIons,
                                           int totalPartitions,
                                           double scoreCutoff,
                                           int topNpeaks,
                                           double minRatio,
                                           bool trimMs1Peaks,
                                           bool trimMsMsPeaks, 
                                           bool useDeltaScore,
                                           bool calculateEValue,
                                           Tolerance *productMassTolerance,
                                           Tolerance *precursorMassTolerance,
                                           Tolerance *deconvolutionMassTolerance,
                                           int maxThreadsToUsePerFile,
                                           DigestionParams *digestionParams,
                                           std::vector<std::tuple<std::string, std::string>> *listOfModsVariable,
                                           std::vector<std::tuple<std::string, std::string>> *listOfModsFixed,
                                           double qValueOutputFilter,
                                           bool assumeOrphanPeaksAreZ1Fragments,
                                           int maxHeterozygousVariants,
                                           int minVariantDepth)
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

		//setMaxThreadsToUsePerFile(maxThreadsToUsePerFile == -1 ? Environment::ProcessorCount > 1 ?
                //                          Environment::ProcessorCount - 1 : 1 : maxThreadsToUsePerFile);
		setMaxThreadsToUsePerFile(maxThreadsToUsePerFile == -1 ? 4 : maxThreadsToUsePerFile);
                
		PpmTolerance tempVar(20);
		setProductMassTolerance((productMassTolerance != nullptr) ? productMassTolerance : &tempVar);
		PpmTolerance tempVar2(5);
		setPrecursorMassTolerance((precursorMassTolerance != nullptr) ? precursorMassTolerance : &tempVar2);
		PpmTolerance tempVar3(4);
		setDeconvolutionMassTolerance((deconvolutionMassTolerance != nullptr) ? deconvolutionMassTolerance : &tempVar3);
                if ( digestionParams != nullptr ) {
                    setDigestionParams(digestionParams);
                }
                else {
                    DigestionParams *tempVar4 = new DigestionParams("trypsin");
                    setDigestionParams(tempVar4);
                }

                if ( listOfModsVariable != nullptr ) {
                    privateListOfModsVariable =   listOfModsVariable;
                }
                else {
                    std::vector<std::tuple<std::string, std::string>> *vs= new std::vector<std::tuple<std::string, std::string>>();
                    vs->push_back(std::make_tuple("Common Variable", "Oxidation on M"));
                    privateListOfModsVariable = vs;
                }

                if ( listOfModsFixed != nullptr ) {
                    privateListOfModsFixed = listOfModsFixed;
                }
                else {
                    std::vector<std::tuple<std::string, std::string>> *vs2 = new std::vector<std::tuple<std::string, std::string>>();
                    vs2->push_back(std::make_tuple("Common Fixed", "Carbamidomethyl on C"));
                    vs2->push_back(std::make_tuple("Common Fixed", "Carbamidomethyl on U"));
                    privateListOfModsFixed = vs2;
                }
                
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
        std::vector<std::tuple<std::string, std::string>>* CommonParameters::getListOfModsVariable() const
        {
            return privateListOfModsVariable;
        }    
        std::vector<std::tuple<std::string, std::string>>* CommonParameters::getListOfModsFixed() const
        {
            return privateListOfModsFixed;
        }        
	bool CommonParameters::getDoPrecursorDeconvolution() const
        {
            return privateDoPrecursorDeconvolution;
        }
	void CommonParameters::setDoPrecursorDeconvolution( bool value )
        {
            privateDoPrecursorDeconvolution = value;
        }

        bool CommonParameters::getUseProvidedPrecursorInfo() const
        {
            return privateUseProvidedPrecursorInfo;
        }
        void CommonParameters::setUseProvidedPrecursorInfo( bool value)
        {
            privateUseProvidedPrecursorInfo = value;
        }

        double CommonParameters::getDeconvolutionIntensityRatio() const
        {
            return privateDeconvolutionIntensityRatio;
        }
        void CommonParameters::setDeconvolutionIntensityRatio( double value)
        {
            privateDeconvolutionIntensityRatio = value;
        }

        int CommonParameters::getDeconvolutionMaxAssumedChargeState() const
        {
            return privateDeconvolutionMaxAssumedChargeState;
        }
        void CommonParameters::setDeconvolutionMaxAssumedChargeState(int value)
        {
            privateDeconvolutionMaxAssumedChargeState = value;
        }

        Tolerance *CommonParameters::getDeconvolutionMassTolerance() const 
        {
            return privateDeconvolutionMassTolerance;
        }
        void CommonParameters::setDeconvolutionMassTolerance(Tolerance *value)
        {
            privateDeconvolutionMassTolerance = value;
        }
        Tolerance *CommonParameters::getPrecursorMassTolerance() const 
        {
            return privatePrecursorMassTolerance;
        }
        void CommonParameters::setPrecursorMassTolerance(Tolerance *value)
        {
            privatePrecursorMassTolerance = value;
        }

        int CommonParameters::getTotalPartitions() const
        {
            return privateTotalPartitions;
        }
        void CommonParameters::setTotalPartitions(int value)
        {
            privateTotalPartitions = value;
        }

        Tolerance *CommonParameters::getProductMassTolerance() const
        {
            return privateProductMassTolerance;
        }
        void CommonParameters::setProductMassTolerance(Tolerance *value)
        {
            // public setter required for calibration task
            privateProductMassTolerance = value;
        }

        bool CommonParameters::getAddCompIons() const
        {
            return privateAddCompIons;
        }
        void CommonParameters::setAddCompIons( bool value)
        {
            privateAddCompIons = value;
        }
    
        double CommonParameters::getScoreCutoff() const
        {
            return privateScoreCutoff;
        }

        void CommonParameters::setScoreCutoff(double value) 
        {
            privateScoreCutoff=value;
        }

        void CommonParameters::getScoreCutoff(double value)
        {
            privateScoreCutoff = value;
        }

        DigestionParams *CommonParameters::getDigestionParams() const
        {
            return privateDigestionParams;
        }
        void CommonParameters::setDigestionParams(DigestionParams *value)
        {
            privateDigestionParams = value;
        }

        bool CommonParameters::getReportAllAmbiguity() const
        {
            return privateReportAllAmbiguity;
        }
        void CommonParameters::setReportAllAmbiguity(bool value)
        {
            privateReportAllAmbiguity = value;
        }

        int CommonParameters::getTopNpeaks() const
        {            
            return privateToppeaks;
        }
        void CommonParameters::setTopNpeaks(int value)
        {
            privateToppeaks = value;
        }

        double CommonParameters::getMinRatio() const
        {
            return privateMinRatio;
        }
        void CommonParameters::setMinRatio( double value)
        {
            privateMinRatio = value;
        }

        bool CommonParameters::getTrimMs1Peaks() const
        {
            return privateTrimMs1Peaks;
        }
        void CommonParameters::setTrimMs1Peaks(bool value)
        {
            privateTrimMs1Peaks = value;
        }

        bool CommonParameters::getTrimMsMsPeaks() const
        {
            return privateTrimMsMsPeaks;
        }
        void CommonParameters::setTrimMsMsPeaks(bool value)
        {
            privateTrimMsMsPeaks = value;
        }

        bool CommonParameters::getUseDeltaScore() const
        {
            return privateUseDeltaScore;
        }
        void CommonParameters::setUseDeltaScore(bool value)
        {
            privateUseDeltaScore = value;
        }

        bool CommonParameters::getCalculateEValue() const
        {
            return privateCalculateEValue;
        }
        void CommonParameters::setCalculateEValue(bool value)
        {
            privateCalculateEValue = value;
        }

        double CommonParameters::getQValueOutputFilter() const
        {
            return privateQValueOutputFilter;
        }
        void CommonParameters::setQValueOutputFilter(double value)
        {
            privateQValueOutputFilter = value;
        }

        DissociationType CommonParameters::getDissociationType() const
        {
            return privateDissociationType;
        }
        void CommonParameters::setDissociationType(DissociationType value)
        {
            privateDissociationType = value;
        }

        bool CommonParameters::getAssumeOrphanPeaksAreZ1Fragments() const
        {
            return privateAssumeOrphanPeaksAreZ1Fragments;
        }
        void CommonParameters::setAssumeOrphanPeaksAreZ1Fragments(bool value)
        {
            privateAssumeOrphanPeaksAreZ1Fragments = value;
        }

        int CommonParameters::getMaxHeterozygousVariants() const
        {
            return privateMaxHeterozygousVariants;
        }
        void CommonParameters::setMaxHeterozygousVariants(int value)
        {
            privateMaxHeterozygousVariants = value;
        }

        int CommonParameters::getMinVariantDepth() const
        {
            return privateMinVariantDepth;
        }
        void CommonParameters::setMinVariantDepth(int value)
        {
            privateMinVariantDepth = value;
        }

}

