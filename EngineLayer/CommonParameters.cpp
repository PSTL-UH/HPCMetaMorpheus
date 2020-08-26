#include "CommonParameters.h"
#include "Proteomics/ProteolyticDigestion/CleavageSpecificity.h"
#include "stringhelper.h"

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
                
		auto tempVar = new PpmTolerance(20);
		setProductMassTolerance((productMassTolerance != nullptr) ? productMassTolerance : tempVar);

		auto tempVar2 = new PpmTolerance(5);
		setPrecursorMassTolerance((precursorMassTolerance != nullptr) ? precursorMassTolerance : tempVar2);

		auto tempVar3 = new PpmTolerance(4);
		setDeconvolutionMassTolerance((deconvolutionMassTolerance != nullptr) ? deconvolutionMassTolerance : tempVar3);
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

        CommonParameters::CommonParameters(std::string  &tomlFile, const std::string &taskDescriptor)
        {
            CommonParameters();
            
            setTaskDescriptor(taskDescriptor);            
            Toml trw;
            toml::Value toml_value = trw.tomlReadFile(tomlFile);
            toml::Value* fileParameters = trw.getValue(toml_value, "CommonParameters");
            toml::Table tomlTable = fileParameters->as<toml::Table>();

            // parse toml file and set the values
            for (auto const& keyValuePair : tomlTable)
            {
                // we're using the name of the variable here and not a fixed string
                // in case the variable name changes at some point
                if (keyValuePair.first == "MaxThreadsToUsePerFile") {
                    setMaxThreadsToUsePerFile (keyValuePair.second.as<int>() );
                }
                else if (keyValuePair.first == "ListOfModsFixed" ) {
                    // private static List<(string, string)> GetModsFromString(string value)
                    // {
                    //      return value.Split(new string[] { "\t\t" }, StringSplitOptions.RemoveEmptyEntries).Select(
                    //                                    b => (b.Split('\t').First(), b.Split('\t').Last())).ToList();
                    //}
                    std::string tempstr = keyValuePair.second.as<std::string>();
                    std::vector<std::string> svec = StringHelper::split (tempstr, "\t\t" );
                    auto vs= new std::vector<std::tuple<std::string, std::string>>();
                    for ( auto s: svec ) {
                        std::vector<std::string> st = StringHelper::split(s, '\t' );
                        vs->push_back(std::make_tuple(st[0], st[1]));
                    }
                    privateListOfModsFixed = vs;                    
                }
                else if (keyValuePair.first == "ListOfModsVariable" ) {
                    std::string tempstr = keyValuePair.second.as<std::string>();
                    std::vector<std::string> svec = StringHelper::split (tempstr, "\t\t" );
                    auto vs= new std::vector<std::tuple<std::string, std::string>>();
                    for ( auto s: svec ) {
                        std::vector<std::string> st = StringHelper::split(s, '\t' );
                        vs->push_back(std::make_tuple(st[0], st[1]));
                    }
                    privateListOfModsVariable = vs;                    
                }
                else if (keyValuePair.first == "DoPrecursorDeconvolution" ) {
                    setDoPrecursorDeconvolution(keyValuePair.second.as<bool>() );
                }
                else if (keyValuePair.first == "UseProvidedPrecursorInfo" ) {
                    setUseProvidedPrecursorInfo(keyValuePair.second.as<bool>() );
                }
                else if (keyValuePair.first == "DeconvolutionIntensityRatio" ) {
                    setDeconvolutionIntensityRatio(keyValuePair.second.as<double>() );
                }
                else if (keyValuePair.first == "DeconvolutionMaxAssumedChargeState" ) {
                    setDeconvolutionMaxAssumedChargeState(keyValuePair.second.as<int>() );
                }
                else if (keyValuePair.first == "DeconvolutionMassTolerance" ) {
                    Tolerance* t = Tolerance::ParseToleranceString(keyValuePair.second.as<std::string>());
                    setDeconvolutionMassTolerance(t);
                }
                else if (keyValuePair.first == "TotalPartitions" ) {
                    setTotalPartitions(keyValuePair.second.as<int>() );
                }
                else if (keyValuePair.first == "ProductMassTolerance") {
                    Tolerance* t = Tolerance::ParseToleranceString(keyValuePair.second.as<std::string>());
                    setProductMassTolerance(t);
                }
                else if (keyValuePair.first == "PrecursorMassTolerance") {
                    Tolerance* t = Tolerance::ParseToleranceString(keyValuePair.second.as<std::string>());
                    setPrecursorMassTolerance(t);
                }
                else if (keyValuePair.first == "AddCompIons" ) {
                    setAddCompIons(keyValuePair.second.as<bool>() );
                }
                else if (keyValuePair.first == "ScoreCutoff" ) {
                    setScoreCutoff(keyValuePair.second.as<double>() );
                }
                else if (keyValuePair.first == "ReportAllAmbiguity" ) {
                    setReportAllAmbiguity(keyValuePair.second.as<bool>() );
                }
                else if (keyValuePair.first == "TopNpeaks" ) {
                    setTopNpeaks(keyValuePair.second.as<int>() );
                }
                else if (keyValuePair.first == "MinRatio" ) {
                    setMinRatio(keyValuePair.second.as<double>() );
                }
                else if (keyValuePair.first == "TrimMs1Peaks" ) {
                    setTrimMs1Peaks(keyValuePair.second.as<bool>() );
                }
                else if (keyValuePair.first == "TrimMsMsPeaks" ) {
                    setTrimMsMsPeaks(keyValuePair.second.as<bool>() );
                }
                else if (keyValuePair.first == "UseDeltaScore" ) {
                    setUseDeltaScore(keyValuePair.second.as<bool>() );
                }
                else if (keyValuePair.first == "CalculateEValue" ) {
                    setCalculateEValue(keyValuePair.second.as<bool>() );
                }
                else if (keyValuePair.first == "QValueOutputFilter" ) {
                    setQValueOutputFilter(keyValuePair.second.as<double>() );
                }                
                else if (keyValuePair.first == "DissociationType") {
                    std::string type = keyValuePair.second.as<std::string>();
                    auto d = MassSpectrometry::GetDissocationType::GetDissocationTypeFromString(type);
                    setDissociationType(d);
                }
                else if (keyValuePair.first == "AssumeOrphanPeaksAreZ1Fragments" ) {
                    setAssumeOrphanPeaksAreZ1Fragments(keyValuePair.second.as<bool>() );
                }
                else if (keyValuePair.first == "MaxHeterozygousVariants" ) {
                    setMaxHeterozygousVariants(keyValuePair.second.as<int>() );
                }
                else if (keyValuePair.first == "MinVariantDepth" ) {
                    setMinVariantDepth(keyValuePair.second.as<int>() );
                }
                
            }


            if ( privateListOfModsVariable == nullptr ) {
                auto vs= new std::vector<std::tuple<std::string, std::string>>();
                vs->push_back(std::make_tuple("Common Variable", "Oxidation on M"));
                privateListOfModsVariable = vs;
            }
            
            if ( privateListOfModsFixed == nullptr ) {
                auto vs = new std::vector<std::tuple<std::string, std::string>>();
                vs->push_back(std::make_tuple("Common Fixed", "Carbamidomethyl on C"));
                vs->push_back(std::make_tuple("Common Fixed", "Carbamidomethyl on U"));
                privateListOfModsFixed = vs;
            }
            
            
            //And now DigestionParams
            toml::Value* fileParameters2 = trw.getValue(toml_value, "CommonParameters.DigestionParams");
            toml::Table tomlTable2 = fileParameters2->as<toml::Table>();
            DigestionParams *dp = new DigestionParams("trypsin");
            
            // parse toml file and set the values
            for (auto const& keyValuePair : tomlTable2)
            {
                if ( keyValuePair.first == "MaxMissedCleavages" ) {
                    dp->setMaxMissedCleavages(keyValuePair.second.as<int>());
                }
                else if (keyValuePair.first == "InitiatorMethionineBehavior" ) {
                    auto var = keyValuePair.second.as<std::string>();
                    dp->setInitiatorMethionineBehavior(InitiatorMethionineBehaviorFromString(var));
                }
                else if (keyValuePair.first == "MinPeptideLength" ) {
                    dp->setMinPeptideLength(keyValuePair.second.as<int>());
                }
                else if (keyValuePair.first == "MaxPeptideLength" ) {
                    dp->setMaxPeptideLength(keyValuePair.second.as<int>());
                }
                else if (keyValuePair.first == "MaxModificationIsoforms" ) {
                    dp->setMaxModificationIsoforms(keyValuePair.second.as<int>());
                }
                else if (keyValuePair.first == "MaxModsForPeptide") {
                    dp->setMaxModsForPeptide(keyValuePair.second.as<int>());
                }
                else if (keyValuePair.first == "Protease" ) {
                    std::string name = keyValuePair.second.as<std::string>();
                    std::vector<DigestionMotif*> dm;
                    //create Protease instance
                    auto p = new Protease(name, Proteomics::ProteolyticDigestion::CleavageSpecificity::Unknown,
                                          "", "", dm);
                    dp->setProtease(p);
                }
                else if (keyValuePair.first == "SearchModeType" ) {
                    auto val = keyValuePair.second.as<std::string>();
                    dp->setSearchModeType( ProteolyticDigestion::CleavageSpecificityExtension::ParseString(val));
                }
                else if (keyValuePair.first == "FragmentationTerminus" ) {
                    auto val = keyValuePair.second.as<std::string>();
                    dp->setFragmentationTerminus(FragmentationTerminusFromString(val));
                }
                else if (keyValuePair.first == "SpecificProtease" ) {
                    std::string name = keyValuePair.second.as<std::string>();
                    std::vector<DigestionMotif*> dm;

                    //create Protease instance
                    auto p = new Protease(name, Proteomics::ProteolyticDigestion::CleavageSpecificity::Unknown,
                                          "", "", dm);
                    dp->setSpecificProtease(p);
                }
            }

            setDigestionParams(dp);
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
            if ( privateDigestionParams != nullptr ) {
                free ( privateDigestionParams);
            }
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

