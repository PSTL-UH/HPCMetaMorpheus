#pragma once

#include <string>
#include <vector>

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;

namespace EngineLayer
{
    class CommonParameters
    {
    private:
        std::string privateTaskDescriptor;
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
        int privateToppeaks = 0;
        double privateMinRatio = 0;
        bool privateTrimMs1Peaks = false;
        bool privateTrimMsMsPeaks = false;
        bool privateUseDeltaScore = false;
        bool privateCalculateEValue = false;
        double privateQValueOutputFilter = 0;
        DissociationType privateDissociationType;
        bool privateAssumeOrphanPeaksAreZ1Fragments = false;
        int privateMaxHeterozygousVariants = 0;
        int privateMinVariantDepth = 0;
        std::vector<std::tuple<std::string, std::string>> *privateListOfModsVariable=nullptr;        
        std::vector<std::tuple<std::string, std::string>> *privateListOfModsFixed=nullptr;        

        // this parameterless constructor needs to exist to read the toml.
        // if you can figure out a way to get rid of it, feel free...
    public:
        //CommonParameters();
        
        CommonParameters(const std::string &taskDescriptor = "",
                         DissociationType dissociationType = DissociationType::HCD,
                         bool doPrecursorDeconvolution = true,
                         bool useProvidedPrecursorInfo = true,
                         double deconvolutionIntensityRatio = 3,
                         int deconvolutionMaxAssumedChargeState = 12,
                         bool reportAllAmbiguity = true,
                         bool addCompIons = false,
                         int totalPartitions = 1,
                         double scoreCutoff = 5,
                         int topNpeaks = 200,
                         double minRatio = 0.01,
                         bool trimMs1Peaks = false,
                         bool trimMsMsPeaks = true,
                         bool useDeltaScore = false,
                         bool calculateEValue = false,
                         Tolerance *productMassTolerance = nullptr,
                         Tolerance *precursorMassTolerance = nullptr,
                         Tolerance *deconvolutionMassTolerance = nullptr,
                         int maxThreadsToUsePerFile = -1,
                         DigestionParams *digestionParams = nullptr,
                         std::vector<std::tuple<std::string, std::string>> *listOfModsVariable=nullptr,
                         std::vector<std::tuple<std::string, std::string>> *listOfModsFixed=nullptr,
                         double qValueOutputFilter = 1.0,
                         bool assumeOrphanPeaksAreZ1Fragments = true,
                         int maxHeterozygousVariants = 4,
                         int minVariantDepth = 1);
        
        // Notes:
        // 1) Any new property must not be nullable (such as int?) or else if it is null,
        //    the null setting will not be written to a toml
        //    and the default will override (so it's okay ONLY if the default is null)
        // 2) All setters should be private unless necessary
        
        std::string getTaskDescriptor() const;
        void setTaskDescriptor(const std::string &value);
        int getMaxThreadsToUsePerFile() const;
        void setMaxThreadsToUsePerFile(int value);
        // private *IEnumerable < (std::string, std::string);
        bool getDoPrecursorDeconvolution() const;
        void setDoPrecursorDeconvolution( bool value );
        bool getUseProvidedPrecursorInfo() const;                
        void setUseProvidedPrecursorInfo( bool value);
        double getDeconvolutionIntensityRatio() const;        
        void setDeconvolutionMaxAssumedChargeState(int value);
        int getDeconvolutionMaxAssumedChargeState() const;
        void setDeconvolutionIntensityRatio( double value);
        Tolerance *getDeconvolutionMassTolerance() const ;
        void setDeconvolutionMassTolerance(Tolerance *value);
        Tolerance *getPrecursorMassTolerance() const ;
        void setPrecursorMassTolerance(Tolerance *value);
        int getTotalPartitions() const;                
        void setTotalPartitions(int value);
        Tolerance *getProductMassTolerance() const;                
        void setProductMassTolerance(Tolerance *value);
        bool getAddCompIons() const;                
        void setAddCompIons( bool value);
        double getScoreCutoff() const;                
        void setScoreCutoff(double);                
        void getScoreCutoff(double value);
        DigestionParams *getDigestionParams() const;                
        void setDigestionParams(DigestionParams *value);
        bool getReportAllAmbiguity() const;                
        void setReportAllAmbiguity(bool value);
        int getTopNpeaks() const;
        void setTopNpeaks(int value);
        double getMinRatio() const;                
        void setMinRatio( double value);
        bool getTrimMs1Peaks() const;                
        void setTrimMs1Peaks(bool value);
        bool getTrimMsMsPeaks() const;                
        void setTrimMsMsPeaks(bool value);
        bool getUseDeltaScore() const;                
        void setUseDeltaScore(bool value);
        bool getCalculateEValue() const;                
        void setCalculateEValue(bool value);
        double getQValueOutputFilter() const;                
        void setQValueOutputFilter(double value);
        DissociationType getDissociationType() const;
        void setDissociationType(DissociationType value);
        bool getAssumeOrphanPeaksAreZ1Fragments() const;                
        void setAssumeOrphanPeaksAreZ1Fragments(bool value);
        int getMaxHeterozygousVariants() const;
        void setMaxHeterozygousVariants(int value);
        int getMinVariantDepth() const;                
        void setMinVariantDepth(int value);                                
    };
}
