#pragma once

#include "../MetaMorpheusTask.h"
#include <string>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_filesystem.h"

#include "XLSearchParameters.h"
#include "../DbForTask.h"
#include "../MyTaskResults.h"

#include "../FileSpecificParameters.h"

#include "../../EngineLayer/MetaMorpheusEngine.h"
using namespace EngineLayer;

#include "../../EngineLayer/CrosslinkSearch/CrosslinkSearchEngine.h"
using namespace EngineLayer::CrosslinkSearch;

#include "../../EngineLayer/Indexing/IndexingEngine.h"
using namespace EngineLayer::Indexing;

#include "../../EngineLayer/FdrAnalysis/FdrAnalysisEngine.h"
using namespace EngineLayer::FdrAnalysis;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace Proteomics::Fragmentation;

#include "MzLibUtil.h"
using namespace MzLibUtil;

namespace TaskLayer
{
    class XLSearchTask : public MetaMorpheusTask
    {
    private:
        TaskLayer::XlSearchParameters *privateXlSearchParameters=nullptr;
        // need access to proteinList for reconstruction of CrosslinkSpectralMatches in the MPI
        // version, hence they are stored as part of the class now, not just internally in the
        // RunSpecific() method
        std::vector<Protein*> proteinList;

        static constexpr int AVG_PSMS_SERIALIZED_SIZE = 1280; //10 lines a 128 bytes.
        
    public:
        XLSearchTask();
        XLSearchTask(std::string tomlFile );
        ~XLSearchTask() {
            delete privateXlSearchParameters;
        }
        
        TaskLayer::XlSearchParameters *getXlSearchParameters() const;
        void setXlSearchParameters(TaskLayer::XlSearchParameters *value);
        void writeTomlConfig(std::string &filename, std::ofstream &tomlFd );
        std::vector<Protein*> getProteinList () const;
        
    protected:
        MyTaskResults *RunSpecific(const std::string &OutputFolder, std::vector<DbForTask*> &dbFilenameList,
                                   std::vector<std::string> &currentRawFileList, const std::string &taskId,
                                   std::vector<FileSpecificParameters*> &fileSettingsList) override;

    private:
        //Calculate the FDR of single peptide FP/TP
        void SingleFDRAnalysis(std::vector<CrosslinkSpectralMatch*> &items, std::vector<std::string> &taskIds);
        
        //Calculate the FDR of crosslinked peptide FP/TP
        void DoCrosslinkFdrAnalysis(std::vector<CrosslinkSpectralMatch*> &csms);
        
        //Generate user defined crosslinker

        void Gather_Psms ( std::vector<CrosslinkSpectralMatch*> &allPsms,
                           std::vector<Protein *> &proteinList,
                           MPI_Comm comm);

    public:
        static Crosslinker *GenerateUserDefinedCrosslinker(TaskLayer::XlSearchParameters *xlSearchParameters);
        
        
        void WritePsmCrossToTsv(std::vector<CrosslinkSpectralMatch*> &items, const std::string &filePath, int writeType);
        
        void WriteCrosslinkToTxtForPercolator(std::vector<CrosslinkSpectralMatch*> &items,
                                              const std::string &outputFolder,
                                              const std::string &fileName,
                                              Crosslinker *crosslinker,
                                              std::vector<std::string> &nestedIds);

        void WritePepXML_xl(std::vector<CrosslinkSpectralMatch*> &items,
                            std::vector<Protein*> &proteinList,
                            const std::string &databasePath,
                            std::vector<Modification*> &variableModifications,
                            std::vector<Modification*> &fixedModifications,
                            std::vector<std::string> &localizeableModificationTypes,
                            const std::string &outputFolder,
                            const std::string &fileName,
                            std::vector<std::string> &nestedIds);
    };
}
