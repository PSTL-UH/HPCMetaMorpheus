#pragma once

#include "../MetaMorpheusTask.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include "stringhelper.h"
#include "tangible_filesystem.h"

#include "PostSearchAnalysisParameters.h"
#include "../../EngineLayer/ProteinParsimony/ProteinGroup.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../MyTaskResults.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"

#include "../../EngineLayer/EngineLayer.h"
using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;

#include "../../EngineLayer/HistogramAnalysis/BinTreeStructure.h"
using namespace EngineLayer::HistogramAnalysis;

#include "../../EngineLayer/Localization/LocalizationEngineResults.h"
using namespace EngineLayer::Localization;

#include "../../EngineLayer/ModificationAnalysis/ModificationAnalysisResults.h"
using namespace EngineLayer::ModificationAnalysis;

#include "FlashLFQ/FlashLFQResults.h"
using namespace FlashLFQ;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;
//using namespace MathNet::Numerics::Distributions;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{
	class PostSearchAnalysisTask : public MetaMorpheusTask
	{
	private:
		PostSearchAnalysisParameters *privateParameters;
		std::vector<EngineLayer::ProteinGroup*> privateProteinGroups;
		//std::vector<IGrouping<std::string, PeptideSpectralMatch*>*> privatePsmsGroupedByFile;
		std::vector<std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>>*> privatePsmsGroupedByFile;

	public:
		PostSearchAnalysisParameters *getParameters() const;
		void setParameters(PostSearchAnalysisParameters *value);
	private:
		std::vector<EngineLayer::ProteinGroup*> getProteinGroups() const;
		void setProteinGroups(const std::vector<EngineLayer::ProteinGroup*> &value);

		//std::vector<IGrouping<std::string, PeptideSpectralMatch*>*> getPsmsGroupedByFile() const;
		std::vector<std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>>*> getPsmsGroupedByFile() const;

                //void setPsmsGroupedByFile(const std::vector<IGrouping<std::string, PeptideSpectralMatch*>*> &value);
                void setPsmsGroupedByFile(const std::vector<std::unordered_map<std::string, std::vector<PeptideSpectralMatch*>>*> &value);

	public:
		PostSearchAnalysisTask();

		MyTaskResults *Run();

	protected:
		MyTaskResults *RunSpecific(const std::string &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::string> &currentRawFileList, const std::string &taskId, std::vector<FileSpecificParameters*> &fileSettingsList) override;

		/// <summary>
		/// Calculate estimated false-discovery rate (FDR) for peptide spectral matches (PSMs)
		/// </summary>
	private:
		void CalculatePsmFdr();

		void ProteinAnalysis();

		void DoMassDifferenceLocalizationAnalysis();

		void QuantificationAnalysis();

		void HistogramAnalysis();

		void WritePsmResults();

		void WriteProteinResults();

		void WriteQuantificationResults();

		void WritePrunedDatabase();

		static int GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications *peptideWithSetModifications);

		static void WriteTree(BinTreeStructure *myTreeStructure, const std::string &writtenFile);

		void WritePsmsForPercolator(std::vector<PeptideSpectralMatch*> &psmList, const std::string &writtenFileForPercolator, double qValueCutoff);

		void WriteProteinGroupsToTsv(std::vector<EngineLayer::ProteinGroup*> &proteinGroups, const std::string &filePath, std::vector<std::string> &nestedIds, double qValueCutoff);

		void WritePeptideQuantificationResultsToTsv(FlashLfqResults *flashLFQResults, const std::string &outputFolder, const std::string &fileName, std::vector<std::string> &nestedIds);

		void WritePeakQuantificationResultsToTsv(FlashLfqResults *flashLFQResults, const std::string &outputFolder, const std::string &fileName, std::vector<std::string> &nestedIds);
	};
}
