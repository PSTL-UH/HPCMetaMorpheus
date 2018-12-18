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

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class PostSearchAnalysisParameters; }
namespace EngineLayer { class ProteinGroup; }
namespace EngineLayer { class PeptideSpectralMatch; }
namespace TaskLayer { class MyTaskResults; }
namespace TaskLayer { class DbForTask; }
namespace TaskLayer { class FileSpecificParameters; }

using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::HistogramAnalysis;
using namespace EngineLayer::Localization;
using namespace EngineLayer::ModificationAnalysis;
using namespace FlashLFQ;
using namespace MassSpectrometry;
using namespace MathNet::Numerics::Distributions;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{
	class PostSearchAnalysisTask : public MetaMorpheusTask
	{
	private:
		PostSearchAnalysisParameters *privateParameters;
		std::vector<EngineLayer::ProteinGroup*> privateProteinGroups;
		std::vector<IGrouping<std::wstring, PeptideSpectralMatch*>*> privatePsmsGroupedByFile;

	public:
		PostSearchAnalysisParameters *getParameters() const;
		void setParameters(PostSearchAnalysisParameters *value);
	private:
		std::vector<EngineLayer::ProteinGroup*> getProteinGroups() const;
		void setProteinGroups(const std::vector<EngineLayer::ProteinGroup*> &value);
		std::vector<IGrouping<std::wstring, PeptideSpectralMatch*>*> getPsmsGroupedByFile() const;
		void setPsmsGroupedByFile(const std::vector<IGrouping<std::wstring, PeptideSpectralMatch*>*> &value);

	public:
		PostSearchAnalysisTask();

		MyTaskResults *Run();

	protected:
		MyTaskResults *RunSpecific(const std::wstring &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::wstring> &currentRawFileList, const std::wstring &taskId, std::vector<FileSpecificParameters*> &fileSettingsList) override;

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

		static void WriteTree(BinTreeStructure *myTreeStructure, const std::wstring &writtenFile);

		void WritePsmsForPercolator(std::vector<PeptideSpectralMatch*> &psmList, const std::wstring &writtenFileForPercolator, double qValueCutoff);

		void WriteProteinGroupsToTsv(std::vector<EngineLayer::ProteinGroup*> &proteinGroups, const std::wstring &filePath, std::vector<std::wstring> &nestedIds, double qValueCutoff);

		void WritePeptideQuantificationResultsToTsv(FlashLfqResults *flashLFQResults, const std::wstring &outputFolder, const std::wstring &fileName, std::vector<std::wstring> &nestedIds);

		void WritePeakQuantificationResultsToTsv(FlashLfqResults *flashLFQResults, const std::wstring &outputFolder, const std::wstring &fileName, std::vector<std::wstring> &nestedIds);
	};
}
