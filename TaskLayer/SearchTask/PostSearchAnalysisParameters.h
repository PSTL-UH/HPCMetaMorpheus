#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class MyTaskResults; }
namespace TaskLayer { class SearchParameters; }
namespace EngineLayer { class PeptideSpectralMatch; }
namespace TaskLayer { class FileSpecificParameters; }
namespace TaskLayer { class MyFileManager; }
namespace TaskLayer { class DbForTask; }

using namespace EngineLayer;
using namespace FlashLFQ;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{
	class PostSearchAnalysisParameters
	{
	private:
		MyTaskResults *privateSearchTaskResults;
		std::string privateSearchTaskId;
		TaskLayer::SearchParameters *privateSearchParameters;
		std::vector<Protein*> privateProteinList;
		std::vector<Modification*> privateVariableModifications;
		std::unordered_set<DigestionParams*> privateListOfDigestionParams;
		std::vector<PeptideSpectralMatch*> privateAllPsms;
		FlashLfqResults *privateFlashLfqResults;
		std::vector<Modification*> privateFixedModifications;
		int privateNumNotches = 0;
		std::string privateOutputFolder;
		std::string privateIndividualResultsOutputFolder;
		std::vector<FileSpecificParameters*> privateFileSettingsList;
		std::unordered_map<std::string, std::vector<int>> privateNumMs2SpectraPerFile;
		TaskLayer::MyFileManager *privateMyFileManager;
		std::vector<DbForTask*> privateDatabaseFilenameList;
		std::vector<std::string> privateCurrentRawFileList;

	public:
		MyTaskResults *getSearchTaskResults() const;
		void setSearchTaskResults(MyTaskResults *value);
		std::string getSearchTaskId() const;
		void setSearchTaskId(const std::string &value);
		TaskLayer::SearchParameters *getSearchParameters() const;
		void setSearchParameters(TaskLayer::SearchParameters *value);
		std::vector<Protein*> getProteinList() const;
		void setProteinList(const std::vector<Protein*> &value);
		std::vector<Modification*> getVariableModifications() const;
		void setVariableModifications(const std::vector<Modification*> &value);
		std::unordered_set<DigestionParams*> getListOfDigestionParams() const;
		void setListOfDigestionParams(const std::unordered_set<DigestionParams*> &value);
		std::vector<PeptideSpectralMatch*> getAllPsms() const;
		void setAllPsms(const std::vector<PeptideSpectralMatch*> &value);
		FlashLfqResults *getFlashLfqResults() const;
		void setFlashLfqResults(FlashLfqResults *value);
		std::vector<Modification*> getFixedModifications() const;
		void setFixedModifications(const std::vector<Modification*> &value);
		int getNumNotches() const;
		void setNumNotches(int value);
		std::string getOutputFolder() const;
		void setOutputFolder(const std::string &value);
		std::string getIndividualResultsOutputFolder() const;
		void setIndividualResultsOutputFolder(const std::string &value);
		std::vector<FileSpecificParameters*> getFileSettingsList() const;
		void setFileSettingsList(const std::vector<FileSpecificParameters*> &value);
		std::unordered_map<std::string, std::vector<int>> getNumMs2SpectraPerFile() const;
		void setNumMs2SpectraPerFile(const std::unordered_map<std::string, std::vector<int>> &value);
		TaskLayer::MyFileManager *getMyFileManager() const;
		void setMyFileManager(TaskLayer::MyFileManager *value);
		std::vector<DbForTask*> getDatabaseFilenameList() const;
		void setDatabaseFilenameList(const std::vector<DbForTask*> &value);
		std::vector<std::string> getCurrentRawFileList() const;
		void setCurrentRawFileList(const std::vector<std::string> &value);
	};
}
