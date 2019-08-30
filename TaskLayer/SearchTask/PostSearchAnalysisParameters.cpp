#include "PostSearchAnalysisParameters.h"
#include "../MyTaskResults.h"
#include "SearchParameters.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../FileSpecificParameters.h"
#include "../MyFileManager.h"
#include "../DbForTask.h"

using namespace EngineLayer;
using namespace FlashLFQ;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{

	MyTaskResults *PostSearchAnalysisParameters::getSearchTaskResults() const
	{
		return privateSearchTaskResults;
	}

	void PostSearchAnalysisParameters::setSearchTaskResults(MyTaskResults *value)
	{
		privateSearchTaskResults = value;
	}

	std::string PostSearchAnalysisParameters::getSearchTaskId() const
	{
		return privateSearchTaskId;
	}

	void PostSearchAnalysisParameters::setSearchTaskId(const std::string &value)
	{
		privateSearchTaskId = value;
	}

	TaskLayer::SearchParameters *PostSearchAnalysisParameters::getSearchParameters() const
	{
		return privateSearchParameters;
	}

	void PostSearchAnalysisParameters::setSearchParameters(TaskLayer::SearchParameters *value)
	{
		privateSearchParameters = value;
	}

	std::vector<Protein*> PostSearchAnalysisParameters::getProteinList() const
	{
		return privateProteinList;
	}

	void PostSearchAnalysisParameters::setProteinList(const std::vector<Protein*> &value)
	{
		privateProteinList = value;
	}

	std::vector<Modification*> PostSearchAnalysisParameters::getVariableModifications() const
	{
		return privateVariableModifications;
	}

	void PostSearchAnalysisParameters::setVariableModifications(const std::vector<Modification*> &value)
	{
		privateVariableModifications = value;
	}

	std::unordered_set<DigestionParams*> PostSearchAnalysisParameters::getListOfDigestionParams() const
	{
		return privateListOfDigestionParams;
	}

	void PostSearchAnalysisParameters::setListOfDigestionParams(const std::unordered_set<DigestionParams*> &value)
	{
		privateListOfDigestionParams = value;
	}

	std::vector<PeptideSpectralMatch*> PostSearchAnalysisParameters::getAllPsms() const
	{
		return privateAllPsms;
	}

	void PostSearchAnalysisParameters::setAllPsms(const std::vector<PeptideSpectralMatch*> &value)
	{
		privateAllPsms = value;
	}

	FlashLfqResults *PostSearchAnalysisParameters::getFlashLfqResults() const
	{
		return privateFlashLfqResults;
	}

	void PostSearchAnalysisParameters::setFlashLfqResults(FlashLfqResults *value)
	{
		privateFlashLfqResults = value;
	}

	std::vector<Modification*> PostSearchAnalysisParameters::getFixedModifications() const
	{
		return privateFixedModifications;
	}

	void PostSearchAnalysisParameters::setFixedModifications(const std::vector<Modification*> &value)
	{
		privateFixedModifications = value;
	}

	int PostSearchAnalysisParameters::getNumNotches() const
	{
		return privateNumNotches;
	}

	void PostSearchAnalysisParameters::setNumNotches(int value)
	{
		privateNumNotches = value;
	}

	std::string PostSearchAnalysisParameters::getOutputFolder() const
	{
		return privateOutputFolder;
	}

	void PostSearchAnalysisParameters::setOutputFolder(const std::string &value)
	{
		privateOutputFolder = value;
	}

	std::string PostSearchAnalysisParameters::getIndividualResultsOutputFolder() const
	{
		return privateIndividualResultsOutputFolder;
	}

	void PostSearchAnalysisParameters::setIndividualResultsOutputFolder(const std::string &value)
	{
		privateIndividualResultsOutputFolder = value;
	}

	std::vector<FileSpecificParameters*> PostSearchAnalysisParameters::getFileSettingsList() const
	{
		return privateFileSettingsList;
	}

	void PostSearchAnalysisParameters::setFileSettingsList(const std::vector<FileSpecificParameters*> &value)
	{
		privateFileSettingsList = value;
	}

	std::unordered_map<std::string, std::vector<int>> PostSearchAnalysisParameters::getNumMs2SpectraPerFile() const
	{
		return privateNumMs2SpectraPerFile;
	}

	void PostSearchAnalysisParameters::setNumMs2SpectraPerFile(const std::unordered_map<std::string, std::vector<int>> &value)
	{
		privateNumMs2SpectraPerFile = value;
	}

	TaskLayer::MyFileManager *PostSearchAnalysisParameters::getMyFileManager() const
	{
		return privateMyFileManager;
	}

	void PostSearchAnalysisParameters::setMyFileManager(TaskLayer::MyFileManager *value)
	{
		privateMyFileManager = value;
	}

	std::vector<DbForTask*> PostSearchAnalysisParameters::getDatabaseFilenameList() const
	{
		return privateDatabaseFilenameList;
	}

	void PostSearchAnalysisParameters::setDatabaseFilenameList(const std::vector<DbForTask*> &value)
	{
		privateDatabaseFilenameList = value;
	}

	std::vector<std::string> PostSearchAnalysisParameters::getCurrentRawFileList() const
	{
		return privateCurrentRawFileList;
	}

	void PostSearchAnalysisParameters::setCurrentRawFileList(const std::vector<std::string> &value)
	{
		privateCurrentRawFileList = value;
	}
}
