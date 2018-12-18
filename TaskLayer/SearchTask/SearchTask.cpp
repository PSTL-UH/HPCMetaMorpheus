#include "SearchTask.h"
#include "SearchParameters.h"
#include "../../EngineLayer/CommonParameters.h"
#include "../../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/SinglePpmAroundZeroMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/SingleAbsoluteAroundZeroMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/DotMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/IntervalMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/OpenMassDiffAcceptor.h"
#include "../../EngineLayer/MetaMorpheusException.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"
#include "../MyTaskResults.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../../EngineLayer/ProteinScoringAndFdr/FdrCategory.h"
#include "../MyFileManager.h"
#include "../../EngineLayer/GlobalVariables.h"
#include "../../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "PostSearchAnalysisParameters.h"
#include "PostSearchAnalysisTask.h"

using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::Indexing;
using namespace EngineLayer::ModernSearch;
using namespace EngineLayer::NonSpecificEnzymeSearch;
using namespace FlashLFQ;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{

	SearchTask::SearchTask() : MetaMorpheusTask(MyTask::Search)
	{
		EngineLayer::CommonParameters tempVar();
		setCommonParameters(&tempVar);

		TaskLayer::SearchParameters tempVar2();
		setSearchParameters(&tempVar2);
	}

	TaskLayer::SearchParameters *SearchTask::getSearchParameters() const
	{
		return privateSearchParameters;
	}

	void SearchTask::setSearchParameters(TaskLayer::SearchParameters *value)
	{
		privateSearchParameters = value;
	}

	MassDiffAcceptor *SearchTask::GetMassDiffAcceptor(Tolerance *precursorMassTolerance, MassDiffAcceptorType massDiffAcceptorType, const std::wstring &customMdac)
	{
		switch (massDiffAcceptorType)
		{
			case MassDiffAcceptorType::Exact:
				if (dynamic_cast<PpmTolerance*>(precursorMassTolerance) != nullptr)
				{
					return new SinglePpmAroundZeroSearchMode(precursorMassTolerance->Value);
				}
				else
				{
					return new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance->Value);
				}

			case MassDiffAcceptorType::OneMM:
				return new DotMassDiffAcceptor(L"1mm", {0, 1.0029}, precursorMassTolerance);

			case MassDiffAcceptorType::TwoMM:
				return new DotMassDiffAcceptor(L"2mm", {0, 1.0029, 2.0052}, precursorMassTolerance);

			case MassDiffAcceptorType::ThreeMM:
				return new DotMassDiffAcceptor(L"3mm", {0, 1.0029, 2.0052, 3.0077}, precursorMassTolerance);

			case MassDiffAcceptorType::ModOpen:
				return new IntervalMassDiffAcceptor(L"-187andUp", std::vector<DoubleRange*> {new DoubleRange(-187, std::numeric_limits<double>::infinity())});

			case MassDiffAcceptorType::Open:
				return new OpenSearchMode();

			case MassDiffAcceptorType::Custom:
				return ParseSearchMode(customMdac);

			default:
				throw MetaMorpheusException(L"Unknown MassDiffAcceptorType");
		}
	}

	MyTaskResults *SearchTask::RunSpecific(const std::wstring &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::wstring> &currentRawFileList, const std::wstring &taskId, std::vector<FileSpecificParameters*> &fileSettingsList)
	{
		// disable quantification if a .mgf is being used
		if (getSearchParameters()->getDoQuantification() && currentRawFileList.Any([&] (std::any x)
		{
			Path::GetExtension(x).Equals(L".mgf", StringComparison::OrdinalIgnoreCase);
		}))
		{
			getSearchParameters()->setDoQuantification(false);
		}

		std::vector<Modification> variableModifications;
		std::vector<Modification> fixedModifications;
		std::vector<string> localizeableModificationTypes;
		LoadModifications(taskId, variableModifications, fixedModifications, localizeableModificationTypes);

		// load proteins
		std::vector<Protein*> proteinList = LoadProteins(taskId, dbFilenameList, getSearchParameters()->getSearchTarget(), getSearchParameters()->getDecoyType(), localizeableModificationTypes, getCommonParameters());

		// write prose settings
		ProseCreatedWhileRunning->append(L"The following search settings were used: ");
		ProseCreatedWhileRunning->append(L"protease = " + getCommonParameters()->getDigestionParams()->Protease + L"; ");
		ProseCreatedWhileRunning->append(L"maximum missed cleavages = " + getCommonParameters()->getDigestionParams()->MaxMissedCleavages + L"; ");
		ProseCreatedWhileRunning->append(L"minimum peptide length = " + getCommonParameters()->getDigestionParams()->MinPeptideLength + L"; ");
		ProseCreatedWhileRunning->append(getCommonParameters()->getDigestionParams()->MaxPeptideLength == std::numeric_limits<int>::max() ? L"maximum peptide length = unspecified; " : L"maximum peptide length = " + getCommonParameters()->getDigestionParams()->MaxPeptideLength + L"; ");
		ProseCreatedWhileRunning->append(L"initiator methionine behavior = " + getCommonParameters()->getDigestionParams()->InitiatorMethionineBehavior + L"; ");
		ProseCreatedWhileRunning->append(L"fixed modifications = " + std::wstring::Join(L", ", fixedModifications->Select([&] (std::any m)
		{
			m::IdWithMotif;
		})) + L"; ");
		ProseCreatedWhileRunning->append(L"variable modifications = " + std::wstring::Join(L", ", variableModifications->Select([&] (std::any m)
		{
			m::IdWithMotif;
		})) + L"; ");
		ProseCreatedWhileRunning->append(L"max mods per peptide = " + getCommonParameters()->getDigestionParams()->MaxModsForPeptide + L"; ");
		ProseCreatedWhileRunning->append(L"max modification isoforms = " + getCommonParameters()->getDigestionParams()->MaxModificationIsoforms + L"; ");
		ProseCreatedWhileRunning->append(L"precursor mass tolerance = " + getCommonParameters()->getPrecursorMassTolerance() + L"; ");
		ProseCreatedWhileRunning->append(L"product mass tolerance = " + getCommonParameters()->getProductMassTolerance() + L"; ");
		ProseCreatedWhileRunning->append(L"report PSM ambiguity = " + StringHelper::toString(getCommonParameters()->getReportAllAmbiguity()) + L". ");
		ProseCreatedWhileRunning->append(L"The combined search database contained " + proteinList.size()([&] (std::any p)
		{
			!p::IsDecoy;
		}) + L" non-decoy protein entries including " + proteinList.size()([&] (std::any p)
		{
			p::IsContaminant;
		}) + L" contaminant sequences. ");

		// start the search task
		MyTaskResults = new MyTaskResults(this);
		std::vector<PeptideSpectralMatch*> allPsms;

		//generate an array to store category specific fdr values (for speedy semi/nonspecific searches)
		int numFdrCategories = static_cast<int>(Enum::GetValues(FdrCategory::typeid)->Cast<FdrCategory>().Last() + 1); //+1 because it starts at zero
		std::vector<std::vector<PeptideSpectralMatch*>> allCategorySpecificPsms(numFdrCategories);
		for (int i = 0; i < numFdrCategories; i++)
		{
			allCategorySpecificPsms[i] = std::vector<PeptideSpectralMatch*>();
		}

		FlashLfqResults *flashLfqResults = nullptr;

		MyFileManager *myFileManager = new MyFileManager(getSearchParameters()->getDisposeOfFileWhenDone());

		auto fileSpecificCommonParams = fileSettingsList.Select([&] (std::any b)
		{
			SetAllFileSpecificCommonParams(getCommonParameters(), b);
		});

		int completedFiles = 0;
		std::any indexLock = std::any();
		std::any psmLock = std::any();

		Status(L"Searching files...", taskId);
		Status(L"Searching files...", std::vector<std::wstring> {taskId, L"Individual Spectra Files"});

		std::unordered_map<std::wstring, std::vector<int>> numMs2SpectraPerFile;
		for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.size(); spectraFileIndex++)
		{
			if (GlobalVariables::getStopLoops())
			{
				break;
			}

			auto origDataFile = currentRawFileList[spectraFileIndex];

			// mark the file as in-progress
			StartingDataFile(origDataFile, std::vector<std::wstring> {taskId, L"Individual Spectra Files", origDataFile});

			EngineLayer::CommonParameters *combinedParams = SetAllFileSpecificCommonParams(getCommonParameters(), fileSettingsList[spectraFileIndex]);

			MassDiffAcceptor *massDiffAcceptor = GetMassDiffAcceptor(combinedParams->getPrecursorMassTolerance(), getSearchParameters()->getMassDiffAcceptorType(), getSearchParameters()->getCustomMdac());

			auto thisId = std::vector<std::wstring> {taskId, L"Individual Spectra Files", origDataFile};
			NewCollection(FileSystem::getFileName(origDataFile), thisId);
			Status(L"Loading spectra file...", thisId);
			MsDataFile *myMsDataFile = myFileManager->LoadFile(origDataFile, std::make_optional(combinedParams->getTopNpeaks()), std::make_optional(combinedParams->getMinRatio()), combinedParams->getTrimMs1Peaks(), combinedParams->getTrimMsMsPeaks(), combinedParams);
			Status(L"Getting ms2 scans...", thisId);
			std::vector<Ms2ScanWithSpecificMass*> arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy([&] (std::any b)
			{
				b::PrecursorMass;
			})->ToArray();
			numMs2SpectraPerFile.emplace(Path::GetFileNameWithoutExtension(origDataFile), std::vector<int> { myMsDataFile->GetAllScansList()->Count([&] (std::any p)
			{
			delete myFileManager;
				return p->MsnOrder == 2;
			}), arrayOfMs2ScansSortedByMass.size() });
			myFileManager->DoneWithFile(origDataFile);

			std::vector<PeptideSpectralMatch*> fileSpecificPsms(arrayOfMs2ScansSortedByMass.size());

			// modern search
			if (getSearchParameters()->getSearchType() == SearchType::Modern)
			{
				for (int currentPartition = 0; currentPartition < combinedParams->getTotalPartitions(); currentPartition++)
				{
					std::vector<PeptideWithSetModifications*> peptideIndex;
					std::vector<Protein*> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.size() / combinedParams->getTotalPartitions(), ((currentPartition + 1) * proteinList.size() / combinedParams->getTotalPartitions()) - (currentPartition * proteinList.size() / combinedParams->getTotalPartitions()));

					Status(L"Getting fragment dictionary...", std::vector<std::wstring> {taskId});
					auto indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, currentPartition, getSearchParameters()->getDecoyType(), combinedParams, getSearchParameters()->getMaxFragmentSize(), false, dbFilenameList.Select([&] (std::any p)
					{
						new FileInfo(p::FilePath);
					}).ToList(), std::vector<std::wstring> {taskId});
					std::vector<std::vector<int>> fragmentIndex;
					std::vector<std::vector<int>> precursorIndex;
					{
						std::lock_guard<std::mutex> lock(indexLock);
						GenerateIndexes(indexEngine, dbFilenameList, peptideIndex, fragmentIndex, precursorIndex, proteinList, GlobalVariables::getAllModsKnown().ToList(), taskId);
					}

					Status(L"Searching files...", taskId);

					ModernSearchEngine tempVar(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition, combinedParams, massDiffAcceptor, getSearchParameters()->getMaximumMassThatFragmentIonScoreIsDoubled(), thisId);
					(&tempVar)->Run();

					ProgressEventArgs tempVar2(100, L"Done with search " + std::to_wstring(currentPartition + 1) + L"/" + std::to_wstring(combinedParams->getTotalPartitions()) + L"!", thisId);
					ReportProgress(&tempVar2);

//C# TO C++ CONVERTER TODO TASK: A 'delete indexEngine' statement was not added since indexEngine was passed to a method or constructor. Handle memory management manually.
				}
			}
			// nonspecific search
			else if (getSearchParameters()->getSearchType() == SearchType::NonSpecific)
			{
				std::vector<std::vector<PeptideSpectralMatch*>> fileSpecificPsmsSeparatedByFdrCategory(numFdrCategories); //generate an array of all possible locals
				for (int i = 0; i < numFdrCategories; i++) //only add if we're using for FDR, else ignore it as null.
				{
					fileSpecificPsmsSeparatedByFdrCategory[i] = std::vector<PeptideSpectralMatch*>(arrayOfMs2ScansSortedByMass.size());
				}

				std::vector<EngineLayer::CommonParameters*> paramsToUse = {combinedParams};
				if (combinedParams->getDigestionParams()->SearchModeType == CleavageSpecificity::Semi) //if semi, we need to do both N and C to hit everything
				{
					paramsToUse.clear();
					std::vector<FragmentationTerminus*> terminiToUse = {FragmentationTerminus::N, FragmentationTerminus::C};
					for (auto terminus : terminiToUse) //set both termini
					{
						paramsToUse.push_back(combinedParams->CloneWithNewTerminus(std::make_optional(terminus)));
					}
				}
				for (auto paramToUse : paramsToUse)
				{
					for (int currentPartition = 0; currentPartition < paramToUse->getTotalPartitions(); currentPartition++)
					{
						std::vector<PeptideWithSetModifications*> peptideIndex;

						std::vector<Protein*> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.size() / paramToUse->getTotalPartitions(), ((currentPartition + 1) * proteinList.size() / paramToUse->getTotalPartitions()) - (currentPartition * proteinList.size() / paramToUse->getTotalPartitions()));

						std::vector<std::vector<int>> fragmentIndex(1);
						std::vector<std::vector<int>> precursorIndex(1);

						Status(L"Getting fragment dictionary...", std::vector<std::wstring> {taskId});
						auto indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, currentPartition, getSearchParameters()->getDecoyType(), paramToUse, getSearchParameters()->getMaxFragmentSize(), true, dbFilenameList.Select([&] (std::any p)
						{
							new FileInfo(p::FilePath);
						}).ToList(), std::vector<std::wstring> {taskId});
						{
							std::lock_guard<std::mutex> lock(indexLock);
							GenerateIndexes(indexEngine, dbFilenameList, peptideIndex, fragmentIndex, precursorIndex, proteinList, GlobalVariables::getAllModsKnown().ToList(), taskId);
						}

						Status(L"Searching files...", taskId);

						NonSpecificEnzymeSearchEngine tempVar3(fileSpecificPsmsSeparatedByFdrCategory, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, precursorIndex, currentPartition, paramToUse, massDiffAcceptor, getSearchParameters()->getMaximumMassThatFragmentIonScoreIsDoubled(), thisId);
						(&tempVar3)->Run();

						ProgressEventArgs tempVar4(100, L"Done with search " + std::to_wstring(currentPartition + 1) + L"/" + std::to_wstring(paramToUse->getTotalPartitions()) + L"!", thisId);
						ReportProgress(&tempVar4);

//C# TO C++ CONVERTER TODO TASK: A 'delete indexEngine' statement was not added since indexEngine was passed to a method or constructor. Handle memory management manually.
					}
				}
				{
					std::lock_guard<std::mutex> lock(psmLock);
					for (int i = 0; i < allCategorySpecificPsms.size(); i++)
					{
						if (allCategorySpecificPsms[i].size() > 0)
						{
							allCategorySpecificPsms[i].insert(allCategorySpecificPsms[i].end(), fileSpecificPsmsSeparatedByFdrCategory[i].begin(), fileSpecificPsmsSeparatedByFdrCategory[i].end());
						}
					}
				}
			}
			// classic search
			else
			{
				Status(L"Starting search...", thisId);
				ClassicSearchEngine tempVar5(fileSpecificPsms, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, massDiffAcceptor, combinedParams, thisId);
				(&tempVar5)->Run();

				ProgressEventArgs tempVar6(100, L"Done with search!", thisId);
				ReportProgress(&tempVar6);
			}

			{
				std::lock_guard<std::mutex> lock(psmLock);
				allPsms.insert(allPsms.end(), fileSpecificPsms.begin(), fileSpecificPsms.end());
			}

			completedFiles++;
			FinishedDataFile(origDataFile, std::vector<std::wstring> {taskId, L"Individual Spectra Files", origDataFile});
			ProgressEventArgs tempVar7(completedFiles / currentRawFileList.size(), L"Searching...", new std::vector<std::wstring> {taskId, L"Individual Spectra Files"});
			ReportProgress(&tempVar7);
		}

		ProgressEventArgs tempVar8(100, L"Done with all searches!", new std::vector<std::wstring> {taskId, L"Individual Spectra Files"});
		ReportProgress(&tempVar8);

		int numNotches = GetNumNotches(getSearchParameters()->getMassDiffAcceptorType(), getSearchParameters()->getCustomMdac());
		//resolve category specific fdrs (for speedy semi and nonspecific
		if (getSearchParameters()->getSearchType() == SearchType::NonSpecific)
		{
			allPsms = NonSpecificEnzymeSearchEngine::ResolveFdrCategorySpecificPsms(allCategorySpecificPsms, numNotches, taskId, getCommonParameters());
		}

		PostSearchAnalysisParameters *parameters = new PostSearchAnalysisParameters();
		parameters->setSearchTaskResults(MyTaskResults);
		parameters->setSearchTaskId(taskId);
		parameters->setSearchParameters(getSearchParameters());
		parameters->setProteinList(proteinList);
		parameters->setAllPsms(allPsms);
		parameters->setFixedModifications(fixedModifications);
		parameters->setVariableModifications(variableModifications);
		parameters->setListOfDigestionParams(std::unordered_set<DigestionParams*>(fileSpecificCommonParams->Select([&] (std::any p)
		{
			p::DigestionParams;
		})));
		parameters->setCurrentRawFileList(currentRawFileList);
		parameters->setMyFileManager(myFileManager);
		parameters->setNumNotches(numNotches);
		parameters->setOutputFolder(OutputFolder);
		parameters->setIndividualResultsOutputFolder(FileSystem::combine(OutputFolder, L"Individual File Results"));
		parameters->setFlashLfqResults(flashLfqResults);
		parameters->setFileSettingsList(fileSettingsList);
		parameters->setNumMs2SpectraPerFile(numMs2SpectraPerFile);
		parameters->setDatabaseFilenameList(dbFilenameList);
		PostSearchAnalysisTask *postProcessing = new PostSearchAnalysisTask();
		postProcessing->setParameters(parameters);
		postProcessing->setCommonParameters(getCommonParameters());

		delete postProcessing;
//C# TO C++ CONVERTER TODO TASK: A 'delete parameters' statement was not added since parameters was assigned to another object. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete myFileManager' statement was not added since myFileManager was assigned to another object. Handle memory management manually.
		return postProcessing->Run();
	}

	int SearchTask::GetNumNotches(MassDiffAcceptorType massDiffAcceptorType, const std::wstring &customMdac)
	{
		switch (massDiffAcceptorType)
		{
			case MassDiffAcceptorType::Exact:
				return 1;
			case MassDiffAcceptorType::OneMM:
				return 2;
			case MassDiffAcceptorType::TwoMM:
				return 3;
			case MassDiffAcceptorType::ThreeMM:
				return 4;
			case MassDiffAcceptorType::ModOpen:
				return 1;
			case MassDiffAcceptorType::Open:
				return 1;
			case MassDiffAcceptorType::Custom:
				return ParseSearchMode(customMdac)->NumNotches;

			default:
				throw MetaMorpheusException(L"Unknown mass difference acceptor type");
		}
	}

	MassDiffAcceptor *SearchTask::ParseSearchMode(const std::wstring &text)
	{
		MassDiffAcceptor *massDiffAcceptor = nullptr;

		try
		{
			auto split = StringHelper::split(text, L' ');

//C# TO C++ CONVERTER NOTE: The following 'switch' operated on a string variable and was converted to C++ 'if-else' logic:
//			switch (split[1])
//ORIGINAL LINE: case "dot":
			if (split[1] == L"dot")
			{
					std::vector<double> massShifts = Array::ConvertAll(StringHelper::split(split[4], L','), double::Parse);
					std::wstring newString = StringHelper::replace(split[2], L"�", L"");
					double toleranceValue = std::stod(newString, CultureInfo::InvariantCulture);
//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
					if (split[3].ToUpperInvariant() == L"PPM")
					{
						PpmTolerance tempVar(toleranceValue);
						massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, &tempVar);
					}
//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
					else if (split[3].ToUpperInvariant() == L"DA")
					{
						AbsoluteTolerance tempVar2(toleranceValue);
						massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, &tempVar2);
					}


			}
//ORIGINAL LINE: case "interval":
			else if (split[1] == L"interval")
			{
					std::vector<DoubleRange*> doubleRanges = Array::ConvertAll(StringHelper::split(split[2], L';'), [&] (std::any b)
					{
						new DoubleRange(std::stod(b->Trim(std::vector<wchar_t> {L'[', L']'})->Split(L',')[0], CultureInfo::InvariantCulture), std::stod(b->Trim(std::vector<wchar_t> {L'[', L']'})->Split(L',')[1], CultureInfo::InvariantCulture));
					});
					massDiffAcceptor = new IntervalMassDiffAcceptor(split[0], doubleRanges);

			}
//ORIGINAL LINE: case "OpenSearch":
			else if (split[1] == L"OpenSearch")
			{
					massDiffAcceptor = new OpenSearchMode();

			}
//ORIGINAL LINE: case "daltonsAroundZero":
			else if (split[1] == L"daltonsAroundZero")
			{
					massDiffAcceptor = new SingleAbsoluteAroundZeroSearchMode(std::stod(split[2], CultureInfo::InvariantCulture));

			}
//ORIGINAL LINE: case "ppmAroundZero":
			else if (split[1] == L"ppmAroundZero")
			{
					massDiffAcceptor = new SinglePpmAroundZeroSearchMode(std::stod(split[2], CultureInfo::InvariantCulture));

			}
			else
			{
					delete massDiffAcceptor;
					throw MetaMorpheusException(L"Unrecognized search mode type: " + split[1]);
			}
		}
		catch (const std::runtime_error &e)
		{
			delete massDiffAcceptor;
			throw MetaMorpheusException(L"Could not parse search mode string: " + e.what());
		}

//C# TO C++ CONVERTER TODO TASK: A 'delete massDiffAcceptor' statement was not added since massDiffAcceptor was used in a 'return' or 'throw' statement.
		return massDiffAcceptor;
	}
}
