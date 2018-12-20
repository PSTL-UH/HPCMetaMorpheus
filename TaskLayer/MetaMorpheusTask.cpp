#include "MetaMorpheusTask.h"
#include "../EngineLayer/CommonParameters.h"
#include "MyTaskResults.h"
#include "../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "FileSpecificParameters.h"
#include "DbForTask.h"
#include "../EngineLayer/MetaMorpheusEngine.h"
#include "../EngineLayer/EventArgs/SingleEngineFinishedEventArgs.h"
#include "../EngineLayer/GlobalVariables.h"
#include "../EngineLayer/MetaMorpheusException.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../EngineLayer/EventArgs/SingleFileEventArgs.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "EventArgs/SingleTaskEventArgs.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Nett;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{

TomlSettings *const MetaMorpheusTask::tomlConfig = TomlSettings::Create([&] (std::any cfg)
{
			cfg::ConfigureType<Tolerance*>([&] (std::any type)
			{
			type::WithConversionFor<TomlString*>([&] (std::any convert)
			{
			convert::FromToml([&] (std::any tmlString)
			{
			Tolerance::ParseToleranceString(tmlString->Value);
			});
			});
			}).ConfigureType<PpmTolerance*>([&] (std::any type)
			{
			type::WithConversionFor<TomlString*>([&] (std::any convert)
			{
			convert::ToToml([&] (std::any custom)
			{
			custom.ToString();
			});
			});
		}).ConfigureType<AbsoluteTolerance*>([&] (std::any type)
		{
			type::WithConversionFor<TomlString*>([&] (std::any convert)
			{
			convert::ToToml([&] (std::any custom)
			{
			custom.ToString();
			});
			});
		}).ConfigureType<Protease*>([&] (std::any type)
		{
			type::WithConversionFor<TomlString*>([&] (std::any convert)
			{
			convert::ToToml([&] (std::any custom)
			{
			custom.ToString();
			}).FromToml([&] (std::any tmlString)
			{
			ProteaseDictionary::Dictionary[tmlString->Value];
		});
			});
		}).ConfigureType<std::vector<std::wstring>>([&] (std::any type)
		{
			type::WithConversionFor<TomlString*>([&] (std::any convert)
			{
			convert::ToToml([&] (std::any custom)
			{
			std::wstring::Join(L"\t", custom);
			}).FromToml([&] (std::any tmlString)
			{
			GetModsTypesFromString(tmlString->Value);
		});
			});
		}).ConfigureType<std::vector<(std::wstring, std::wstring)*>>([&] (std::any type)
		{
			type::WithConversionFor<TomlString*>([&] (std::any convert)
			{
			convert::ToToml([&] (std::any custom)
			{
			std::wstring::Join(L"\t\t", custom->Select([&] (std::any b)
			{
			return b::Item1 + L"\t" + b::Item2;
			}));
			}).FromToml([&] (std::any tmlString)
			{
			GetModsFromString(tmlString->Value);
		});
			});
		});
});

	MetaMorpheusTask::MetaMorpheusTask(MyTask taskType)
	{
		this->setTaskType(taskType);
	}

	MyTask MetaMorpheusTask::getTaskType() const
	{
		return privateTaskType;
	}

	void MetaMorpheusTask::setTaskType(MyTask value)
	{
		privateTaskType = value;
	}

	EngineLayer::CommonParameters *MetaMorpheusTask::getCommonParameters() const
	{
		return privateCommonParameters;
	}

	void MetaMorpheusTask::setCommonParameters(EngineLayer::CommonParameters *value)
	{
		privateCommonParameters = value;
	}

const std::wstring MetaMorpheusTask::IndexFolderName = L"DatabaseIndex";

	std::vector<Ms2ScanWithSpecificMass*> MetaMorpheusTask::GetMs2Scans(MsDataFile *myMSDataFile, const std::wstring &fullFilePath, EngineLayer::CommonParameters *commonParameters)
	{
		auto ms2Scans = myMSDataFile->GetAllScansList().Where([&] (std::any x)
		{
			return x::MsnOrder > 1;
		})->ToArray();
		std::vector<std::vector<Ms2ScanWithSpecificMass*>> scansWithPrecursors(ms2Scans.size());

		ParallelOptions *tempVar = new ParallelOptions();
		tempVar->MaxDegreeOfParallelism = commonParameters->getMaxThreadsToUsePerFile();
		Parallel::ForEach(Partitioner::Create(0, ms2Scans.size()), tempVar, [&] (partitionRange, loopState)
		{
				for (int i = partitionRange::Item1; i < partitionRange::Item2; i++)
				{
					if (GlobalVariables::getStopLoops())
					{
						break;
					}
    
					MsDataScan *ms2scan = ms2Scans[i];
    
					std::vector<(double, int)*> precursors;
					if (ms2scan->OneBasedPrecursorScanNumber.HasValue)
					{
						auto precursorSpectrum = myMSDataFile->GetOneBasedScan(ms2scan->OneBasedPrecursorScanNumber->Value);
    
						try
						{
							ms2scan->RefineSelectedMzAndIntensity(precursorSpectrum->MassSpectrum);
						}
						catch (const MzLibException &ex)
						{
							Warn(L"Could not get precursor ion for MS2 scan #" + ms2scan->OneBasedScanNumber + L"; " + ex->Message);
							continue;
						}
    
						if (ms2scan->SelectedIonMonoisotopicGuessMz.HasValue)
						{
							ms2scan->ComputeMonoisotopicPeakIntensity(precursorSpectrum->MassSpectrum);
						}
    
						if (commonParameters->getDoPrecursorDeconvolution())
						{
							for (auto envelope : ms2scan->GetIsolatedMassesAndCharges(precursorSpectrum->MassSpectrum, 1, commonParameters->getDeconvolutionMaxAssumedChargeState(), commonParameters->getDeconvolutionMassTolerance()->Value, commonParameters->getDeconvolutionIntensityRatio()))
							{
								auto monoPeakMz = envelope->monoisotopicMass.ToMz(envelope->charge);
								precursors.Add((monoPeakMz, envelope->charge));
							}
						}
					}
    
					if (commonParameters->getUseProvidedPrecursorInfo() && ms2scan->SelectedIonChargeStateGuess.HasValue)
					{
						auto precursorCharge = ms2scan->SelectedIonChargeStateGuess->Value;
						if (ms2scan->SelectedIonMonoisotopicGuessMz.HasValue)
						{
							double precursorMZ = ms2scan->SelectedIonMonoisotopicGuessMz->Value;
							if (!precursors.Any([&] (std::any b)
							{
				commonParameters->getDeconvolutionMassTolerance()->Within(precursorMZ.ToMass(precursorCharge), b::Item1->ToMass(b::Item2));
							}))
							{
								precursors.Add((precursorMZ, precursorCharge));
							}
						}
						else
						{
							double precursorMZ = ms2scan->SelectedIonMZ->Value;
							if (!precursors.Any([&] (std::any b)
							{
				commonParameters->getDeconvolutionMassTolerance()->Within(precursorMZ.ToMass(precursorCharge), b::Item1->ToMass(b::Item2));
							}))
							{
								precursors.Add((precursorMZ, precursorCharge));
							}
						}
					}
    
					scansWithPrecursors[i] = std::vector<Ms2ScanWithSpecificMass*>();
					std::vector<IsotopicEnvelope*> neutralExperimentalFragments = Ms2ScanWithSpecificMass::GetNeutralExperimentalFragments(ms2scan, commonParameters);
    
					for (auto precursor : precursors)
					{
						Ms2ScanWithSpecificMass tempVar2(ms2scan, precursor->Item1, precursor->Item2, fullFilePath, commonParameters, neutralExperimentalFragments);
						scansWithPrecursors[i].Add(&tempVar2);
					}
				}
		});

		return scansWithPrecursors.SelectMany([&] (std::any p)
		{
			return p;
		});
	}

	EngineLayer::CommonParameters *MetaMorpheusTask::SetAllFileSpecificCommonParams(EngineLayer::CommonParameters *commonParams, FileSpecificParameters *fileSpecificParams)
	{
		if (fileSpecificParams == nullptr)
		{
			return commonParams;
		}

		// set file-specific digestion parameters
		Protease tempVar = fileSpecificParams.getProtease();
		Protease *protease = (tempVar != nullptr) ? tempVar : commonParams->getDigestionParams()->Protease;
		Nullable<int> tempVar2 = fileSpecificParams.getMinPeptideLength();
		int minPeptideLength = tempVar2 ? tempVar2 : commonParams->getDigestionParams()->MinPeptideLength;
		Nullable<int> tempVar3 = fileSpecificParams.getMaxPeptideLength();
		int maxPeptideLength = tempVar3 ? tempVar3 : commonParams->getDigestionParams()->MaxPeptideLength;
		Nullable<int> tempVar4 = fileSpecificParams.getMaxMissedCleavages();
		int maxMissedCleavages = tempVar4 ? tempVar4 : commonParams->getDigestionParams()->MaxMissedCleavages;
		Nullable<int> tempVar5 = fileSpecificParams.getMaxModsForPeptide();
		int maxModsForPeptide = tempVar5 ? tempVar5 : commonParams->getDigestionParams()->MaxModsForPeptide;
		DigestionParams *fileSpecificDigestionParams = new DigestionParams(protease: protease->Name, maxMissedCleavages: maxMissedCleavages, minPeptideLength: minPeptideLength, maxPeptideLength: maxPeptideLength, maxModsForPeptides: maxModsForPeptide, maxModificationIsoforms: commonParams->getDigestionParams()->MaxModificationIsoforms, initiatorMethionineBehavior: commonParams->getDigestionParams()->InitiatorMethionineBehavior, fragmentationTerminus: commonParams->getDigestionParams()->FragmentationTerminus, searchModeType: commonParams->getDigestionParams()->SearchModeType);

		// set the rest of the file-specific parameters
		Tolerance tempVar6 = fileSpecificParams.getPrecursorMassTolerance();
		Tolerance *precursorMassTolerance = (tempVar6 != nullptr) ? tempVar6 : commonParams->getPrecursorMassTolerance();
		Tolerance tempVar7 = fileSpecificParams.getProductMassTolerance();
		Tolerance *productMassTolerance = (tempVar7 != nullptr) ? tempVar7 : commonParams->getProductMassTolerance();
		Nullable<DissociationType*> tempVar8 = fileSpecificParams.getDissociationType();
		DissociationType *dissociationType = tempVar8 ? tempVar8 : commonParams->getDissociationType();


		EngineLayer::CommonParameters *returnParams = new CommonParameters(commonParams->getTaskDescriptor(), dissociationType, commonParams->getDoPrecursorDeconvolution(), commonParams->getUseProvidedPrecursorInfo(), commonParams->getDeconvolutionIntensityRatio(), commonParams->getDeconvolutionMaxAssumedChargeState(), commonParams->getReportAllAmbiguity(), commonParams->getAddCompIons(), commonParams->getTotalPartitions(), commonParams->getScoreCutoff(), commonParams->getTopNpeaks(), commonParams->getMinRatio(), commonParams->getTrimMs1Peaks(), commonParams->getTrimMsMsPeaks(), commonParams->getUseDeltaScore(), commonParams->getCalculateEValue(), productMassTolerance, precursorMassTolerance, commonParams->getDeconvolutionMassTolerance(), commonParams->getMaxThreadsToUsePerFile(), fileSpecificDigestionParams, commonParams->ListOfModsVariable, commonParams->ListOfModsFixed, commonParams->getQValueOutputFilter(), commonParams->getAssumeOrphanPeaksAreZ1Fragments());

//C# TO C++ CONVERTER TODO TASK: A 'delete returnParams' statement was not added since returnParams was used in a 'return' or 'throw' statement.
		delete fileSpecificDigestionParams;
		return returnParams;
	}

	MyTaskResults *MetaMorpheusTask::RunTask(const std::wstring &output_folder, std::vector<DbForTask*> &currentProteinDbFilenameList, std::vector<std::wstring> &currentRawDataFilepathList, const std::wstring &displayName)
	{
		StartingSingleTask(displayName);

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		auto tomlFileName = FileSystem::combine(Directory::GetParent(output_folder)->ToString(), L"Task Settings", displayName + L"config.toml");
		Toml::WriteFile(this, tomlFileName, tomlConfig);
		FinishedWritingFile(tomlFileName, std::vector<std::wstring> {displayName});

		MetaMorpheusEngine::FinishedSingleEngineHandler->addListener(L"SingleEngineHandlerInTask", [&] (std::any sender, SingleEngineFinishedEventArgs* e) {SingleEngineHandlerInTask(sender, e);});
		try
		{
			auto stopWatch = new Stopwatch();
			stopWatch->Start();

			std::vector<FileSpecificParameters*> fileSettingsList(currentRawDataFilepathList.size());
			for (int i = 0; i < currentRawDataFilepathList.size(); i++)
			{
				if (GlobalVariables::getStopLoops())
				{
					break;
				}
				std::wstring rawFilePath = currentRawDataFilepathList[i];
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				std::wstring directory = Directory::GetParent(rawFilePath)->ToString();
				std::wstring fileSpecificTomlPath = FileSystem::combine(directory, Path::GetFileNameWithoutExtension(rawFilePath)) + L".toml";
				if (FileSystem::fileExists(fileSpecificTomlPath))
				{
					TomlTable *fileSpecificSettings = Toml::ReadFile(fileSpecificTomlPath, tomlConfig);
					try
					{
						fileSettingsList[i] = new FileSpecificParameters(fileSpecificSettings);
					}
					catch (const MetaMorpheusException &e)
					{
						// file-specific toml has already been validated in the GUI when the spectra files were added, so...
						// probably the only time you can get here is if the user modifies the file-specific parameter file in the middle of a run...
						Warn(L"Problem parsing the file-specific toml " + FileSystem::getFileName(fileSpecificTomlPath) + L"; " + e->what() + L"; is the toml from an older version of MetaMorpheus?");
					}
				}
			}

			RunSpecific(output_folder, currentProteinDbFilenameList, currentRawDataFilepathList, displayName, fileSettingsList);
			stopWatch->Stop();
			MyTaskResults->Time = stopWatch->Elapsed;
			auto resultsFileName = FileSystem::combine(output_folder, L"results.txt");
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter file = new StreamWriter(resultsFileName))
			{
				StreamWriter file = StreamWriter(resultsFileName);
				file.WriteLine(L"MetaMorpheus: version " + GlobalVariables::getMetaMorpheusVersion());
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				file.Write(MyTaskResults->ToString());
			}
			FinishedWritingFile(resultsFileName, std::vector<std::wstring> {displayName});
			FinishedSingleTask(displayName);

			delete stopWatch;
		}
		catch (const std::runtime_error &e)
		{
			MetaMorpheusEngine::FinishedSingleEngineHandler->removeListener(L"SingleEngineHandlerInTask");
			auto resultsFileName = FileSystem::combine(output_folder, L"results.txt");
			e.Data->Add(L"folder", output_folder);
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter file = new StreamWriter(resultsFileName))
			{
				StreamWriter file = StreamWriter(resultsFileName);
				file.WriteLine(GlobalVariables::getMetaMorpheusVersion() == L"1.0.0.0" ? L"MetaMorpheus: Not a release version" : L"MetaMorpheus: version " + GlobalVariables::getMetaMorpheusVersion());
				file.WriteLine(SystemInfo::CompleteSystemInfo()); //OS, OS Version, .Net Version, RAM, processor count, MSFileReader .dll versions X3
				file.Write(L"e: " + e);
				file.Write(L"e.Message: " + e.what());
				file.Write(L"e.InnerException: " + e.InnerException);
				file.Write(L"e.Source: " + e.Source);
				file.Write(L"e.StackTrace: " + e.StackTrace);
				file.Write(L"e.TargetSite: " + e.TargetSite);
			}
			throw;
		}

		{
			auto proseFilePath = FileSystem::combine(output_folder, L"prose.txt");
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter file = new StreamWriter(proseFilePath))
			{
				StreamWriter file = StreamWriter(proseFilePath);
				file.Write(L"The data analysis was performed using MetaMorpheus version " + GlobalVariables::getMetaMorpheusVersion() + L", available at " + L"https://github.com/smith-chem-wisc/MetaMorpheus.");
				file.Write(ProseCreatedWhileRunning->toString());
				file.Write(SystemInfo::SystemProse()->Replace(L"\r\n", L"") + L" ");
				file.WriteLine(L"The total time to perform the " + getTaskType() + L" task on " + std::to_wstring(currentRawDataFilepathList.size()) + L" spectra file(s) was " + std::wstring::Format(L"{0:0.00}", MyTaskResults->Time.TotalMinutes) + L" minutes.");
				file.WriteLine();
				file.WriteLine(L"Published works using MetaMorpheus software are encouraged to cite: Solntsev, S. K.; Shortreed, M. R.; Frey, B. L.; Smith, L. M. Enhanced Global Post-translational Modification Discovery with MetaMorpheus. Journal of Proteome Research. 2018, 17 (5), 1844-1851.");

				file.WriteLine();
				file.WriteLine(L"Spectra files: ");
				file.WriteLine(std::wstring::Join(L"\r\n", currentRawDataFilepathList.Select([&] (std::any b)
				{
					return L'\t' + b;
				})));
				file.WriteLine(L"Databases:");
				file.Write(std::wstring::Join(L"\r\n", currentProteinDbFilenameList.Select([&] (std::any b)
				{
					return L'\t' + (b::IsContaminant ? L"Contaminant " : L"") + b::FilePath;
				})));
			}
			FinishedWritingFile(proseFilePath, std::vector<std::wstring> {displayName});
		}

		MetaMorpheusEngine::FinishedSingleEngineHandler->removeListener(L"SingleEngineHandlerInTask");
		return MyTaskResults;
	}

	std::vector<Protein*> MetaMorpheusTask::LoadProteins(const std::wstring &taskId, std::vector<DbForTask*> &dbFilenameList, bool searchTarget, DecoyType *decoyType, std::vector<std::wstring> &localizeableModificationTypes, EngineLayer::CommonParameters *commonParameters)
	{
		Status(L"Loading proteins...", std::vector<std::wstring> {taskId});
		int emptyProteinEntries = 0;
		std::vector<Protein*> proteinList;
		for (auto db : dbFilenameList)
		{
			int emptyProteinEntriesForThisDb = 0;
			Dictionary<std::wstring, Modification*> unknownModifications;
			auto dbProteinList = LoadProteinDb(db->getFilePath(), searchTarget, decoyType, localizeableModificationTypes, db->getIsContaminant(), unknownModifications, emptyProteinEntriesForThisDb, commonParameters);
			proteinList = proteinList.Concat(dbProteinList)->ToList();
			emptyProteinEntries += emptyProteinEntriesForThisDb;
		}
		if (!proteinList.Any())
		{
			Warn(L"Warning: No protein entries were found in the database");
		}
		else if (emptyProteinEntries > 0)
		{
			Warn(L"Warning: " + std::to_wstring(emptyProteinEntries) + L" empty protein entries ignored");
		}
		return proteinList;
	}

	std::vector<Protein*> MetaMorpheusTask::LoadProteinDb(const std::wstring &fileName, bool generateTargets, DecoyType *decoyType, std::vector<std::wstring> &localizeableModificationTypes, bool isContaminant, std::unordered_map<std::wstring, Modification*> &um, int &emptyEntriesCount, EngineLayer::CommonParameters *commonParameters)
	{
		std::vector<std::wstring> dbErrors;
		std::vector<Protein*> proteinList;

//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
		std::wstring theExtension = Path::GetExtension(fileName).ToLowerInvariant();
		bool compressed = StringHelper::endsWith(theExtension, L"gz"); // allows for .bgz and .tgz, too which are used on occasion
//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
		theExtension = compressed ? Path::GetExtension(Path::GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

		if (theExtension == L".fasta" || theExtension == L".fa")
		{
			um.clear();
			proteinList = ProteinDbLoader::LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, ProteinDbLoader::UniprotAccessionRegex, ProteinDbLoader::UniprotFullNameRegex, ProteinDbLoader::UniprotFullNameRegex, ProteinDbLoader::UniprotGeneNameRegex, ProteinDbLoader::UniprotOrganismRegex, dbErrors, commonParameters->getMaxThreadsToUsePerFile());
		}
		else
		{
			std::vector<std::wstring> modTypesToExclude = GlobalVariables::getAllModTypesKnown().Where([&] (std::any b)
			{
				!std::find(localizeableModificationTypes.begin(), localizeableModificationTypes.end(), b) != localizeableModificationTypes.end();
			}).ToList();
			proteinList = ProteinDbLoader::LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables::getAllModsKnown(), isContaminant, modTypesToExclude, um, commonParameters->getMaxThreadsToUsePerFile(), commonParameters->getMaxHeterozygousVariants(), commonParameters->getMinVariantDepth());
		}

		emptyEntriesCount = proteinList.size()([&] (std::any p)
		{
			return p::BaseSequence->Length == 0;
		});
		return proteinList.Where([&] (std::any p)
		{
			return p::BaseSequence->Length > 0;
		}).ToList();
	}

	void MetaMorpheusTask::LoadModifications(const std::wstring &taskId, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<std::wstring> &localizableModificationTypes)
	{
		// load modifications
		Status(L"Loading modifications...", taskId);
		variableModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			getCommonParameters()->ListOfModsVariable->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();
		fixedModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			getCommonParameters()->ListOfModsFixed->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();
		localizableModificationTypes = GlobalVariables::getAllModTypesKnown().ToList();

		auto recognizedVariable = variableModifications.Select([&] (std::any p)
		{
			p::IdWithMotif;
		});
		auto recognizedFixed = fixedModifications.Select([&] (std::any p)
		{
			p::IdWithMotif;
		});
		auto unknownMods = getCommonParameters()->ListOfModsVariable->Select([&] (std::any p)
		{
			p::Item2;
		}).Except(recognizedVariable).ToList();
		unknownMods.AddRange(getCommonParameters()->ListOfModsFixed->Select([&] (std::any p)
		{
			p::Item2;
		}).Except(recognizedFixed));
		for (auto unrecognizedMod : unknownMods)
		{
			Warn(L"Unrecognized mod " + unrecognizedMod + L"; are you using an old .toml?");
		}
	}

	void MetaMorpheusTask::WritePsmsToTsv(std::vector<PeptideSpectralMatch*> &psms, const std::wstring &filePath, IReadOnlyDictionary<std::wstring, int> *modstoWritePruned)
	{
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(filePath))
		{
			StreamWriter output = StreamWriter(filePath);
			output.WriteLine(PeptideSpectralMatch::GetTabSeparatedHeader());
			for (auto psm : psms)
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				output.WriteLine(psm->ToString(modstoWritePruned));
			}
		}
	}

	void MetaMorpheusTask::ReportProgress(ProgressEventArgs *v)
	{
		OutProgressHandler +== nullptr ? nullptr : OutProgressHandler::Invoke(this, v);
	}

	void MetaMorpheusTask::FinishedWritingFile(const std::wstring &path, std::vector<std::wstring> &nestedIDs)
	{
		SingleFileEventArgs tempVar(path, nestedIDs);
		FinishedWritingFileHandler +== nullptr ? nullptr : FinishedWritingFileHandler::Invoke(this, &tempVar);
	}

	void MetaMorpheusTask::StartingDataFile(const std::wstring &v, std::vector<std::wstring> &nestedIDs)
	{
		StringEventArgs tempVar(v, nestedIDs);
		StartingDataFileHandler +== nullptr ? nullptr : StartingDataFileHandler::Invoke(this, &tempVar);
	}

	void MetaMorpheusTask::FinishedDataFile(const std::wstring &v, std::vector<std::wstring> &nestedIDs)
	{
		StringEventArgs tempVar(v, nestedIDs);
		FinishedDataFileHandler +== nullptr ? nullptr : FinishedDataFileHandler::Invoke(this, &tempVar);
	}

	void MetaMorpheusTask::Status(const std::wstring &v, const std::wstring &id)
	{
		StringEventArgs tempVar(v, new std::vector<std::wstring> {id});
		OutLabelStatusHandler +== nullptr ? nullptr : OutLabelStatusHandler::Invoke(this, &tempVar);
	}

	void MetaMorpheusTask::Status(const std::wstring &v, std::vector<std::wstring> &nestedIds)
	{
		StringEventArgs tempVar(v, nestedIds);
		OutLabelStatusHandler +== nullptr ? nullptr : OutLabelStatusHandler::Invoke(this, &tempVar);
	}

	void MetaMorpheusTask::Warn(const std::wstring &v)
	{
		StringEventArgs tempVar(v, nullptr);
		WarnHandler +== nullptr ? nullptr : WarnHandler::Invoke(nullptr, &tempVar);
	}

	void MetaMorpheusTask::Log(const std::wstring &v, std::vector<std::wstring> &nestedIds)
	{
		StringEventArgs tempVar(v, nestedIds);
		LogHandler +== nullptr ? nullptr : LogHandler::Invoke(this, &tempVar);
	}

	void MetaMorpheusTask::NewCollection(const std::wstring &displayName, std::vector<std::wstring> &nestedIds)
	{
		StringEventArgs tempVar(displayName, nestedIds);
		NewCollectionHandler +== nullptr ? nullptr : NewCollectionHandler::Invoke(this, &tempVar);
	}

	std::vector<std::wstring> MetaMorpheusTask::GetModsTypesFromString(const std::wstring &value)
	{
		return value.Split({L"\t"}, StringSplitOptions::RemoveEmptyEntries).ToList();
	}

	void MetaMorpheusTask::SingleEngineHandlerInTask(std::any sender, SingleEngineFinishedEventArgs *e)
	{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MyTaskResults->AddResultText(e->ToString());
	}

	void MetaMorpheusTask::FinishedSingleTask(const std::wstring &displayName)
	{
		SingleTaskEventArgs tempVar(displayName);
		FinishedSingleTaskHandler +== nullptr ? nullptr : FinishedSingleTaskHandler::Invoke(this, &tempVar);
	}

	void MetaMorpheusTask::StartingSingleTask(const std::wstring &displayName)
	{
		SingleTaskEventArgs tempVar(displayName);
		StartingSingleTaskHander +== nullptr ? nullptr : StartingSingleTaskHander::Invoke(this, &tempVar);
	}

	std::vector<std::type_info> MetaMorpheusTask::GetSubclassesAndItself(std::type_info type)
	{
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
		yield return type;
	}

	bool MetaMorpheusTask::SameSettings(const std::wstring &pathToOldParamsFile, IndexingEngine *indexEngine)
	{
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamReader reader = new StreamReader(pathToOldParamsFile))
		{
			StreamReader reader = StreamReader(pathToOldParamsFile);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			if (reader.ReadToEnd() == indexEngine->ToString())
			{
				return true;
			}
		}
		return false;
	}

	void MetaMorpheusTask::WritePeptideIndex(std::vector<PeptideWithSetModifications*> &peptideIndex, const std::wstring &peptideIndexFile)
	{
		auto messageTypes = GetSubclassesAndItself(std::vector<PeptideWithSetModifications*>::typeid);
		auto ser = new NetSerializer::Serializer(messageTypes);

//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.Create(peptideIndexFile))
		{
			auto file = File::Create(peptideIndexFile);
			ser->Serialize(file, peptideIndex);
		}

		delete ser;
	}

	void MetaMorpheusTask::WriteFragmentIndexNetSerializer(std::vector<std::vector<int>&> &fragmentIndex, const std::wstring &fragmentIndexFile)
	{
		auto messageTypes = GetSubclassesAndItself(std::vector<std::vector<int>>::typeid);
		auto ser = new NetSerializer::Serializer(messageTypes);

//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.Create(fragmentIndexFile))
		{
			auto file = File::Create(fragmentIndexFile);
			ser->Serialize(file, fragmentIndex);
		}

		delete ser;
	}

	std::wstring MetaMorpheusTask::GetExistingFolderWithIndices(IndexingEngine *indexEngine, std::vector<DbForTask*> &dbFilenameList)
	{
		for (auto database : dbFilenameList)
		{
			std::wstring baseDir = FileSystem::getDirectoryName(database->getFilePath());
			DirectoryInfo *indexDirectory = new DirectoryInfo(FileSystem::combine(baseDir, IndexFolderName));

			if (!FileSystem::directoryExists(indexDirectory->FullName))
			{
				delete indexDirectory;
				return L"";
			}

			// all directories in the same directory as the protein database
			std::vector<DirectoryInfo*> directories = indexDirectory->GetDirectories();

			// look in each subdirectory to find indexes folder
			for (auto possibleFolder : directories)
			{
				std::wstring result = CheckFiles(indexEngine, possibleFolder);

				if (result != L"")
				{
					delete indexDirectory;
					return result;
				}
			}

			delete indexDirectory;
		}

		return L"";
	}

	std::wstring MetaMorpheusTask::CheckFiles(IndexingEngine *indexEngine, DirectoryInfo *folder)
	{
		if (FileSystem::fileExists(FileSystem::combine(folder->FullName, L"indexEngine.params")) && FileSystem::fileExists(FileSystem::combine(folder->FullName, L"peptideIndex.ind")) && FileSystem::fileExists(FileSystem::combine(folder->FullName, L"fragmentIndex.ind")) && (FileSystem::fileExists(FileSystem::combine(folder->FullName, L"precursorIndex.ind")) || !indexEngine->GeneratePrecursorIndex) && SameSettings(FileSystem::combine(folder->FullName, L"indexEngine.params"), indexEngine))
		{
			return folder->FullName;
		}
		return L"";
	}

	void MetaMorpheusTask::WriteIndexEngineParams(IndexingEngine *indexEngine, const std::wstring &fileName)
	{
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(fileName))
		{
			StreamWriter output = StreamWriter(fileName);
			output.Write(indexEngine);
		}
	}

	std::wstring MetaMorpheusTask::GenerateOutputFolderForIndices(std::vector<DbForTask*> &dbFilenameList)
	{
		auto pathToIndexes = FileSystem::combine(FileSystem::getDirectoryName(dbFilenameList.front().FilePath), IndexFolderName);
		if (!FileSystem::fileExists(pathToIndexes))
		{
			FileSystem::createDirectory(pathToIndexes);
		}
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		auto folder = FileSystem::combine(pathToIndexes, DateTime::Now.ToString(L"yyyy-MM-dd-HH-mm-ss", CultureInfo::InvariantCulture));
		FileSystem::createDirectory(folder);
		return folder;
	}

	void MetaMorpheusTask::GenerateIndexes(IndexingEngine *indexEngine, std::vector<DbForTask*> &dbFilenameList, std::vector<PeptideWithSetModifications*> &peptideIndex, std::vector<std::vector<int>> &fragmentIndex, std::vector<std::vector<int>> &precursorIndex, std::vector<Protein*> &allKnownProteins, std::vector<Modification*> &allKnownModifications, const std::wstring &taskId)
	{
		std::wstring pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);
		if (pathToFolderWithIndices == L"")
		{
			auto output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
			Status(L"Writing params...", std::vector<std::wstring> {taskId});
			auto paramsFile = FileSystem::combine(output_folderForIndices, L"indexEngine.params");
			WriteIndexEngineParams(indexEngine, paramsFile);
			FinishedWritingFile(paramsFile, std::vector<std::wstring> {taskId});

			Status(L"Running Index Engine...", std::vector<std::wstring> {taskId});
			auto indexResults = static_cast<IndexingResults*>(indexEngine->Run());
			peptideIndex = indexResults->getPeptideIndex();
			fragmentIndex = indexResults->getFragmentIndex();
			precursorIndex = indexResults->getPrecursorIndex();

			Status(L"Writing peptide index...", std::vector<std::wstring> {taskId});
			auto peptideIndexFile = FileSystem::combine(output_folderForIndices, L"peptideIndex.ind");
			WritePeptideIndex(peptideIndex, peptideIndexFile);
			FinishedWritingFile(peptideIndexFile, std::vector<std::wstring> {taskId});

			Status(L"Writing fragment index...", std::vector<std::wstring> {taskId});
			auto fragmentIndexFile = FileSystem::combine(output_folderForIndices, L"fragmentIndex.ind");
			WriteFragmentIndexNetSerializer(fragmentIndex, fragmentIndexFile);
			FinishedWritingFile(fragmentIndexFile, std::vector<std::wstring> {taskId});

			if (indexEngine->GeneratePrecursorIndex)
			{
				Status(L"Writing precursor index...", std::vector<std::wstring> {taskId});
				auto precursorIndexFile = FileSystem::combine(output_folderForIndices, L"precursorIndex.ind");
				WriteFragmentIndexNetSerializer(precursorIndex, precursorIndexFile);
				FinishedWritingFile(precursorIndexFile, std::vector<std::wstring> {taskId});
			}
		}
		else
		{
			Status(L"Reading peptide index...", std::vector<std::wstring> {taskId});
			auto messageTypes = GetSubclassesAndItself(std::vector<PeptideWithSetModifications*>::typeid);
			auto ser = new NetSerializer::Serializer(messageTypes);
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "peptideIndex.ind")))
			{
				auto file = File::OpenRead(FileSystem::combine(pathToFolderWithIndices, L"peptideIndex.ind"));
				peptideIndex = static_cast<std::vector<PeptideWithSetModifications*>>(ser->Deserialize(file));
			}

			// populate dictionaries of known proteins for deserialization
			std::unordered_map<std::wstring, Protein*> proteinDictionary;

			for (auto protein : allKnownProteins)
			{
				if (proteinDictionary.find(protein->Accession) == proteinDictionary.end())
				{
					proteinDictionary.emplace(protein->Accession, protein);
				}
				else if (proteinDictionary[protein->Accession]->BaseSequence != protein->BaseSequence)
				{
					delete ser;
					throw MetaMorpheusException(StringHelper::formatSimple(L"The protein database contained multiple proteins with accession {0} ! This is not allowed for index-based searches (modern, non-specific, crosslink searches)", protein->Accession));
				}
			}

			// get non-serialized information for the peptides (proteins, mod info)
			for (auto peptide : peptideIndex)
			{
				peptide->SetNonSerializedPeptideInfo(GlobalVariables::getAllModsKnownDictionary(), proteinDictionary);
			}

			Status(L"Reading fragment index...", std::vector<std::wstring> {taskId});
			messageTypes = GetSubclassesAndItself(std::vector<std::vector<int>>::typeid);
			ser = new NetSerializer::Serializer(messageTypes);
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "fragmentIndex.ind")))
			{
				auto file = File::OpenRead(FileSystem::combine(pathToFolderWithIndices, L"fragmentIndex.ind"));
				fragmentIndex = static_cast<std::vector<std::vector<int>>>(ser->Deserialize(file));
			}

			if (indexEngine->GeneratePrecursorIndex)
			{
				Status(L"Reading precursor index...", std::vector<std::wstring> {taskId});
				messageTypes = GetSubclassesAndItself(std::vector<std::vector<int>>::typeid);
				ser = new NetSerializer::Serializer(messageTypes);
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "precursorIndex.ind")))
				{
					auto file = File::OpenRead(FileSystem::combine(pathToFolderWithIndices, L"precursorIndex.ind"));
					precursorIndex = static_cast<std::vector<std::vector<int>>>(ser->Deserialize(file));
				}
			}

			delete ser;
		}
	}
}
