#include "GPTMDTask.h"
#include "GPTMDParameters.h"
#include "../../EngineLayer/CommonParameters.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"
#include "../MyTaskResults.h"
#include "../../EngineLayer/GlobalVariables.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../../EngineLayer/PrecursorSearchModes/DotMassDiffAcceptor.h"
#include "../../EngineLayer/PrecursorSearchModes/MassDiffAcceptor.h"
#include "../MyFileManager.h"
#include "../../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../../EngineLayer/EventArgs/ProgressEventArgs.h"

using namespace EngineLayer;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::Gptmd;
#if defined(NETFRAMEWORK)
using namespace IO::Thermo;
#endif
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace UsefulProteomicsDatabases;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{

	GptmdTask::GptmdTask() : MetaMorpheusTask(MyTask::Gptmd)
	{
		EngineLayer::CommonParameters tempVar();
		setCommonParameters(&tempVar);
		EngineLayer::GptmdParameters tempVar2();
		setGptmdParameters(&tempVar2);
	}

	EngineLayer::GptmdParameters *GptmdTask::getGptmdParameters() const
	{
		return privateGptmdParameters;
	}

	void GptmdTask::setGptmdParameters(EngineLayer::GptmdParameters *value)
	{
		privateGptmdParameters = value;
	}

	MyTaskResults *GptmdTask::RunSpecific(const std::string &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::string> &currentRawFileList, const std::string &taskId, std::vector<FileSpecificParameters*> &fileSettingsList)
	{
		std::vector<Modification> variableModifications;
		std::vector<Modification> fixedModifications;
		std::vector<string> localizeableModificationTypes;
		LoadModifications(taskId, variableModifications, fixedModifications, localizeableModificationTypes);

		// TODO: print error messages loading GPTMD mods
		std::vector<Modification*> gptmdModifications = GlobalVariables::getAllModsKnown().OfType<Modification*>().Where([&] (std::any b)
		{
			getGptmdParameters()->ListOfModsGptmd->Contains((b::ModificationType, b::IdWithMotif));
		}).ToList();
		std::vector<std::tuple<double, double>> combos = LoadCombos(gptmdModifications).ToList();

		// load proteins
		std::vector<Protein*> proteinList = LoadProteins(taskId, dbFilenameList, true, DecoyType::Reverse, localizeableModificationTypes, getCommonParameters());

		std::vector<PeptideSpectralMatch*> allPsms;

		auto numRawFiles = currentRawFileList.size();

		// write prose settings
		ProseCreatedWhileRunning->append("The following G-PTM-D settings were used: ");
		ProseCreatedWhileRunning->append("protease = " + getCommonParameters()->getDigestionParams()->Protease + "; ");
		ProseCreatedWhileRunning->append("maximum missed cleavages = " + getCommonParameters()->getDigestionParams()->MaxMissedCleavages + "; ");
		ProseCreatedWhileRunning->append("minimum peptide length = " + getCommonParameters()->getDigestionParams()->MinPeptideLength + "; ");
		ProseCreatedWhileRunning->append(getCommonParameters()->getDigestionParams()->MaxPeptideLength == std::numeric_limits<int>::max() ? "maximum peptide length = unspecified; " : "maximum peptide length = " + getCommonParameters()->getDigestionParams()->MaxPeptideLength + "; ");
		ProseCreatedWhileRunning->append("initiator methionine behavior = " + getCommonParameters()->getDigestionParams()->InitiatorMethionineBehavior + "; ");
		ProseCreatedWhileRunning->append("max modification isoforms = " + getCommonParameters()->getDigestionParams()->MaxModificationIsoforms + "; ");
		ProseCreatedWhileRunning->append("fixed modifications = " + std::string::Join(", ", fixedModifications->Select([&] (std::any m)
		{
			m::IdWithMotif;
		})) + "; ");
		ProseCreatedWhileRunning->append("variable modifications = " + std::string::Join(", ", variableModifications->Select([&] (std::any m)
		{
			m::IdWithMotif;
		})) + "; ");
		ProseCreatedWhileRunning->append("G-PTM-D modifications count = " + std::to_string(gptmdModifications.size()) + "; ");

		// temporary search type for writing prose
		// the actual search type is technically file-specific but we don't allow file-specific notches, so it's safe to do this
		MassDiffAcceptor *tempSearchMode = new DotMassDiffAcceptor("", GetAcceptableMassShifts(fixedModifications, variableModifications, gptmdModifications, combos), getCommonParameters()->getPrecursorMassTolerance());
		ProseCreatedWhileRunning->append("precursor mass tolerance(s) = {" + tempSearchMode->ToProseString() + "}; ");

		ProseCreatedWhileRunning->append("product mass tolerance = " + getCommonParameters()->getProductMassTolerance() + ". ");
		ProseCreatedWhileRunning->append("The combined search database contained " + proteinList.size()([&] (std::any p)
		{
			!p::IsDecoy;
		}) + " non-decoy protein entries including " + proteinList.Where([&] (std::any p)
		{
			p::IsContaminant;
		})->Count() + " contaminant sequences. ");

		// start the G-PTM-D task
		Status("Running G-PTM-D...", std::vector<std::string> {taskId});
		MyTaskResults = new MyTaskResults(this);
		MyTaskResults->NewDatabases = std::vector<DbForTask*>();
		auto fileSpecificCommonParams = fileSettingsList.Select([&] (std::any b)
		{
			SetAllFileSpecificCommonParams(getCommonParameters(), b);
		});
		std::unordered_set<DigestionParams*> ListOfDigestionParams = std::unordered_set<DigestionParams*>(fileSpecificCommonParams->Select([&] (std::any p)
		{
			p::DigestionParams;
		}));

		MyFileManager *myFileManager = new MyFileManager(true);

		std::any lock1 = std::any();
		std::any lock2 = std::any();

		for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.size(); spectraFileIndex++)
		{
			// Stop if canceled
			if (GlobalVariables::getStopLoops())
			{
				break;
			}

			auto origDataFile = currentRawFileList[spectraFileIndex];

			// mark the file as in-progress
			StartingDataFile(origDataFile, std::vector<std::string> {taskId, "Individual Spectra Files", origDataFile});

			EngineLayer::CommonParameters *combinedParams = SetAllFileSpecificCommonParams(getCommonParameters(), fileSettingsList[spectraFileIndex]);
			MassDiffAcceptor *searchMode = new DotMassDiffAcceptor("", GetAcceptableMassShifts(fixedModifications, variableModifications, gptmdModifications, combos), combinedParams->getPrecursorMassTolerance());

			NewCollection(FileSystem::getFileName(origDataFile), std::vector<std::string> {taskId, "Individual Spectra Files", origDataFile});

			Status("Loading spectra file...", std::vector<std::string> {taskId, "Individual Spectra Files", origDataFile});
			MsDataFile *myMsDataFile = myFileManager->LoadFile(origDataFile, std::make_optional(combinedParams->getTopNpeaks()), std::make_optional(combinedParams->getMinRatio()), combinedParams->getTrimMs1Peaks(), combinedParams->getTrimMsMsPeaks(), combinedParams);
			Status("Getting ms2 scans...", std::vector<std::string> {taskId, "Individual Spectra Files", origDataFile});
			std::vector<Ms2ScanWithSpecificMass*> arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy([&] (std::any b)
			{
				b::PrecursorMass;
			})->ToArray();
			myFileManager->DoneWithFile(origDataFile);
			std::vector<PeptideSpectralMatch*> allPsmsArray(arrayOfMs2ScansSortedByMass.size());
			ClassicSearchEngine tempVar(allPsmsArray, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, searchMode, combinedParams, new std::vector<std::string> {taskId, "Individual Spectra Files", origDataFile});
			(&tempVar)->Run();
			allPsms.AddRange(allPsmsArray.Where([&] (std::any p)
			{
//C# TO C++ CONVERTER TODO TASK: A 'delete searchMode' statement was not added since searchMode was passed to a method or constructor. Handle memory management manually.
			delete myFileManager;
			delete tempSearchMode;
				return p != nullptr;
			}));
			FinishedDataFile(origDataFile, std::vector<std::string> {taskId, "Individual Spectra Files", origDataFile});
			ProgressEventArgs tempVar2(100, "Done!", new std::vector<std::string> {taskId, "Individual Spectra Files", origDataFile});
			ReportProgress(&tempVar2);

//C# TO C++ CONVERTER TODO TASK: A 'delete searchMode' statement was not added since searchMode was passed to a method or constructor. Handle memory management manually.
		}
		ProgressEventArgs tempVar3(100, "Done!", new std::vector<std::string> {taskId, "Individual Spectra Files"});
		ReportProgress(&tempVar3);

		allPsms = allPsms.OrderByDescending([&] (std::any b)
		{
			b::Score;
		}).ThenBy([&] (std::any b)
		{
			b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
		}).GroupBy([&] (std::any b)
		{
		delete myFileManager;
		delete tempSearchMode;
			return std::tuple < std::string;
		}, int, std::optional<double>>(b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass))->Select([&] (std::any b)
		{
			b::First();
		}).ToList();

		FdrAnalysisEngine tempVar4(allPsms, tempSearchMode->NumNotches, getCommonParameters(), new std::vector<std::string> {taskId});
		(&tempVar4)->Run();

		auto writtenFile = FileSystem::combine(OutputFolder, "GPTMD_Candidates.psmtsv");
		WritePsmsToTsv(allPsms, writtenFile, std::unordered_map<std::string, int>());
		FinishedWritingFile(writtenFile, std::vector<std::string> {taskId});

		// get file-specific precursor mass tolerances for the GPTMD engine
		auto filePathToPrecursorMassTolerance = std::unordered_map<std::string, Tolerance*>();
		for (int i = 0; i < currentRawFileList.size(); i++)
		{
			std::string filePath = currentRawFileList[i];
			Tolerance *fileTolerance = getCommonParameters()->getPrecursorMassTolerance();
			if (fileSettingsList[i] != nullptr && fileSettingsList[i]->getPrecursorMassTolerance() != nullptr)
			{
				fileTolerance = fileSettingsList[i]->getPrecursorMassTolerance();
			}
			filePathToPrecursorMassTolerance.emplace(filePath, fileTolerance);
		}

		// run GPTMD engine
		GptmdEngine tempVar5(allPsms, gptmdModifications, combos, filePathToPrecursorMassTolerance, getCommonParameters(), new std::vector<std::string> {taskId});
		auto gptmdResults = static_cast<GptmdResults*>((&tempVar5)->Run());

		// Stop if canceled
		if (GlobalVariables::getStopLoops())
		{
			delete myFileManager;
			delete tempSearchMode;
			return MyTaskResults;
		}

		// write GPTMD databases
		if (dbFilenameList.Any([&] (std::any b)
		{
			!b::IsContaminant;
		}))
		{
			std::vector<std::string> databaseNames;
			for (auto nonContaminantDb : dbFilenameList.Where([&] (std::any p)
			{
				!p::IsContaminant;
			}))
			{
				auto dbName = Path::GetFileNameWithoutExtension(nonContaminantDb::FilePath);
//C# TO C++ CONVERTER TODO TASK: There is no direct native C++ equivalent to this .NET String method:
				auto theExtension = Path::GetExtension(nonContaminantDb::FilePath).ToLowerInvariant();
				bool compressed = StringHelper::endsWith(theExtension, "gz");
				databaseNames.push_back(compressed ? Path::GetFileNameWithoutExtension(dbName) : dbName);
			}
			std::string outputXMLdbFullName = FileSystem::combine(OutputFolder, std::string::Join("-", databaseNames) + "GPTMD.xml");

			auto newModsActuallyWritten = ProteinDbWriter::WriteXmlDatabase(gptmdResults->getMods(), proteinList.Where([&] (std::any b)
			{
			delete myFileManager;
			delete tempSearchMode;
				return !b::IsDecoy && !b::IsContaminant;
			}).ToList(), outputXMLdbFullName);

			FinishedWritingFile(outputXMLdbFullName, std::vector<std::string> {taskId});

			DbForTask tempVar6(outputXMLdbFullName, false);
			MyTaskResults->NewDatabases.push_back(&tempVar6);
			MyTaskResults->AddNiceText("Modifications added: " + newModsActuallyWritten->Select([&] (std::any b)
			{
				b->Value;
			}).Sum());
			MyTaskResults->AddNiceText("Mods types and counts:");
			MyTaskResults->AddNiceText(std::string::Join("\r\n", newModsActuallyWritten->OrderByDescending([&] (std::any b)
			{
				b->Value;
			})->Select([&] (std::any b)
			{
			delete myFileManager;
			delete tempSearchMode;
				return "\t" + b::Key + "\t" + b->Value;
			})));
		}
		if (dbFilenameList.Any([&] (std::any b)
		{
			b::IsContaminant;
		}))
		{
			// do NOT use this code (Path.GetFilenameWithoutExtension) because GPTMD on .xml.gz will result in .xml.xml file type being written
			//string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "GPTMD.xml");
			std::vector<std::string> databaseNames;
			for (auto contaminantDb : dbFilenameList.Where([&] (std::any p)
			{
				p::IsContaminant;
			}))
			{
				auto dbName = FileSystem::getFileName(contaminantDb::FilePath);
				int indexOfFirstDot = (int)dbName.find(".");
				databaseNames.push_back(dbName.substr(0, indexOfFirstDot));
			}
			std::string outputXMLdbFullNameContaminants = FileSystem::combine(OutputFolder, std::string::Join("-", databaseNames) + "GPTMD.xml");

			auto newModsActuallyWritten = ProteinDbWriter::WriteXmlDatabase(gptmdResults->getMods(), proteinList.Where([&] (std::any b)
			{
			delete myFileManager;
			delete tempSearchMode;
				return !b::IsDecoy && b::IsContaminant;
			}).ToList(), outputXMLdbFullNameContaminants);

			FinishedWritingFile(outputXMLdbFullNameContaminants, std::vector<std::string> {taskId});

			DbForTask tempVar7(outputXMLdbFullNameContaminants, true);
			MyTaskResults->NewDatabases.push_back(&tempVar7);
			MyTaskResults->AddNiceText("Contaminant modifications added: " + newModsActuallyWritten->Select([&] (std::any b)
			{
				b->Value;
			}).Sum());
			MyTaskResults->AddNiceText("Mods types and counts:");
			MyTaskResults->AddNiceText(std::string::Join("\r\n", newModsActuallyWritten->OrderByDescending([&] (std::any b)
			{
				b->Value;
			})->Select([&] (std::any b)
			{
			delete myFileManager;
			delete tempSearchMode;
				return "\t" + b::Key + "\t" + b->Value;
			})));
		}

		delete myFileManager;
		delete tempSearchMode;
		return MyTaskResults;
	}

	std::vector<std::tuple<double, double>> GptmdTask::LoadCombos(std::vector<Modification*> &modificationsThatCanBeCombined)
	{
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamReader r = new StreamReader(Path.Combine(GlobalVariables.DataDir, "Data", "TangibleTempVerbatimOpenTagcombos.txtTangibleTempVerbatimCloseTag")))
		{
			StreamReader r = StreamReader(FileSystem::combine(GlobalVariables::getDataDir(), "Data", "combos.txt"));
			while (r.Peek() >= 0)
			{
				auto line = StringHelper::split(r.ReadLine(), L' ');
				auto mass1 = std::stod(line[0]);
				auto mass2 = std::stod(line[1]);
				if (modificationsThatCanBeCombined.Where([&] (std::any b)
				{
					return b->ValidModification == true;
				}).Any([&] (std::any b)
				{
					return std::abs(static_cast<double>(b::MonoisotopicMass) - mass1) < tolForComboLoading;
				}) && modificationsThatCanBeCombined.Where([&] (std::any b)
				{
					return b->ValidModification == true;
				}).Any([&] (std::any b)
				{
					return std::abs(static_cast<double>(b::MonoisotopicMass) - mass2) < tolForComboLoading;
				}))yield return std::tuple<double, double>(mass1, mass2);
				{
				}
			}
		}
	}

	std::vector<double> GptmdTask::GetAcceptableMassShifts(std::vector<Modification*> &fixedMods, std::vector<Modification*> &variableMods, std::vector<Modification*> &gptmdMods, std::vector<std::tuple<double, double>> &combos)
	{
		std::vector<double> gptmdNotches = gptmdMods.Where([&] (std::any b)
		{
			return b->ValidModification == true;
		})->Select([&] (std::any b)
		{
			(double)b.MonoisotopicMass;
		});
		std::vector<double> gptmdMinusOtherModsNotches = GetObservedMasses(variableMods.Concat(fixedMods), gptmdMods);
		std::vector<double> multipleGptmdNotches = combos.Select([&] (std::any b)
		{
			return b::Item1 + b::Item2;
		});
		std::vector<double> zeroNotch = {0};

		std::vector<double> allNotches = gptmdNotches.Concat(gptmdMinusOtherModsNotches)->Concat(multipleGptmdNotches)->Concat(zeroNotch);
		return allNotches.GroupBy([&] (std::any b)
		{
			std::round(b * std::pow(10, 5)) / std::pow(10, 5);
		})->Select([&] (std::any b)
		{
			b::FirstOrDefault();
		}).OrderBy([&] (std::any b)
		{
			return b;
		});
	}

	std::vector<double> GptmdTask::GetObservedMasses(std::vector<Modification*> &enumerable, std::vector<Modification*> &gptmdModifications)
	{
		for (auto modOnPeptide : enumerable.Where([&] (std::any b)
		{
			return b->ValidModification == true;
		}))
		{
			for (auto modToLocalize : gptmdModifications.Where([&] (std::any b)
			{
				return b->ValidModification == true;
			}))
			{
				if (modOnPeptide::Target->Equals(modToLocalize::Target))
				{
//C# TO C++ CONVERTER TODO TASK: C++ does not have an equivalent to the C# 'yield' keyword:
					yield return static_cast<double>(modToLocalize::MonoisotopicMass) - static_cast<double>(modOnPeptide::MonoisotopicMass);
				}
			}
		}
	}
}
