#include "TaskLayer.XLSearchTask.h"
#include "XLSearchParameters.h"
#include "../../EngineLayer/CommonParameters.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"
#include "../MyTaskResults.h"
#include "../MyFileManager.h"
#include "../../EngineLayer/Ms2ScanWithSpecificMass.h"
#include "../../EngineLayer/GlobalVariables.h"
#include "../../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"

using namespace EngineLayer;
using namespace EngineLayer::CrosslinkSearch;
using namespace EngineLayer::Indexing;
using namespace MassSpectrometry;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace MzLibUtil;
using namespace EngineLayer::FdrAnalysis;
using namespace Proteomics::Fragmentation;

namespace TaskLayer
{

	XLSearchTask::XLSearchTask() : MetaMorpheusTask(MyTask::XLSearch)
	{
		EngineLayer::CommonParameters tempVar(, , = true, = true, = 3, = 12, = true, = false, = 1, 3, , , , , , , , new PpmTolerance(10));
		setCommonParameters(&tempVar);

		TaskLayer::XlSearchParameters tempVar2();
		setXlSearchParameters(&tempVar2);
	}

	TaskLayer::XlSearchParameters *XLSearchTask::getXlSearchParameters() const
	{
		return privateXlSearchParameters;
	}

	void XLSearchTask::setXlSearchParameters(TaskLayer::XlSearchParameters *value)
	{
		privateXlSearchParameters = value;
	}

	MyTaskResults *XLSearchTask::RunSpecific(const std::string &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::string> &currentRawFileList, const std::string &taskId, std::vector<FileSpecificParameters*> &fileSettingsList)
	{
		MyTaskResults = new MyTaskResults(this);
		std::vector<CrosslinkSpectralMatch*> allPsms;

		std::vector<Modification> variableModifications;
		std::vector<Modification> fixedModifications;
		std::vector<string> localizeableModificationTypes;
		LoadModifications(taskId, variableModifications, fixedModifications, localizeableModificationTypes);

		// load proteins
		std::vector<Protein*> proteinList = LoadProteins(taskId, dbFilenameList, true, getXlSearchParameters()->getDecoyType(), localizeableModificationTypes, getCommonParameters());

		auto crosslinker = new Crosslinker();
		crosslinker = crosslinker->SelectCrosslinker(getXlSearchParameters()->getCrosslinkerType());
		if (getXlSearchParameters()->getCrosslinkerType() == CrosslinkerType::UserDefined)
		{
			crosslinker = GenerateUserDefinedCrosslinker(getXlSearchParameters());
		}

		MyFileManager *myFileManager = new MyFileManager(true);

		auto fileSpecificCommonParams = fileSettingsList.Select([&] (std::any b)
		{
			SetAllFileSpecificCommonParams(getCommonParameters(), b);
		});
		std::unordered_set<DigestionParams*> ListOfDigestionParams = std::unordered_set<DigestionParams*>(fileSpecificCommonParams->Select([&] (std::any p)
		{
			p::DigestionParams;
		}));

		int completedFiles = 0;
		std::any indexLock = std::any();
		std::any psmLock = std::any();

		Status("Searching files...", taskId);

		ProseCreatedWhileRunning->append("The following crosslink discovery were used: ");
		ProseCreatedWhileRunning->append("crosslinker name = " + crosslinker->getCrosslinkerName() + "; ");
		ProseCreatedWhileRunning->append("crosslinker type = " + StringHelper::toString(crosslinker->getCleavable()) + "; ");
		ProseCreatedWhileRunning->append("crosslinker mass = " + std::to_string(crosslinker->getTotalMass()) + "; ");
		ProseCreatedWhileRunning->append("crosslinker modification site(s) = " + crosslinker->getCrosslinkerModSites() + "; ");

		ProseCreatedWhileRunning->append("protease = " + getCommonParameters()->getDigestionParams()->Protease + "; ");
		ProseCreatedWhileRunning->append("maximum missed cleavages = " + getCommonParameters()->getDigestionParams()->MaxMissedCleavages + "; ");
		ProseCreatedWhileRunning->append("minimum peptide length = " + getCommonParameters()->getDigestionParams()->MinPeptideLength + "; ");
		ProseCreatedWhileRunning->append(getCommonParameters()->getDigestionParams()->MaxPeptideLength == std::numeric_limits<int>::max() ? "maximum peptide length = unspecified; " : "maximum peptide length = " + getCommonParameters()->getDigestionParams()->MaxPeptideLength + "; ");
		ProseCreatedWhileRunning->append("initiator methionine behavior = " + getCommonParameters()->getDigestionParams()->InitiatorMethionineBehavior + "; ");
		ProseCreatedWhileRunning->append("max modification isoforms = " + getCommonParameters()->getDigestionParams()->MaxModificationIsoforms + "; ");

		ProseCreatedWhileRunning->append("fixed modifications = " + std::string::Join(", ", fixedModifications->Select([&] (std::any m)
		{
			m::IdWithMotif;
		}) + "; "));
		ProseCreatedWhileRunning->append("variable modifications = " + std::string::Join(", ", variableModifications->Select([&] (std::any m)
		{
			m::IdWithMotif;
		})) + "; ");

		ProseCreatedWhileRunning->append("parent mass tolerance(s) = " + getCommonParameters()->getPrecursorMassTolerance() + "; ");
		ProseCreatedWhileRunning->append("product mass tolerance = " + getCommonParameters()->getProductMassTolerance() + "; ");
		ProseCreatedWhileRunning->append("The combined search database contained " + std::to_string(proteinList.size()) + " total entries including " + proteinList.Where([&] (std::any p)
		{
			p::IsContaminant;
		})->Count() + " contaminant sequences. ");

		for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.size(); spectraFileIndex++)
		{
			auto origDataFile = currentRawFileList[spectraFileIndex];
			EngineLayer::CommonParameters *combinedParams = SetAllFileSpecificCommonParams(getCommonParameters(), fileSettingsList[spectraFileIndex]);

			auto thisId = std::vector<std::string> {taskId, "Individual Spectra Files", origDataFile};
			NewCollection(FileSystem::getFileName(origDataFile), thisId);

			Status("Loading spectra file...", thisId);
			MsDataFile *myMsDataFile = myFileManager->LoadFile(origDataFile, std::make_optional(combinedParams->getTopNpeaks()), std::make_optional(combinedParams->getMinRatio()), combinedParams->getTrimMs1Peaks(), combinedParams->getTrimMsMsPeaks(), combinedParams);

			Status("Getting ms2 scans...", thisId);
			std::vector<Ms2ScanWithSpecificMass*> arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy([&] (std::any b)
			{
				b::PrecursorMass;
			})->ToArray();

			std::vector<CrosslinkSpectralMatch*> newPsms(arrayOfMs2ScansSortedByMass.size());
			for (int currentPartition = 0; currentPartition < getCommonParameters()->getTotalPartitions(); currentPartition++)
			{
				std::vector<PeptideWithSetModifications*> peptideIndex;
				std::vector<Protein*> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.size()() / combinedParams->getTotalPartitions(), ((currentPartition + 1) * proteinList.size()() / combinedParams->getTotalPartitions()) - (currentPartition * proteinList.size()() / combinedParams->getTotalPartitions()));

				Status("Getting fragment dictionary...", std::vector<std::string> {taskId});
				auto indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, currentPartition, UsefulProteomicsDatabases::DecoyType::Reverse, combinedParams, 30000.0, false, dbFilenameList.Select([&] (std::any p)
				{
					new FileInfo(p::FilePath);
				}).ToList(), std::vector<std::string> {taskId});
				std::vector<std::vector<int>> fragmentIndex;
				std::vector<std::vector<int>> precursorIndex;

				GenerateIndexes(indexEngine, dbFilenameList, peptideIndex, fragmentIndex, precursorIndex, proteinList, GlobalVariables::getAllModsKnown().ToList(), taskId);

				Status("Searching files...", taskId);
				CrosslinkSearchEngine tempVar(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition, combinedParams, crosslinker, getXlSearchParameters()->getRestrictToTopNHits(), getXlSearchParameters()->getCrosslinkSearchTopNum(), getXlSearchParameters()->getXlQuench_H2O(), getXlSearchParameters()->getXlQuench_NH2(), getXlSearchParameters()->getXlQuench_Tris(), thisId);
				(&tempVar)->Run();

				ProgressEventArgs tempVar2(100, "Done with search " + std::to_string(currentPartition + 1) + "/" + std::to_string(getCommonParameters()->getTotalPartitions()) + "!", thisId);
				ReportProgress(&tempVar2);

//C# TO C++ CONVERTER TODO TASK: A 'delete indexEngine' statement was not added since indexEngine was passed to a method or constructor. Handle memory management manually.
			}

			allPsms.AddRange(newPsms.Where([&] (std::any p)
			{
			delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
				return p != nullptr;
			}));

			completedFiles++;
			ProgressEventArgs tempVar3(completedFiles / currentRawFileList.size(), "Searching...", new std::vector<std::string> {taskId, "Individual Spectra Files"});
			ReportProgress(&tempVar3);
		}

		ProgressEventArgs tempVar4(100, "Done with all searches!", new std::vector<std::string> {taskId, "Individual Spectra Files"});
		ReportProgress(&tempVar4);

		allPsms = allPsms.OrderByDescending([&] (std::any p)
		{
			p::XLTotalScore;
		}).ToList();

		auto allPsmsXL = allPsms.Where([&] (std::any p)
		{
		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
			return p->CrossType == PsmCrossType::Cross;
		}).ToList();

		// inter-crosslinks; different proteins are linked
		auto interCsms = allPsmsXL.Where([&] (std::any p)
		{
			!p::ProteinAccession->Equals(p::BetaPeptide::ProteinAccession);
		}).ToList();
		for (auto item : interCsms)
		{
			item.CrossType = PsmCrossType::Inter;
		}

		// intra-crosslinks; crosslinks within a protein
		auto intraCsms = allPsmsXL.Where([&] (std::any p)
		{
			p::ProteinAccession->Equals(p::BetaPeptide::ProteinAccession);
		}).ToList();
		for (auto item : intraCsms)
		{
			item.CrossType = PsmCrossType::Intra;
		}

		// calculate FDR
		DoCrosslinkFdrAnalysis(interCsms);
		DoCrosslinkFdrAnalysis(intraCsms);
		SingleFDRAnalysis(allPsms, std::vector<std::string> {taskId});

		// calculate protein crosslink residue numbers
		for (auto csm : allPsmsXL)
		{
			// alpha peptide crosslink residue in the protein
			csm.XlProteinPos = csm.OneBasedStartResidueInProtein->Value + csm.LinkPositions[0] - 1;

			// beta crosslink residue in protein
			csm.BetaPeptide->XlProteinPos = csm.BetaPeptide::OneBasedStartResidueInProtein->Value + csm.BetaPeptide::LinkPositions[0] - 1;
		}

		// write interlink CSMs
		if (interCsms.Any())
		{
			std::string file = FileSystem::combine(OutputFolder, "XL_Interlinks.tsv");
			WritePsmCrossToTsv(interCsms, file, 2);
			FinishedWritingFile(file, std::vector<std::string> {taskId});
		}
		MyTaskResults->AddNiceText("Target inter-crosslinks within 1% FDR: " + interCsms.size()([&] (std::any p)
		{
		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
			return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy && !p::BetaPeptide::IsDecoy;
		}));

		if (getXlSearchParameters()->getWriteOutputForPercolator())
		{
			auto interPsmsXLPercolator = interCsms.Where([&] (std::any p)
			{
			delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
				return p::Score >= 2 && p::BetaPeptide::Score >= 2;
			}).OrderBy([&] (std::any p)
			{
				p::ScanNumber;
			}).ToList();
			WriteCrosslinkToTxtForPercolator(interPsmsXLPercolator, OutputFolder, "XL_Interlinks_Percolator", crosslinker, std::vector<std::string> {taskId});
		}

		// write intralink CSMs
		if (intraCsms.Any())
		{
			std::string file = FileSystem::combine(OutputFolder, "XL_Intralinks.tsv");
			WritePsmCrossToTsv(intraCsms, file, 2);
			FinishedWritingFile(file, std::vector<std::string> {taskId});
		}
		MyTaskResults->AddNiceText("Target intra-crosslinks within 1% FDR: " + intraCsms.size()([&] (std::any p)
		{
		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
			return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy && !p::BetaPeptide::IsDecoy;
		}));

		if (getXlSearchParameters()->getWriteOutputForPercolator())
		{
			auto intraPsmsXLPercolator = intraCsms.Where([&] (std::any p)
			{
			delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
				return p::Score >= 2 && p::BetaPeptide::Score >= 2;
			}).OrderBy([&] (std::any p)
			{
				p::ScanNumber;
			}).ToList();
			WriteCrosslinkToTxtForPercolator(intraPsmsXLPercolator, OutputFolder, "XL_Intralinks_Percolator", crosslinker, std::vector<std::string> {taskId});
		}

		// write single peptides
		auto singlePsms = allPsms.Where([&] (std::any p)
		{
		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
			return p->CrossType == PsmCrossType::Single;
		}).ToList();
		if (singlePsms.Any())
		{
			std::string writtenFileSingle = FileSystem::combine(OutputFolder, std::string("SinglePeptides") + ".tsv");
			WritePsmCrossToTsv(singlePsms, writtenFileSingle, 1);
			FinishedWritingFile(writtenFileSingle, std::vector<std::string> {taskId});
		}
		MyTaskResults->AddNiceText("Target single peptides within 1% FDR: " + singlePsms.size()([&] (std::any p)
		{
		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
			return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy;
		}));

		// write loops
		auto loopPsms = allPsms.Where([&] (std::any p)
		{
		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
			return p->CrossType == PsmCrossType::Loop;
		}).ToList();
		if (loopPsms.Any())
		{
			std::string writtenFileLoop = FileSystem::combine(OutputFolder, std::string("Looplinks") + ".tsv");
			WritePsmCrossToTsv(loopPsms, writtenFileLoop, 1);
			FinishedWritingFile(writtenFileLoop, std::vector<std::string> {taskId});
		}
		MyTaskResults->AddNiceText("Target loop-linked peptides within 1% FDR: " + loopPsms.size()([&] (std::any p)
		{
		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
			return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy;
		}));

		// write deadends
		auto deadendPsms = allPsms.Where([&] (std::any p)
		{
		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
			return p->CrossType == PsmCrossType::DeadEnd || p->CrossType == PsmCrossType::DeadEndH2O || p->CrossType == PsmCrossType::DeadEndNH2 || p->CrossType == PsmCrossType::DeadEndTris;
		}).ToList();
		if (deadendPsms.Any())
		{
			std::string writtenFileDeadend = FileSystem::combine(OutputFolder, std::string("Deadends") + ".tsv");
			WritePsmCrossToTsv(deadendPsms, writtenFileDeadend, 1);
			FinishedWritingFile(writtenFileDeadend, std::vector<std::string> {taskId});
		}
		MyTaskResults->AddNiceText("Target deadend peptides within 1% FDR: " + deadendPsms.size()([&] (std::any p)
		{
		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
			return p::FdrInfo::QValue <= 0.01 && !p::IsDecoy;
		}));

		// write pepXML
		if (getXlSearchParameters()->getWritePepXml())
		{
			std::vector<CrosslinkSpectralMatch*> writeToXml;
			writeToXml.AddRange(intraCsms.Where([&] (std::any p)
			{
			delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
				return !p::IsDecoy && !p::BetaPeptide::IsDecoy && p::FdrInfo::QValue <= 0.05;
			}));
			writeToXml.AddRange(interCsms.Where([&] (std::any p)
			{
			delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
				return !p::IsDecoy && !p::BetaPeptide::IsDecoy && p::FdrInfo::QValue <= 0.05;
			}));
			writeToXml.AddRange(singlePsms.Where([&] (std::any p)
			{
			delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
				return !p::IsDecoy && p::FdrInfo::QValue <= 0.05;
			}));
			writeToXml.AddRange(loopPsms.Where([&] (std::any p)
			{
			delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
				return !p::IsDecoy && p::FdrInfo::QValue <= 0.05;
			}));
			writeToXml.AddRange(deadendPsms.Where([&] (std::any p)
			{
			delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
				return !p::IsDecoy && p::FdrInfo::QValue <= 0.05;
			}));
			writeToXml = writeToXml.OrderBy([&] (std::any p)
			{
				p::ScanNumber;
			}).ToList();

			for (auto fullFilePath : currentRawFileList)
			{
				std::string fileNameNoExtension = Path::GetFileNameWithoutExtension(fullFilePath);
				WritePepXML_xl(writeToXml.Where([&] (std::any p)
				{
				delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
					return p->FullFilePath == fullFilePath;
				}).ToList(), proteinList, dbFilenameList[0]->getFilePath(), variableModifications, fixedModifications, localizeableModificationTypes, OutputFolder, fileNameNoExtension, std::vector<std::string> {taskId});
			}
		}

		delete myFileManager;
//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was passed to a method or constructor. Handle memory management manually.
		return MyTaskResults;
	}

	void XLSearchTask::SingleFDRAnalysis(std::vector<CrosslinkSpectralMatch*> &items, std::vector<std::string> &taskIds)
	{
		// calculate single PSM FDR
		std::vector<PeptideSpectralMatch*> psms = items.Where([&] (std::any p)
		{
			return p->CrossType == PsmCrossType::Single;
		})->Select([&] (std::any p)
		{
			dynamic_cast<PeptideSpectralMatch*>(p);
		}).ToList();
		FdrAnalysisEngine tempVar(psms, 0, getCommonParameters(), taskIds);
		(&tempVar)->Run();

		// calculate loop PSM FDR
		psms = items.Where([&] (std::any p)
		{
			return p->CrossType == PsmCrossType::Loop;
		})->Select([&] (std::any p)
		{
			dynamic_cast<PeptideSpectralMatch*>(p);
		}).ToList();
		FdrAnalysisEngine tempVar2(psms, 0, getCommonParameters(), taskIds);
		(&tempVar2)->Run();

		// calculate deadend FDR
		psms = items.Where([&] (std::any p)
		{
			return p->CrossType == PsmCrossType::DeadEnd || p->CrossType == PsmCrossType::DeadEndH2O || p->CrossType == PsmCrossType::DeadEndNH2 || p->CrossType == PsmCrossType::DeadEndTris;
		})->Select([&] (std::any p)
		{
			dynamic_cast<PeptideSpectralMatch*>(p);
		}).ToList();
		FdrAnalysisEngine tempVar3(psms, 0, getCommonParameters(), taskIds);
		(&tempVar3)->Run();
	}

	void XLSearchTask::DoCrosslinkFdrAnalysis(std::vector<CrosslinkSpectralMatch*> &csms)
	{
		int cumulativeTarget = 0;
		int cumulativeDecoy = 0;

		for (int i = 0; i < csms.size(); i++)
		{
			auto csm = csms[i];
			if (csm->getIsDecoy() || csm->getBetaPeptide()->getIsDecoy())
			{
				cumulativeDecoy++;
			}
			else
			{
				cumulativeTarget++;
			}

			double qValue = std::min(1, static_cast<double>(cumulativeDecoy) / cumulativeTarget);
			csm->SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, 0, 0, 0, 0, 0, 0, false);
		}

		double qValueThreshold = 1.0;
		for (int i = csms.size() - 1; i >= 0; i--)
		{
			CrosslinkSpectralMatch *csm = csms[i];

			// threshold q-values
			if (csm->getFdrInfo()->getQValue() > qValueThreshold)
			{
				csm->getFdrInfo()->setQValue(qValueThreshold);
			}
			else if (csm->getFdrInfo()->getQValue() < qValueThreshold)
			{
				qValueThreshold = csm->getFdrInfo()->getQValue();
			}
		}
	}

	Crosslinker *XLSearchTask::GenerateUserDefinedCrosslinker(TaskLayer::XlSearchParameters *xlSearchParameters)
	{
		Nullable<double> tempVar = xlSearchParameters.getCrosslinkerTotalMass();
		Nullable<double> tempVar2 = xlSearchParameters.getCrosslinkerShortMass();
		Nullable<double> tempVar3 = xlSearchParameters.getCrosslinkerLongMass();
		Nullable<double> tempVar4 = xlSearchParameters.getCrosslinkerLoopMass();
		Nullable<double> tempVar5 = xlSearchParameters.getCrosslinkerDeadEndMassH2O();
		Nullable<double> tempVar6 = xlSearchParameters.getCrosslinkerDeadEndMassNH2();
		Nullable<double> tempVar7 = xlSearchParameters.getCrosslinkerDeadEndMassTris();
		auto crosslinker = new Crosslinker(xlSearchParameters->getCrosslinkerResidues(), xlSearchParameters->getCrosslinkerResidues2(), xlSearchParameters->getCrosslinkerName(), xlSearchParameters->getIsCleavable(), tempVar ? tempVar : NAN, tempVar2 ? tempVar2 : NAN, tempVar3 ? tempVar3 : NAN, tempVar4 ? tempVar4 : NAN, tempVar5 ? tempVar5 : NAN, tempVar6 ? tempVar6 : NAN, tempVar7 ? tempVar7 : NAN);

//C# TO C++ CONVERTER TODO TASK: A 'delete crosslinker' statement was not added since crosslinker was used in a 'return' or 'throw' statement.
		return crosslinker;
	}

	void XLSearchTask::WritePsmCrossToTsv(std::vector<CrosslinkSpectralMatch*> &items, const std::string &filePath, int writeType)
	{
		if (items.empty())
		{
			return;
		}

//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(filePath))
		{
			StreamWriter output = StreamWriter(filePath);
			std::string header = "";
			switch (writeType)
			{
				case 1:
					header = CrosslinkSpectralMatch::GetTabSepHeaderSingle();
					break;
				case 2:
					header = CrosslinkSpectralMatch::GetTabSepHeaderCross();
					break;
				default:
					break;
			}
			output.WriteLine(header);
			for (auto heh : items)
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				output.WriteLine(heh->ToString());
			}
		}
	}

	void XLSearchTask::WriteCrosslinkToTxtForPercolator(std::vector<CrosslinkSpectralMatch*> &items, const std::string &outputFolder, const std::string &fileName, Crosslinker *crosslinker, std::vector<std::string> &nestedIds)
	{
		if (items.empty())
		{
			return;
		}
		auto writtenFile = FileSystem::combine(outputFolder, fileName + ".txt");
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(writtenFile))
		{
			StreamWriter output = StreamWriter(writtenFile);
			output.WriteLine(std::string("SpecId\tLabel\tScannr\tScore\tdScore\tNormRank\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum") + "\tPeptide\tProtein");
			for (auto item : items)
			{
				if (item->getBaseSequence() != "" && item->getBetaPeptide()->getBaseSequence() != "" && item->getProteinAccession() != "" && item->getBetaPeptide()->getProteinAccession() != "")
				{
					std::string x = "T";
					int label = 1;
					if (item->getIsDecoy() || item->getBetaPeptide()->getIsDecoy())
					{
						x = "D";
						label = -1;
					}
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					output.WriteLine(x + "-" + item->getScanNumber().ToString(CultureInfo::InvariantCulture) + "-" + item->getScanRetentionTime().ToString(CultureInfo::InvariantCulture) + "\t" + label.ToString(CultureInfo::InvariantCulture) + "\t" + item->getScanNumber().ToString(CultureInfo::InvariantCulture) + "\t" + item->getXLTotalScore().ToString(CultureInfo::InvariantCulture) + "\t" + item->getDeltaScore().ToString(CultureInfo::InvariantCulture) + "\t" + (item->getXlRank()[0] + item->getXlRank()[1]).ToString(CultureInfo::InvariantCulture) + "\t" + item->getScanPrecursorCharge().ToString(CultureInfo::InvariantCulture) + "\t" + item->getScanPrecursorMass().ToString(CultureInfo::InvariantCulture) + "\t" + ((item->getPeptideMonisotopicMass().HasValue && item->getBetaPeptide()->getPeptideMonisotopicMass().HasValue) ? ((item->getScanPrecursorMass() - item->getBetaPeptide()->getPeptideMonisotopicMass().Value - item->getPeptideMonisotopicMass().Value - crosslinker->getTotalMass()) / item->getScanPrecursorMass() * 1E6).ToString(CultureInfo::InvariantCulture) : "---") + "\t" + item->getBetaPeptide()->getBaseSequence().length().ToString(CultureInfo::InvariantCulture) + "\t" + item->getBaseSequence().length().ToString(CultureInfo::InvariantCulture) + "\t" + (item->getBetaPeptide()->getBaseSequence().length() + item->getBaseSequence().length()).ToString(CultureInfo::InvariantCulture) + "\t" + "-." + item->getBaseSequence() + item->getLinkPositions().front().ToString(CultureInfo::InvariantCulture) + "--" + item->getBetaPeptide()->getBaseSequence() + item->getBetaPeptide()->getLinkPositions().front().ToString(CultureInfo::InvariantCulture) + ".-" + "\t" + item->BestMatchingPeptides.First().Peptide.Protein.Accession.ToString(CultureInfo::InvariantCulture) + "(" + item->getXlProteinPos().ToString(CultureInfo::InvariantCulture) + ")" + "\t" + item->getBetaPeptide()->BestMatchingPeptides.First().Peptide.Protein.Accession.ToString(CultureInfo::InvariantCulture) + "(" + item->getBetaPeptide()->getXlProteinPos().ToString(CultureInfo::InvariantCulture) + ")");
				}
			}
		}
		FinishedWritingFile(writtenFile, nestedIds);
	}

	void XLSearchTask::WritePepXML_xl(std::vector<CrosslinkSpectralMatch*> &items, std::vector<Protein*> &proteinList, const std::string &databasePath, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<std::string> &localizeableModificationTypes, const std::string &outputFolder, const std::string &fileName, std::vector<std::string> &nestedIds)
	{
		if (!items.Any())
		{
			return;
		}

		XmlSerializer *_indexedSerializer = new XmlSerializer(pepXML::Generated::msms_pipeline_analysis::typeid);
		auto _pepxml = new pepXML::Generated::msms_pipeline_analysis();

		_pepxml->date = DateTime::Now;
		_pepxml->summary_xml = items[0]->getFullFilePath() + ".pep.XM";

		std::string proteaseC = "";
		std::string proteaseNC = "";
		for (auto x : getCommonParameters()->getDigestionParams()->Protease.DigestionMotifs->Select([&] (std::any m)
		{
			m::InducingCleavage;
		}))
		{
			proteaseC += x;
		}
		for (auto x : getCommonParameters()->getDigestionParams()->Protease.DigestionMotifs->Select([&] (std::any m)
		{
			m::PreventingCleavage;
		}))
		{
			proteaseNC += x;
		}

		Crosslinker tempVar();
		Crosslinker *crosslinker = (&tempVar)->SelectCrosslinker(getXlSearchParameters()->getCrosslinkerType());
		if (getXlSearchParameters()->getCrosslinkerType() == CrosslinkerType::UserDefined)
		{
			crosslinker = GenerateUserDefinedCrosslinker(getXlSearchParameters());
		}

		std::string fileNameNoExtension = Path::GetFileNameWithoutExtension(items[0]->getFullFilePath());
		std::string filePathNoExtension = Path::ChangeExtension(items[0]->getFullFilePath(), "");
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		std::string modSites = crosslinker->getCrosslinkerModSites().ToCharArray().Concat(crosslinker->getCrosslinkerModSites2().ToCharArray())->Distinct().ToString();

		auto para = std::vector<pepXML::Generated::nameValueType*>();
		{
			pepXML::Generated::nameValueType *tempVar2 = new pepXML::Generated::nameValueType();
			tempVar2->name = "threads";
			tempVar2->value = std::to_string(getCommonParameters()->getMaxThreadsToUsePerFile());
			para.push_back(tempVar2);
			pepXML::Generated::nameValueType *tempVar3 = new pepXML::Generated::nameValueType();
			tempVar3->name = "database";
			tempVar3->value = databasePath;
			para.push_back(tempVar3);
			pepXML::Generated::nameValueType *tempVar4 = new pepXML::Generated::nameValueType();
			tempVar4->name = "MS_data_file";
			tempVar4->value = items[0]->getFullFilePath();
			para.push_back(tempVar4);
			pepXML::Generated::nameValueType *tempVar5 = new pepXML::Generated::nameValueType();
			tempVar5->name = "Cross-link precursor Mass Tolerance";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar5->value = getCommonParameters()->getPrecursorMassTolerance()->ToString();
			para.push_back(tempVar5);
			pepXML::Generated::nameValueType *tempVar6 = new pepXML::Generated::nameValueType();
			tempVar6->name = "Cross-linker type";
			tempVar6->value = crosslinker->getCrosslinkerName();
			para.push_back(tempVar6);
			pepXML::Generated::nameValueType *tempVar7 = new pepXML::Generated::nameValueType();
			tempVar7->name = "Cross-linker mass";
			tempVar7->value = std::to_string(crosslinker->getTotalMass());
			para.push_back(tempVar7);
			pepXML::Generated::nameValueType *tempVar8 = new pepXML::Generated::nameValueType();
			tempVar8->name = "Cross-linker cleavable";
			tempVar8->value = StringHelper::toString(crosslinker->getCleavable());
			para.push_back(tempVar8);
			pepXML::Generated::nameValueType *tempVar9 = new pepXML::Generated::nameValueType();
			tempVar9->name = "Cross-linker cleavable long mass";
			tempVar9->value = std::to_string(crosslinker->getCleaveMassLong());
			para.push_back(tempVar9);
			pepXML::Generated::nameValueType *tempVar10 = new pepXML::Generated::nameValueType();
			tempVar10->name = "Cross-linker cleavable short mass";
			tempVar10->value = std::to_string(crosslinker->getCleaveMassShort());
			para.push_back(tempVar10);
			pepXML::Generated::nameValueType *tempVar11 = new pepXML::Generated::nameValueType();
			tempVar11->name = "Cross-linker xl site";
			tempVar11->value = modSites;
			para.push_back(tempVar11);

			pepXML::Generated::nameValueType *tempVar12 = new pepXML::Generated::nameValueType();
			tempVar12->name = "Generate decoy proteins";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar12->value = getXlSearchParameters()->getDecoyType()->ToString();
			para.push_back(tempVar12);
			pepXML::Generated::nameValueType *tempVar13 = new pepXML::Generated::nameValueType();
			tempVar13->name = "MaxMissed Cleavages";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar13->value = getCommonParameters()->getDigestionParams()->MaxMissedCleavages.ToString();
			para.push_back(tempVar13);
			pepXML::Generated::nameValueType *tempVar14 = new pepXML::Generated::nameValueType();
			tempVar14->name = "Protease";
			tempVar14->value = getCommonParameters()->getDigestionParams()->Protease->Name;
			para.push_back(tempVar14);
			pepXML::Generated::nameValueType *tempVar15 = new pepXML::Generated::nameValueType();
			tempVar15->name = "Initiator Methionine";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar15->value = getCommonParameters()->getDigestionParams()->InitiatorMethionineBehavior.ToString();
			para.push_back(tempVar15);
			pepXML::Generated::nameValueType *tempVar16 = new pepXML::Generated::nameValueType();
			tempVar16->name = "Max Modification Isoforms";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar16->value = getCommonParameters()->getDigestionParams()->MaxModificationIsoforms.ToString();
			para.push_back(tempVar16);
			pepXML::Generated::nameValueType *tempVar17 = new pepXML::Generated::nameValueType();
			tempVar17->name = "Min Peptide Len";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar17->value = getCommonParameters()->getDigestionParams()->MinPeptideLength.ToString();
			para.push_back(tempVar17);
			pepXML::Generated::nameValueType *tempVar18 = new pepXML::Generated::nameValueType();
			tempVar18->name = "Max Peptide Len";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar18->value = getCommonParameters()->getDigestionParams()->MaxPeptideLength.ToString();
			para.push_back(tempVar18);
			pepXML::Generated::nameValueType *tempVar19 = new pepXML::Generated::nameValueType();
			tempVar19->name = "Product Mass Tolerance";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar19->value = getCommonParameters()->getProductMassTolerance()->ToString();
			para.push_back(tempVar19);
			pepXML::Generated::nameValueType *tempVar20 = new pepXML::Generated::nameValueType();
			tempVar20->name = "Ions to search";
			tempVar20->value = std::string::Join(", ", DissociationTypeCollection::ProductsFromDissociationType[getCommonParameters()->getDissociationType()]);
			para.push_back(tempVar20);

			for (auto fixedMod : fixedModifications)
			{
				pepXML::Generated::nameValueType *tempVar21 = new pepXML::Generated::nameValueType();
				tempVar21->name = "Fixed Modifications: " + fixedMod->IdWithMotif;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar21->value = fixedMod->MonoisotopicMass.ToString();
				para.push_back(tempVar21);
			}
			for (auto variableMod : variableModifications)
			{
				pepXML::Generated::nameValueType *tempVar22 = new pepXML::Generated::nameValueType();
				tempVar22->name = "Variable Modifications: " + variableMod->IdWithMotif;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar22->value = variableMod->MonoisotopicMass.ToString();
				para.push_back(tempVar22);
			}

			pepXML::Generated::nameValueType *tempVar23 = new pepXML::Generated::nameValueType();
			tempVar23->name = "Localize All Modifications";
			tempVar23->value = "true";
			para.push_back(tempVar23);
		}

		pepXML::Generated::msms_pipeline_analysisMsms_run_summary *tempVar24 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summary();
		tempVar24->base_name = filePathNoExtension;
		tempVar24->raw_data_type = "raw";
		tempVar24->raw_data = ".mzM";
		tempVar24->sample_enzyme = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySample_enzyme();
		tempVar24->sample_enzyme->name = getCommonParameters()->getDigestionParams()->Protease->Name;
		pepXML::Generated::msms_pipeline_analysisMsms_run_summarySample_enzymeSpecificity *tempVar25 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySample_enzymeSpecificity();
		tempVar25->cut = proteaseC;
		tempVar25->no_cut = proteaseNC;
		tempVar24->sample_enzyme->specificity = {tempVar25};
		pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summary *tempVar26 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summary();
		tempVar26->base_name = filePathNoExtension;
		tempVar26->search_engine_version = GlobalVariables::getMetaMorpheusVersion();
		tempVar26->precursor_mass_type = pepXML::Generated::massType::monoisotopic;
		tempVar26->fragment_mass_type = pepXML::Generated::massType::monoisotopic;
		tempVar26->search_id = 1;
		tempVar26->search_database = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summarySearch_database();
		tempVar26->search_database->local_path = databasePath;
		tempVar26->search_database->type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summarySearch_databaseType::AA;
		tempVar26->enzymatic_search_constraint = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySearch_summaryEnzymatic_search_constraint();
		tempVar26->enzymatic_search_constraint->enzyme = getCommonParameters()->getDigestionParams()->Protease->Name;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		tempVar26->enzymatic_search_constraint->max_num_internal_cleavages = getCommonParameters()->getDigestionParams()->MaxMissedCleavages.ToString();
		tempVar26->parameter = para.ToArray();
		tempVar24->search_summary = {tempVar26};
		_pepxml->msms_run_summary = {tempVar24};

		_pepxml->msms_run_summary[0]->spectrum_query = std::vector<pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_query*>(items.size());

		auto searchHits = std::vector<pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit*>();
		for (int i = 0; i < items.size(); i++)
		{
			auto mods = std::vector<pepXML::Generated::modInfoDataTypeMod_aminoacid_mass*>();
			auto alphaPeptide = items[i]->BestMatchingPeptides.First().Peptide;

			for (auto modification : alphaPeptide->AllModsOneIsNterminus)
			{
				auto mod = new pepXML::Generated::modInfoDataTypeMod_aminoacid_mass();
				mod->mass = modification->Value->MonoisotopicMass->Value;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				mod->position = (modification->Key - 1).ToString();
				mods.push_back(mod);

//C# TO C++ CONVERTER TODO TASK: A 'delete mod' statement was not added since mod was passed to a method or constructor. Handle memory management manually.
			}

			if (items[i]->getCrossType() == PsmCrossType::Single)
			{
				auto searchHit = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit();
				searchHit->hit_rank = 1;
				searchHit->peptide = alphaPeptide->BaseSequence;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				searchHit->peptide_prev_aa = alphaPeptide->PreviousAminoAcid.ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				searchHit->peptide_next_aa = alphaPeptide->NextAminoAcid.ToString();
				searchHit->protein = alphaPeptide->Protein.Accession;
				searchHit->num_tot_proteins = 1;
				searchHit->calc_neutral_pep_mass = static_cast<float>(items[i]->getScanPrecursorMass());
				searchHit->massdiff = std::to_string(items[i]->getScanPrecursorMass() - items[i]->getPeptideMonisotopicMass()->Value);
				searchHit->xlink_typeSpecified = true;
				searchHit->xlink_type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type::na;
				searchHit->modification_info = new pepXML::Generated::modInfoDataType();
				searchHit->modification_info->mod_aminoacid_mass = mods.ToArray();
				pepXML::Generated::nameValueType *tempVar27 = new pepXML::Generated::nameValueType();
				tempVar27->name = "xlTotalScore";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar27->value = items[i]->getXLTotalScore().ToString();
				pepXML::Generated::nameValueType *tempVar28 = new pepXML::Generated::nameValueType();
				tempVar28->name = "Qvalue";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar28->value = items[i]->getFdrInfo().getQValue().ToString();
				searchHit->search_score = {tempVar27, tempVar28};
				searchHits.push_back(searchHit);

//C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since searchHit was passed to a method or constructor. Handle memory management manually.
			}
			else if (items[i]->getCrossType() == PsmCrossType::DeadEnd || items[i]->getCrossType() == PsmCrossType::DeadEndH2O || items[i]->getCrossType() == PsmCrossType::DeadEndNH2 || items[i]->getCrossType() == PsmCrossType::DeadEndTris)
			{
				double crosslinkerDeadEndMass = 0;
				switch (items[i]->getCrossType())
				{
					case PsmCrossType::DeadEndNH2:
						crosslinkerDeadEndMass = crosslinker->getDeadendMassNH2();
						break;

					case PsmCrossType::DeadEndTris:
						crosslinkerDeadEndMass = crosslinker->getDeadendMassTris();
						break;

					default:
						crosslinkerDeadEndMass = crosslinker->getDeadendMassH2O();
						break;
				}
				auto mod = new pepXML::Generated::modInfoDataTypeMod_aminoacid_mass();
				mod->mass = crosslinkerDeadEndMass;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				mod->position = items[i]->getLinkPositions().front().ToString();
				mods.push_back(mod);
				auto searchHit = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit();
				searchHit->hit_rank = 1;
				searchHit->peptide = alphaPeptide->BaseSequence;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				searchHit->peptide_prev_aa = alphaPeptide->PreviousAminoAcid.ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				searchHit->peptide_next_aa = alphaPeptide->NextAminoAcid.ToString();
				searchHit->protein = alphaPeptide->Protein.Accession;
				searchHit->num_tot_proteins = 1;
				searchHit->calc_neutral_pep_mass = static_cast<float>(items[i]->getScanPrecursorMass());
				searchHit->massdiff = std::to_string(items[i]->getScanPrecursorMass() - items[i]->getPeptideMonisotopicMass()->Value - crosslinkerDeadEndMass);
				searchHit->xlink_typeSpecified = true;
				searchHit->xlink_type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type::na;
				searchHit->modification_info = new pepXML::Generated::modInfoDataType();
				searchHit->modification_info->mod_aminoacid_mass = mods.ToArray();
				pepXML::Generated::nameValueType *tempVar29 = new pepXML::Generated::nameValueType();
				tempVar29->name = "xlTotalScore";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar29->value = items[i]->getXLTotalScore().ToString();
				pepXML::Generated::nameValueType *tempVar30 = new pepXML::Generated::nameValueType();
				tempVar30->name = "Qvalue";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar30->value = items[i]->getFdrInfo().getQValue().ToString();
				searchHit->search_score = {tempVar29, tempVar30};
				searchHits.push_back(searchHit);

//C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since searchHit was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete mod' statement was not added since mod was passed to a method or constructor. Handle memory management manually.
			}
			else if (items[i]->getCrossType() == PsmCrossType::Inter || items[i]->getCrossType() == PsmCrossType::Intra || items[i]->getCrossType() == PsmCrossType::Cross)
			{
				auto betaPeptide = items[i]->getBetaPeptide().BestMatchingPeptides.First().Peptide;
				auto modsBeta = std::vector<pepXML::Generated::modInfoDataTypeMod_aminoacid_mass*>();

				for (auto mod : betaPeptide->AllModsOneIsNterminus)
				{
					auto modBeta = new pepXML::Generated::modInfoDataTypeMod_aminoacid_mass();
					modBeta->mass = mod->Value->MonoisotopicMass->Value;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					modBeta->position = (mod->Key - 1).ToString();
					modsBeta.push_back(modBeta);

//C# TO C++ CONVERTER TODO TASK: A 'delete modBeta' statement was not added since modBeta was passed to a method or constructor. Handle memory management manually.
				}

				auto alpha = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide();
				alpha->peptide = alphaPeptide->BaseSequence;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				alpha->peptide_prev_aa = alphaPeptide->PreviousAminoAcid.ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				alpha->peptide_next_aa = alphaPeptide->NextAminoAcid.ToString();
				alpha->protein = alphaPeptide->Protein.Accession;
				alpha->num_tot_proteins = 1;
				alpha->calc_neutral_pep_mass = static_cast<float>(items[i]->getPeptideMonisotopicMass()->Value);
				alpha->complement_mass = static_cast<float>(items[i]->getScanPrecursorMass() - alphaPeptide->MonoisotopicMass);
				alpha->designation = "alpha";
				alpha->modification_info = new pepXML::Generated::modInfoDataType();
				alpha->modification_info->mod_aminoacid_mass = mods.ToArray();
				pepXML::Generated::nameValueType *tempVar31 = new pepXML::Generated::nameValueType();
				tempVar31->name = "xlscore";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar31->value = items[i]->getXLTotalScore().ToString();
				pepXML::Generated::nameValueType *tempVar32 = new pepXML::Generated::nameValueType();
				tempVar32->name = "link";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar32->value = items[i]->getLinkPositions().front().ToString();
				alpha->xlink_score = {tempVar31, tempVar32};
				auto beta = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide();
				beta->peptide = betaPeptide->BaseSequence;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				beta->peptide_prev_aa = betaPeptide->PreviousAminoAcid.ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				beta->peptide_next_aa = betaPeptide->NextAminoAcid.ToString();
				beta->protein = betaPeptide->Protein.Accession;
				beta->num_tot_proteins = 1;
				beta->calc_neutral_pep_mass = static_cast<float>(betaPeptide->MonoisotopicMass);
				beta->complement_mass = static_cast<float>(items[i]->getScanPrecursorMass() - betaPeptide->MonoisotopicMass);
				beta->designation = "beta";
				beta->modification_info = new pepXML::Generated::modInfoDataType();
				beta->modification_info->mod_aminoacid_mass = modsBeta.ToArray();
				pepXML::Generated::nameValueType *tempVar33 = new pepXML::Generated::nameValueType();
				tempVar33->name = "xlscore";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar33->value = items[i]->getBetaPeptide().getScore().ToString();
				pepXML::Generated::nameValueType *tempVar34 = new pepXML::Generated::nameValueType();
				tempVar34->name = "link";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar34->value = items[i]->getBetaPeptide().getLinkPositions().front().ToString();
				beta->xlink_score = {tempVar33, tempVar34};
				auto cross = {alpha, beta};
				auto searchHit = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit();
				searchHit->hit_rank = 1;
				searchHit->peptide = "-";
				searchHit->peptide_prev_aa = "-";
				searchHit->peptide_next_aa = "-";
				searchHit->protein = "-";
				searchHit->num_tot_proteins = 1;
				searchHit->calc_neutral_pep_mass = static_cast<float>(items[i]->getScanPrecursorMass());
				searchHit->massdiff = std::to_string(items[i]->getScanPrecursorMass() - betaPeptide->MonoisotopicMass - alphaPeptide->MonoisotopicMass - crosslinker->getTotalMass());
				searchHit->xlink_typeSpecified = true;
				searchHit->xlink_type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type::xl;
				searchHit->xlink = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink();
				searchHit->xlink->identifier = crosslinker->getCrosslinkerName();
				searchHit->xlink->mass = static_cast<float>(crosslinker->getTotalMass());
				searchHit->xlink->linked_peptide = cross;
				pepXML::Generated::nameValueType *tempVar35 = new pepXML::Generated::nameValueType();
				tempVar35->name = "xlTotalScore";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar35->value = items[i]->getXLTotalScore().ToString();
				pepXML::Generated::nameValueType *tempVar36 = new pepXML::Generated::nameValueType();
				tempVar36->name = "Qvalue";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar36->value = items[i]->getFdrInfo().getQValue().ToString();
				searchHit->search_score = {tempVar35, tempVar36};
				searchHits.push_back(searchHit);

//C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since searchHit was passed to a method or constructor. Handle memory management manually.
				delete beta;
				delete alpha;
			}
			else if (items[i]->getCrossType() == PsmCrossType::Loop)
			{
				auto thePeptide = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide();
				pepXML::Generated::nameValueType *tempVar37 = new pepXML::Generated::nameValueType();
				tempVar37->name = "link";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar37->value = items[i]->getLinkPositions().front().ToString();
				pepXML::Generated::nameValueType *tempVar38 = new pepXML::Generated::nameValueType();
				tempVar38->name = "link";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar38->value = items[i]->getLinkPositions()[1].ToString();
				thePeptide->xlink_score = {tempVar37, tempVar38};
				auto cross = {thePeptide};
				auto searchHit = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit();
				searchHit->hit_rank = 1;
				searchHit->peptide = alphaPeptide->BaseSequence;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				searchHit->peptide_prev_aa = alphaPeptide->PreviousAminoAcid.ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				searchHit->peptide_next_aa = alphaPeptide->NextAminoAcid.ToString();
				searchHit->protein = alphaPeptide->Protein.Accession;
				searchHit->num_tot_proteins = 1;
				searchHit->calc_neutral_pep_mass = static_cast<float>(items[i]->getScanPrecursorMass());
				searchHit->massdiff = std::to_string(items[i]->getScanPrecursorMass() - alphaPeptide->MonoisotopicMass - crosslinker->getLoopMass());
				searchHit->xlink_typeSpecified = true;
				searchHit->xlink_type = pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type::loop;
				searchHit->modification_info = new pepXML::Generated::modInfoDataType();
				searchHit->modification_info->mod_aminoacid_mass = mods.ToArray();
				searchHit->xlink = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink();
				searchHit->xlink->identifier = crosslinker->getCrosslinkerName();
				searchHit->xlink->mass = static_cast<float>(crosslinker->getTotalMass());
				searchHit->xlink->linked_peptide = cross;
				pepXML::Generated::nameValueType *tempVar39 = new pepXML::Generated::nameValueType();
				tempVar39->name = "xlTotalScore";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar39->value = items[i]->getXLTotalScore().ToString();
				pepXML::Generated::nameValueType *tempVar40 = new pepXML::Generated::nameValueType();
				tempVar40->name = "Qvalue";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				tempVar40->value = items[i]->getFdrInfo().getQValue().ToString();
				searchHit->search_score = {tempVar39, tempVar40};
				searchHits.push_back(searchHit);

//C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since searchHit was passed to a method or constructor. Handle memory management manually.
				delete thePeptide;
			}
		}

		for (int i = 0; i < items.size(); i++)
		{
			pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_query *tempVar41 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_query();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar41->spectrum = fileNameNoExtension + "." + items[i]->getScanNumber().ToString();
			tempVar41->start_scan = static_cast<unsigned int>(items[i]->getScanNumber());
			tempVar41->end_scan = static_cast<unsigned int>(items[i]->getScanNumber());
			tempVar41->precursor_neutral_mass = static_cast<float>(items[i]->getScanPrecursorMass());
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar41->assumed_charge = items[i]->getScanPrecursorCharge().ToString();
			tempVar41->index = static_cast<unsigned int>(i + 1);
			tempVar41->retention_time_sec = static_cast<float>(items[i]->getScanRetentionTime() * 60);
			pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_result *tempVar42 = new pepXML::Generated::msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_result();
			tempVar42->search_hit = {searchHits[i]};
			tempVar41->search_result = {tempVar42};
			_pepxml->msms_run_summary[0].spectrum_query[i] = tempVar41;
		}

		TextWriter *writer = new StreamWriter(FileSystem::combine(outputFolder, fileName + ".pep.XM"));
		_indexedSerializer->Serialize(writer, _pepxml);
		writer->Close();
		FinishedWritingFile(FileSystem::combine(outputFolder, fileName + ".pep.XM"), nestedIds);

//C# TO C++ CONVERTER TODO TASK: A 'delete writer' statement was not added since writer was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete _pepxml' statement was not added since _pepxml was passed to a method or constructor. Handle memory management manually.
		delete _indexedSerializer;
	}
}
