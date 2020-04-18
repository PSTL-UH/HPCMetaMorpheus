#include "PostSearchAnalysisTask.h"
#include "PostSearchAnalysisParameters.h"
#include "../../EngineLayer/ProteinParsimony/ProteinGroup.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../MyTaskResults.h"
#include "../../EngineLayer/GlobalVariables.h"
#include "../../EngineLayer/EventArgs/ProgressEventArgs.h"
#include "../DbForTask.h"
#include "../FileSpecificParameters.h"
#include "../../EngineLayer/CommonParameters.h"
#include "../../EngineLayer/ProteinParsimony/ProteinParsimonyEngine.h"
#include "../../EngineLayer/ProteinParsimony/ProteinParsimonyResults.h"
#include "../../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.h"
#include "../../EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrResults.h"
#include "../../EngineLayer/MetaMorpheusException.h"
#include "MzIdentMLWriter.h"
#include "../PepXMLWriter.h"

using namespace EngineLayer;
using namespace EngineLayer::FdrAnalysis;
using namespace EngineLayer::HistogramAnalysis;
using namespace EngineLayer::Localization;
using namespace EngineLayer::ModificationAnalysis;
using namespace FlashLFQ;
using namespace MassSpectrometry;
//using namespace MathNet::Numerics::Distributions;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{

	PostSearchAnalysisParameters *PostSearchAnalysisTask::getParameters() const
	{
		return privateParameters;
	}

	void PostSearchAnalysisTask::setParameters(PostSearchAnalysisParameters *value)
	{
		privateParameters = value;
	}

	std::vector<EngineLayer::ProteinGroup*> PostSearchAnalysisTask::getProteinGroups() const
	{
		return privateProteinGroups;
	}

	void PostSearchAnalysisTask::setProteinGroups(const std::vector<EngineLayer::ProteinGroup*> &value)
	{
		privateProteinGroups = value;
	}

	std::vector<IGrouping<std::string, PeptideSpectralMatch*>*> PostSearchAnalysisTask::getPsmsGroupedByFile() const
	{
		return privatePsmsGroupedByFile;
	}

	void PostSearchAnalysisTask::setPsmsGroupedByFile(const std::vector<IGrouping<std::string, PeptideSpectralMatch*>*> &value)
	{
		privatePsmsGroupedByFile = value;
	}

	PostSearchAnalysisTask::PostSearchAnalysisTask() : MetaMorpheusTask(MyTask::Search)
	{
	}

	MyTaskResults *PostSearchAnalysisTask::Run()
	{
		// Stop loop if canceled
		if (GlobalVariables::getStopLoops())
		{
			return getParameters()->getSearchTaskResults();
		}

		if ( getParameters()->getSearchParameters()->getMassDiffAcceptorType() == MassDiffAcceptorType::ModOpen ||
                     getParameters()->getSearchParameters()->getMassDiffAcceptorType() == MassDiffAcceptorType::Open ||
                     getParameters()->getSearchParameters()->getMassDiffAcceptorType() == MassDiffAcceptorType::Custom)
		{
                    // This only makes sense if there is a mass difference that you want to localize.
                    //No use for exact and missed monoisotopic mass searches.
                    getParameters()->getSearchParameters()->setDoLocalizationAnalysis(true);
		}
		else
		{
                    getParameters()->getSearchParameters()->setDoLocalizationAnalysis(false);
		}

		//update all psms with peptide info
                //if it hasn't been done already
		if (getParameters()->getSearchParameters()->getSearchType() != SearchType::NonSpecific) 
		{
                    getParameters()->setAllPsms(getParameters()->getAllPsms().Where([&] (std::any psm){
				return psm != nullptr;
                            }).ToList());
                    std::for_each(getParameters()->AllPsms.begin(), getParameters()->AllPsms.end(), [&] (std::any psm)	{
                            psm::ResolveAllAmbiguities();
			});
                    
                    getParameters()->setAllPsms(getParameters()->getAllPsms().OrderByDescending([&] (std::any b) {
				b::Score;
                            }).ThenBy([&] (std::any b)	{
                                    b::PeptideMonisotopicMass.HasValue ? std::abs(b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value) : std::numeric_limits<double>::max();
                                }).GroupBy([&] (std::any b) {
				(b::FullFilePath, b::ScanNumber, b::PeptideMonisotopicMass);
                                    })->Select([&] (std::any b)  {
                                            b::First();
                                        }).ToList());
                    
                    CalculatePsmFdr();
		}

		DoMassDifferenceLocalizationAnalysis();
		ProteinAnalysis();
		QuantificationAnalysis();

		ProgressEventArgs tempVar(100, "Done!", new std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files"});
		ReportProgress(&tempVar);

		HistogramAnalysis();
		WritePsmResults();
		WriteProteinResults();
		WriteQuantificationResults();
		WritePrunedDatabase();

		return getParameters()->getSearchTaskResults();
	}

	MyTaskResults *PostSearchAnalysisTask::RunSpecific(const std::string &OutputFolder,
                                                           std::vector<DbForTask*> &dbFilenameList,
                                                           std::vector<std::string> &currentRawFileList,
                                                           const std::string &taskId,
                                                           std::vector<FileSpecificParameters*> &fileSettingsList)
	{
		return nullptr;
	}

	void PostSearchAnalysisTask::CalculatePsmFdr()
	{
		// TODO: because FDR is done before parsimony, if a PSM matches to a target and a decoy protein, there may be conflicts between how it's handled in parsimony and the FDR engine here
		// for example, here it may be treated as a decoy PSM, where as in parsimony it will be determined by the parsimony algorithm which is agnostic of target/decoy assignments
		// this could cause weird PSM FDR issues

		Status("Estimating PSM FDR...", getParameters()->getSearchTaskId());
		int massDiffAcceptorNumNotches = getParameters()->getNumNotches();
		FdrAnalysisEngine tempVar(getParameters()->getAllPsms(), massDiffAcceptorNumNotches, getCommonParameters(), new std::vector<std::string> {getParameters()->getSearchTaskId()});
		(&tempVar)->Run();

		// sort by q-value because of group FDR stuff
		// e.g. multiprotease FDR, non/semi-specific protease, etc
		getParameters()->setAllPsms(getParameters()->getAllPsms().OrderBy([&] (std::any p)
		{
			p::FdrInfo::QValue;
		}).ThenByDescending([&] (std::any p)
		{
			p::Score;
		}).ThenBy([&] (std::any p)
		{
			p::FdrInfo::CumulativeTarget;
		}).ToList());

		Status("Done estimating PSM FDR!", getParameters()->getSearchTaskId());
	}

	void PostSearchAnalysisTask::ProteinAnalysis()
	{
		if (!getParameters()->getSearchParameters()->getDoParsimony())
		{
			return;
		}

		Status("Constructing protein groups...", getParameters()->getSearchTaskId());

		// run parsimony
		ProteinParsimonyEngine tempVar(getParameters()->getAllPsms(), getParameters()->getSearchParameters()->getModPeptidesAreDifferent(), getCommonParameters(), new std::vector<std::string> {getParameters()->getSearchTaskId()});
		ProteinParsimonyResults *proteinAnalysisResults = static_cast<ProteinParsimonyResults*>((&tempVar)->Run());

		// score protein groups and calculate FDR
		ProteinScoringAndFdrEngine tempVar2(proteinAnalysisResults->getProteinGroups(), getParameters()->getAllPsms(), getParameters()->getSearchParameters()->getNoOneHitWonders(), getParameters()->getSearchParameters()->getModPeptidesAreDifferent(), true, getCommonParameters(), new std::vector<std::string> {getParameters()->getSearchTaskId()});
		ProteinScoringAndFdrResults *proteinScoringAndFdrResults = static_cast<ProteinScoringAndFdrResults*>((&tempVar2)->Run());

		setProteinGroups(proteinScoringAndFdrResults->SortedAndScoredProteinGroups);

		for (auto psm : getParameters()->getAllPsms())
		{
			psm->ResolveAllAmbiguities();
		}

		Status("Done constructing protein groups!", getParameters()->getSearchTaskId());
	}

	void PostSearchAnalysisTask::DoMassDifferenceLocalizationAnalysis()
	{
		if (getParameters()->getSearchParameters()->getDoLocalizationAnalysis())
		{
			Status("Running mass-difference localization analysis...", getParameters()->getSearchTaskId());
			for (int spectraFileIndex = 0; spectraFileIndex < getParameters()->getCurrentRawFileList().size(); spectraFileIndex++)
			{
				EngineLayer::CommonParameters *combinedParams = SetAllFileSpecificCommonParams(getCommonParameters(), getParameters()->getFileSettingsList()[spectraFileIndex]);

				auto origDataFile = getParameters()->getCurrentRawFileList()[spectraFileIndex];
				Status("Running mass-difference localization analysis...", std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", origDataFile});
				MsDataFile *myMsDataFile = getParameters()->getMyFileManager()->LoadFile(origDataFile, std::make_optional(combinedParams->getTopNpeaks()), std::make_optional(combinedParams->getMinRatio()), combinedParams->getTrimMs1Peaks(), combinedParams->getTrimMsMsPeaks(), combinedParams);
				LocalizationEngine tempVar(getParameters()->getAllPsms().Where([&] (std::any b)
				{
					b::FullFilePath->Equals(origDataFile);
				}).ToList(), myMsDataFile, combinedParams, new std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", origDataFile});
				(&tempVar)->Run();
				getParameters()->getMyFileManager()->DoneWithFile(origDataFile);
				ProgressEventArgs tempVar2(100, "Done with localization analysis!", new std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", origDataFile});
				ReportProgress(&tempVar2);
			}
		}

		// count different modifications observed
		ModificationAnalysisEngine tempVar3(getParameters()->getAllPsms(), getCommonParameters(), new std::vector<std::string> {getParameters()->getSearchTaskId()});
		(&tempVar3)->Run();
	}

	void PostSearchAnalysisTask::QuantificationAnalysis()
	{
		if (!getParameters()->getSearchParameters()->getDoQuantification())
		{
			return;
		}

		// pass quantification parameters to FlashLFQ
		Status("Quantifying...", getParameters()->getSearchTaskId());

		// construct file info for FlashLFQ
		auto spectraFileInfo = std::vector<SpectraFileInfo*>();

		// get experimental design info for normalization
		if (getParameters()->getSearchParameters()->getNormalize())
		{
			std::string assumedExperimentalDesignPath = Directory::GetParent(getParameters()->getCurrentRawFileList().front())->FullName;
			assumedExperimentalDesignPath = FileSystem::combine(assumedExperimentalDesignPath, GlobalVariables::getExperimentalDesignFileName());

			if (FileSystem::fileExists(assumedExperimentalDesignPath))
			{
				auto experimentalDesign = File::ReadAllLines(assumedExperimentalDesignPath).ToDictionary([&] (std::any p)
				{
					p->Split(L'\t')[0];
				}, [&] (std::any p)
				{
					return p;
				});

				for (auto file : getParameters()->getCurrentRawFileList())
				{
					std::string filename = Path::GetFileNameWithoutExtension(file);

					auto expDesignForThisFile = experimentalDesign[filename];
					auto split = expDesignForThisFile.Split(L'\t');

					std::string condition = split[1];
					int biorep = std::stoi(split[2]);
					int fraction = std::stoi(split[3]);
					int techrep = std::stoi(split[4]);

					// experimental design info passed in here for each spectra file
					SpectraFileInfo tempVar(fullFilePathWithExtension: file, condition: condition, biorep: biorep - 1, fraction: fraction - 1, techrep: techrep - 1);
					spectraFileInfo.push_back(&tempVar);

					getParameters()->getMyFileManager()->DoneWithFile(file);
				}
			}
			else
			{
				throw MetaMorpheusException("Could not find experimental design file at location:\n" + assumedExperimentalDesignPath);
			}
		}
		else
		{
			for (auto file : getParameters()->getCurrentRawFileList())
			{
				// experimental design info passed in here for each spectra file
				SpectraFileInfo tempVar2(fullFilePathWithExtension: file, condition: "", biorep: 0, fraction: 0, techrep: 0);
				spectraFileInfo.push_back(&tempVar2);
				getParameters()->getMyFileManager()->DoneWithFile(file);
			}
		}

		// get PSMs to pass to FlashLFQ
		auto unambiguousPsmsBelowOnePercentFdr = getParameters()->getAllPsms().Where([&] (std::any p)
		{
			return p::FdrInfo::QValue <= 0.01 && p::FdrInfo::QValueNotch <= 0.01 && !p::IsDecoy && p::FullSequence != nullptr;
		}).ToList();

		auto psmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy([&] (std::any p)
		{
			p::FullFilePath;
		});

		// pass protein group info for each PSM
		auto psmToProteinGroups = std::unordered_map<PeptideSpectralMatch*, std::vector<FlashLFQ::ProteinGroup*>>();
		if (getProteinGroups().size() > 0)
		{
			for (auto proteinGroup : getProteinGroups())
			{
				auto proteinsOrderedByAccession = proteinGroup->getProteins().OrderBy([&] (std::any p)
				{
					p::Accession;
				});

				auto flashLfqProteinGroup = new FlashLFQ::ProteinGroup(proteinGroup->getProteinGroupName(), std::string::Join("|", proteinsOrderedByAccession->Select([&] (std::any p)
				{
					p::GeneNames->Select([&] (std::any x)
					{
						x::Item2;
					}).FirstOrDefault();
				})), std::string::Join("|", proteinsOrderedByAccession->Select([&] (std::any p)
				{
					p::Organism;
				}).Distinct()));

				for (auto psm : proteinGroup->getAllPsmsBelowOnePercentFDR().Where([&] (std::any v)
				{
				delete flashLfqProteinGroup;
					return v::FullSequence != nullptr;
				}))
				{
					TValue flashLfqProteinGroups;
					std::unordered_map<PeptideSpectralMatch*, std::vector<FlashLFQ::ProteinGroup*>>::const_iterator psmToProteinGroups_iterator = psmToProteinGroups.find(psm);
					if (psmToProteinGroups_iterator != psmToProteinGroups.end())
					{
						flashLfqProteinGroups = psmToProteinGroups_iterator->second;
						flashLfqProteinGroups->Add(flashLfqProteinGroup);
					}
					else
					{
						flashLfqProteinGroups = psmToProteinGroups_iterator->second;
						psmToProteinGroups.emplace(psm, std::vector<FlashLFQ::ProteinGroup*> {flashLfqProteinGroup});
					}
				}

//C# TO C++ CONVERTER TODO TASK: A 'delete flashLfqProteinGroup' statement was not added since flashLfqProteinGroup was passed to a method or constructor. Handle memory management manually.
			}
		}
		else
		{
			// if protein groups were not constructed, just use accession numbers
			auto accessionToPg = std::unordered_map<std::string, FlashLFQ::ProteinGroup*>();
			for (auto psm : unambiguousPsmsBelowOnePercentFdr)
			{
				auto proteins = psm.BestMatchingPeptides->Select([&] (std::any b)
				{
					b::Peptide::Protein;
				}).Distinct();

				for (auto protein : proteins)
				{
					if (accessionToPg.find(protein->Accession) == accessionToPg.end())
					{
						FlashLFQ::ProteinGroup tempVar3(protein->Accession, std::string::Join("|", protein->GeneNames->Select([&] (std::any p)
						{
							p::Item2;
						}).Distinct()), protein->Organism);
						accessionToPg.emplace(protein->Accession, &tempVar3);
					}

					TValue proteinGroups;
					std::unordered_map<PeptideSpectralMatch*, std::vector<FlashLFQ::ProteinGroup*>>::const_iterator psmToProteinGroups_iterator = psmToProteinGroups.find(psm);
					if (psmToProteinGroups_iterator != psmToProteinGroups.end())
					{
						proteinGroups = psmToProteinGroups_iterator->second;
						proteinGroups->Add(accessionToPg[protein->Accession]);
					}
					else
					{
						proteinGroups = psmToProteinGroups_iterator->second;
						psmToProteinGroups.emplace(psm, std::vector<std::vector<FlashLFQ::ProteinGroup*>>(protein->Accession) });
					}
				}
			}
		}

		// some PSMs may not have protein groups (if 2 peptides are required to construct a protein group, some PSMs will be left over)
		// the peptides should still be quantified but not considered for protein quantification
		auto undefinedPg = new FlashLFQ::ProteinGroup("UNDEFINED", "", "");
		//sort the unambiguous psms by protease to make MBR compatible with multiple proteases
		std::unordered_map<Protease*, std::vector<PeptideSpectralMatch*>> proteaseSortedPsms;
		std::unordered_map<Protease*, FlashLfqResults*> proteaseSortedFlashLFQResults;

		for (auto dp : getParameters()->getListOfDigestionParams())
		{
			if (proteaseSortedPsms.find(dp->Protease) == proteaseSortedPsms.end())
			{
				proteaseSortedPsms.emplace(dp->Protease, std::vector<PeptideSpectralMatch*>());
			}
		}
		for (auto psm : unambiguousPsmsBelowOnePercentFdr)
		{
			if (psmToProteinGroups.find(psm) == psmToProteinGroups.end())
			{
				psmToProteinGroups.emplace(psm, std::vector<FlashLFQ::ProteinGroup*> {undefinedPg});
			}

			proteaseSortedPsms[psm.DigestionParams::Protease].push_back(psm);
		}

		// pass PSM info to FlashLFQ
		auto flashLFQIdentifications = std::vector<Identification*>();
		for (auto spectraFile : psmsGroupedByFile)
		{
			auto rawfileinfo = spectraFileInfo.Where([&] (std::any p)
			{
				p::FullFilePathWithExtension->Equals(spectraFile->Key);
			}).First();

			for (auto psm : spectraFile)
			{
				Identification tempVar4(rawfileinfo, psm->BaseSequence, psm->FullSequence, psm->PeptideMonisotopicMass->Value, psm->ScanRetentionTime, psm->ScanPrecursorCharge, psmToProteinGroups[psm]);
				flashLFQIdentifications.push_back(&tempVar4);
			}
		}

		// run FlashLFQ
		auto FlashLfqEngine = new FlashLfqEngine(allIdentifications: flashLFQIdentifications, normalize: getParameters()->getSearchParameters()->getNormalize(), ppmTolerance: getParameters()->getSearchParameters()->getQuantifyPpmTol(), matchBetweenRuns: getParameters()->getSearchParameters()->getMatchBetweenRuns(), silent: true, optionalPeriodicTablePath: GlobalVariables::getElementsLocation(), maxThreads: getCommonParameters()->getMaxThreadsToUsePerFile());

		if (flashLFQIdentifications.Any())
		{
			getParameters()->setFlashLfqResults(FlashLfqEngine->Run());
		}

		//MultiProtease MBR capability code
		//Parameters.FlashLfqResults = null;

		//foreach (var proteasePsms in proteaseSortedPsms)
		//{
		//    var flashLFQIdentifications = new List<Identification>();
		//    var proteasePsmsGroupedByFile = proteasePsms.Value.GroupBy(p => p.FullFilePath);
		//    foreach (var spectraFile in proteasePsmsGroupedByFile)
		//    {
		//        var rawfileinfo = spectraFileInfo.Where(p => p.FullFilePathWithExtension.Equals(spectraFile.Key)).First();

		//        foreach (var psm in spectraFile)
		//        {
		//            flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence,
		//                psm.PeptideMonisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroups[psm]));
		//        }
		//    }

		//    // run FlashLFQ
		//    var FlashLfqEngine = new FlashLFQEngine(
		//        allIdentifications: flashLFQIdentifications,
		//        normalize: Parameters.SearchParameters.Normalize,
		//        ppmTolerance: Parameters.SearchParameters.QuantifyPpmTol,
		//        matchBetweenRuns: Parameters.SearchParameters.MatchBetweenRuns,
		//        silent: true,
		//        optionalPeriodicTablePath: GlobalVariables.ElementsLocation);

		//    if (flashLFQIdentifications.Any())
		//    {
		//        //make specific to protease
		//        var results = FlashLfqEngine.Run();

		//        if (Parameters.FlashLfqResults == null)
		//        {
		//            Parameters.FlashLfqResults = results;
		//        }
		//        else
		//        {
		//            Parameters.FlashLfqResults.MergeResultsWith(results);
		//        }
		//    }
		//}

		// get protein intensity back from FlashLFQ
		if (getProteinGroups().size() > 0 && getParameters()->getFlashLfqResults() != nullptr)
		{
			for (auto proteinGroup : getProteinGroups())
			{
				proteinGroup->setFilesForQuantification(spectraFileInfo);
				proteinGroup->setIntensitiesByFile(std::unordered_map<SpectraFileInfo*, double>());

				for (auto spectraFile : proteinGroup->getFilesForQuantification())
				{
					std::any flashLfqProteinGroup;
					if (getParameters()->getFlashLfqResults()->ProteinGroups.TryGetValue(proteinGroup->getProteinGroupName(), flashLfqProteinGroup))
					{
						proteinGroup->getIntensitiesByFile().emplace(spectraFile, flashLfqProteinGroup::GetIntensity(spectraFile));
					}
					else
					{
						proteinGroup->getIntensitiesByFile().emplace(spectraFile, 0);
					}
				}
			}
		}

		delete FlashLfqEngine;
//C# TO C++ CONVERTER TODO TASK: A 'delete undefinedPg' statement was not added since undefinedPg was passed to a method or constructor. Handle memory management manually.
	}

	void PostSearchAnalysisTask::HistogramAnalysis()
	{
		if (getParameters()->getSearchParameters()->getDoHistogramAnalysis())
		{
			auto limitedpsms_with_fdr = getParameters()->getAllPsms().Where([&] (std::any b)
			{
				(b::FdrInfo::QValue <= 0.01);
			}).ToList();
			if (limitedpsms_with_fdr.Any([&] (std::any b)
			{
				!b::IsDecoy;
			}))
			{
				Status("Running histogram analysis...", std::vector<std::string> {getParameters()->getSearchTaskId()});
				auto myTreeStructure = new BinTreeStructure();
				myTreeStructure->GenerateBins(limitedpsms_with_fdr, getParameters()->getSearchParameters()->getHistogramBinTolInDaltons());
				auto writtenFile = FileSystem::combine(getParameters()->getOutputFolder(), "MassDifferenceHistogram.tsv");
				WriteTree(myTreeStructure, writtenFile);
				FinishedWritingFile(writtenFile, std::vector<std::string> {getParameters()->getSearchTaskId()});

//C# TO C++ CONVERTER TODO TASK: A 'delete myTreeStructure' statement was not added since myTreeStructure was passed to a method or constructor. Handle memory management manually.
			}
		}
	}

	void PostSearchAnalysisTask::WritePsmResults()
	{
		Status("Writing results...", getParameters()->getSearchTaskId());
		std::vector<PeptideSpectralMatch*> filteredPsmListForOutput = getParameters()->getAllPsms().Where([&] (std::any p)
		{
			return p::FdrInfo::QValue <= getCommonParameters()->getQValueOutputFilter() && p::FdrInfo::QValueNotch <= getCommonParameters()->getQValueOutputFilter();
		}).ToList();

		if (!getParameters()->getSearchParameters()->getWriteDecoys())
		{
			filteredPsmListForOutput.RemoveAll([&] (std::any b)
			{
				b::IsDecoy;
			});
		}
		if (!getParameters()->getSearchParameters()->getWriteContaminants())
		{
			filteredPsmListForOutput.RemoveAll([&] (std::any b)
			{
				b::IsContaminant;
			});
		}

		// write PSMs
		std::string writtenFile = FileSystem::combine(getParameters()->getOutputFolder(), "AllPSMs.psmtsv");
		WritePsmsToTsv(filteredPsmListForOutput, writtenFile, getParameters()->getSearchParameters()->getModsToWriteSelection());
		FinishedWritingFile(writtenFile, std::vector<std::string> {getParameters()->getSearchTaskId()});

		// write PSMs for percolator
		writtenFile = FileSystem::combine(getParameters()->getOutputFolder(), "AllPSMs_FormattedForPercolator.tsv");
		WritePsmsForPercolator(filteredPsmListForOutput, writtenFile, getCommonParameters()->getQValueOutputFilter());
		FinishedWritingFile(writtenFile, std::vector<std::string> {getParameters()->getSearchTaskId()});

		// write best (highest-scoring) PSM per peptide
		writtenFile = FileSystem::combine(getParameters()->getOutputFolder(), "AllPeptides.psmtsv");
		std::vector<PeptideSpectralMatch*> peptides = getParameters()->getAllPsms().GroupBy([&] (std::any b)
		{
			b::FullSequence;
		})->Select([&] (std::any b)
		{
			b::FirstOrDefault();
		}).ToList();
		WritePsmsToTsv(filteredPsmListForOutput.GroupBy([&] (std::any b)
		{
			b::FullSequence;
		})->Select([&] (std::any b)
		{
			b::FirstOrDefault();
		}).ToList(), writtenFile, getParameters()->getSearchParameters()->getModsToWriteSelection());
		FinishedWritingFile(writtenFile, std::vector<std::string> {getParameters()->getSearchTaskId()});

		// write summary text
		getParameters()->getSearchTaskResults()->AddNiceText("All target PSMS within 1% FDR: " + getParameters()->getAllPsms().size()([&] (std::any a)
		{
			return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
		}));
		getParameters()->getSearchTaskResults()->AddNiceText("All target peptides within 1% FDR: " + peptides.size()([&] (std::any a)
		{
			return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
		}));
		if (getParameters()->getSearchParameters()->getDoParsimony())
		{
			getParameters()->getSearchTaskResults()->AddNiceText("All target protein groups within 1% FDR: " + std::to_string(getProteinGroups().size()([&] (std::any b)
			{
				return b::QValue <= 0.01 && !b::IsDecoy;
			})) + "\r\n");
		}

		setPsmsGroupedByFile(filteredPsmListForOutput.GroupBy([&] (std::any p)
		{
			p::FullFilePath;
		}));

		for (auto file : getPsmsGroupedByFile())
		{
			// write summary text
			auto psmsForThisFile = file->ToList();
			std::string strippedFileName = Path::GetFileNameWithoutExtension(file->First().FullFilePath);
			auto peptidesForFile = psmsForThisFile.GroupBy([&] (std::any b)
			{
				b::FullSequence;
			})->Select([&] (std::any b)
			{
				b::FirstOrDefault();
			}).ToList();

			getParameters()->getSearchTaskResults()->AddNiceText("MS2 spectra in " + strippedFileName + ": " + std::to_string(getParameters()->getNumMs2SpectraPerFile()[strippedFileName][0]));
			getParameters()->getSearchTaskResults()->AddNiceText("Precursors fragmented in " + strippedFileName + ": " + std::to_string(getParameters()->getNumMs2SpectraPerFile()[strippedFileName][1]));
			getParameters()->getSearchTaskResults()->AddNiceText("Target PSMs within 1% FDR in " + strippedFileName + ": " + psmsForThisFile.size()([&] (std::any a)
			{
				return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
			}));
			getParameters()->getSearchTaskResults()->AddNiceText("Target peptides within 1% FDR in " + strippedFileName + ": " + std::to_string(peptidesForFile.size()([&] (std::any a)
			{
				return a::FdrInfo::QValue <= 0.01 && !a::IsDecoy;
			})) + "\r\n");

			// writes all individual spectra file search results to subdirectory
			if (getParameters()->getCurrentRawFileList().size() > 1)
			{
				// create individual files subdirectory
				FileSystem::createDirectory(getParameters()->getIndividualResultsOutputFolder());

				// write PSMs
				writtenFile = FileSystem::combine(getParameters()->getIndividualResultsOutputFolder(), strippedFileName + "_PSMs.psmtsv");
				WritePsmsToTsv(psmsForThisFile, writtenFile, getParameters()->getSearchParameters()->getModsToWriteSelection());
				FinishedWritingFile(writtenFile, std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", file->First().FullFilePath});

				// write PSMs for percolator
				writtenFile = FileSystem::combine(getParameters()->getIndividualResultsOutputFolder(), strippedFileName + "_PSMsFormattedForPercolator.tsv");
				WritePsmsForPercolator(psmsForThisFile, writtenFile, getCommonParameters()->getQValueOutputFilter());
				FinishedWritingFile(writtenFile, std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", file->First().FullFilePath});

				// write best (highest-scoring) PSM per peptide
				writtenFile = FileSystem::combine(getParameters()->getIndividualResultsOutputFolder(), strippedFileName + "_Peptides.psmtsv");
				WritePsmsToTsv(peptidesForFile, writtenFile, getParameters()->getSearchParameters()->getModsToWriteSelection());
				FinishedWritingFile(writtenFile, std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", file->First().FullFilePath});
			}
		}
	}

	void PostSearchAnalysisTask::WriteProteinResults()
	{
		if (getParameters()->getSearchParameters()->getDoParsimony())
		{
			// write protein groups to tsv
			std::string writtenFile = FileSystem::combine(getParameters()->getOutputFolder(), "AllProteinGroups.tsv");
			WriteProteinGroupsToTsv(getProteinGroups(), writtenFile, {getParameters()->getSearchTaskId()}, getCommonParameters()->getQValueOutputFilter());

			// write all individual file results to subdirectory
			// local protein fdr, global parsimony, global psm fdr
			if (getParameters()->getCurrentRawFileList().size() > 1 || getParameters()->getSearchParameters()->getWriteMzId() || getParameters()->getSearchParameters()->getWritePepXml())
			{
				FileSystem::createDirectory(getParameters()->getIndividualResultsOutputFolder());

				for (auto fullFilePath : getPsmsGroupedByFile().Select([&] (std::any v)
				{
					v::Key;
				}))
				{
					std::string strippedFileName = Path::GetFileNameWithoutExtension(fullFilePath);

					std::vector<PeptideSpectralMatch*> psmsForThisFile = getPsmsGroupedByFile().Where([&] (std::any p)
					{
						return p->Key == fullFilePath;
					}).SelectMany([&] (std::any g)
					{
						return g;
					}).ToList();
					auto subsetProteinGroupsForThisFile = getProteinGroups().Select([&] (std::any p)
					{
						p::ConstructSubsetProteinGroup(fullFilePath);
					}).ToList();

					ProteinScoringAndFdrEngine tempVar(subsetProteinGroupsForThisFile, psmsForThisFile, getParameters()->getSearchParameters()->getNoOneHitWonders(), getParameters()->getSearchParameters()->getModPeptidesAreDifferent(), false, getCommonParameters(), new std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", fullFilePath});
					ProteinScoringAndFdrResults *subsetProteinScoringAndFdrResults = static_cast<ProteinScoringAndFdrResults*>((&tempVar)->Run());

					subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults->SortedAndScoredProteinGroups;

					getParameters()->getSearchTaskResults()->AddNiceText("Target protein groups within 1 % FDR in " + strippedFileName + ": " + subsetProteinGroupsForThisFile.size()([&] (std::any b)
					{
						return b::QValue <= 0.01 && !b::IsDecoy;
					}));

					// write individual spectra file protein groups results to tsv
					if (getParameters()->getCurrentRawFileList().size() > 1)
					{
						writtenFile = FileSystem::combine(getParameters()->getIndividualResultsOutputFolder(), strippedFileName + "_ProteinGroups.tsv");
						WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, writtenFile, {getParameters()->getSearchTaskId(), "Individual Spectra Files", fullFilePath}, getCommonParameters()->getQValueOutputFilter());
					}

					// write mzID
					if (getParameters()->getSearchParameters()->getWriteMzId())
					{
						Status("Writing mzID...", std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", fullFilePath});

						auto mzidFilePath = FileSystem::combine(getParameters()->getIndividualResultsOutputFolder(), strippedFileName + ".mzID");
						MzIdentMLWriter::WriteMzIdentMl(psmsForThisFile, subsetProteinGroupsForThisFile, getParameters()->getVariableModifications(), getParameters()->getFixedModifications(), {getCommonParameters()->getDigestionParams()->Protease}, getCommonParameters()->getQValueOutputFilter(), getCommonParameters()->getProductMassTolerance(), getCommonParameters()->getPrecursorMassTolerance(), getCommonParameters()->getDigestionParams()->MaxMissedCleavages, mzidFilePath);

						FinishedWritingFile(mzidFilePath, std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", fullFilePath});
					}

					// write pepXML
					if (getParameters()->getSearchParameters()->getWritePepXml())
					{
						Status("Writing pepXML...", std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", fullFilePath});

						auto pepXMLFilePath = FileSystem::combine(getParameters()->getIndividualResultsOutputFolder(), strippedFileName + ".pep.XM");
						PepXMLWriter::WritePepXml(psmsForThisFile, getParameters()->getDatabaseFilenameList(), getParameters()->getVariableModifications(), getParameters()->getFixedModifications(), getCommonParameters(), pepXMLFilePath, getCommonParameters()->getQValueOutputFilter());

						FinishedWritingFile(pepXMLFilePath, std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", fullFilePath});
					}

					ProgressEventArgs tempVar2(100, "Done!", new std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", fullFilePath});
					ReportProgress(&tempVar2);
				}
			}
		}
	}

	void PostSearchAnalysisTask::WriteQuantificationResults()
	{
		if (getParameters()->getSearchParameters()->getDoQuantification() && getParameters()->getFlashLfqResults() != nullptr)
		{
			// write peaks
			WritePeakQuantificationResultsToTsv(getParameters()->getFlashLfqResults(), getParameters()->getOutputFolder(), "AllQuantifiedPeaks", std::vector<std::string> {getParameters()->getSearchTaskId()});

			// write peptide quant results
			WritePeptideQuantificationResultsToTsv(getParameters()->getFlashLfqResults(), getParameters()->getOutputFolder(), "AllQuantifiedPeptides", std::vector<std::string> {getParameters()->getSearchTaskId()});

			// write individual results
			if (getParameters()->getCurrentRawFileList().size() > 1)
			{
				for (auto file : getParameters()->getFlashLfqResults()->Peaks)
				{
					WritePeakQuantificationResultsToTsv(getParameters()->getFlashLfqResults(), getParameters()->getIndividualResultsOutputFolder(), file->Key->FilenameWithoutExtension + "_QuantifiedPeaks", std::vector<std::string> {getParameters()->getSearchTaskId(), "Individual Spectra Files", file->Key->FullFilePathWithExtension});
				}
			}
		}
	}

	void PostSearchAnalysisTask::WritePrunedDatabase()
	{
		if (getParameters()->getSearchParameters()->getWritePrunedDatabase())
		{
			Status("Writing Pruned Database...", std::vector<std::string> {getParameters()->getSearchTaskId()});
			std::unordered_set<Modification*> modificationsToWriteIfBoth;
			std::unordered_set<Modification*> modificationsToWriteIfInDatabase;
			std::unordered_set<Modification*> modificationsToWriteIfObserved;

			auto confidentPsms = getParameters()->getAllPsms().Where([&] (std::any b)
			{
				return b::FdrInfo::QValueNotch <= 0.01 && b::FdrInfo::QValue <= 0.01 && !b::IsDecoy && b::BaseSequence != nullptr;
			}).ToList();
			auto proteinToConfidentBaseSequences = std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>();

			// associate all confident PSMs with all possible proteins they could be digest products of (before or after parsimony)
			for (auto psm : confidentPsms)
			{
				auto myPepsWithSetMods = psm->BestMatchingPeptides->Select([&] (std::any p)
				{
					p::Peptide;
				});

				for (auto peptide : myPepsWithSetMods)
				{
					TValue myPepList;
					std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>::const_iterator proteinToConfidentBaseSequences_iterator = proteinToConfidentBaseSequences.find(peptide.Protein.NonVariantProtein);
					if (proteinToConfidentBaseSequences_iterator != proteinToConfidentBaseSequences.end())
					{
						myPepList = proteinToConfidentBaseSequences_iterator->second;
						myPepList->Add(peptide);
					}
					else
					{
						myPepList = proteinToConfidentBaseSequences_iterator->second;
						proteinToConfidentBaseSequences.emplace(peptide->Protein.NonVariantProtein, std::vector<PeptideWithSetModifications*> {peptide});
					}
				}
			}

			// Add user mod selection behavours to Pruned DB
			for (auto modType : getParameters()->getSearchParameters()->getModsToWriteSelection())
			{
				for (Modification *mod : GlobalVariables::getAllModsKnown().Where([&] (std::any b)
				{
					b::ModificationType->Equals(modType.Key);
				}))
				{
					if (modType.Value == 1) // Write if observed and in database
					{
						modificationsToWriteIfBoth.insert(mod);
					}
					if (modType.Value == 2) // Write if in database
					{
						modificationsToWriteIfInDatabase.insert(mod);
					}
					if (modType.Value == 3) // Write if observed
					{
						modificationsToWriteIfObserved.insert(mod);
					}
				}
			}

			//generates dictionary of proteins with only localized modifications
			auto ModPsms = getParameters()->getAllPsms().Where([&] (std::any b)
			{
				return b::FdrInfo::QValueNotch <= 0.01 && b::FdrInfo::QValue <= 0.01 && !b::IsDecoy && b::FullSequence != nullptr;
			}).ToList();
			auto proteinToConfidentModifiedSequences = std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>();

			for (auto psm : ModPsms)
			{
				auto myPepsWithSetMods = psm->BestMatchingPeptides->Select([&] (std::any p)
				{
					p::Peptide;
				});

				for (auto peptide : myPepsWithSetMods)
				{
					TValue myPepList;
					std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>::const_iterator proteinToConfidentModifiedSequences_iterator = proteinToConfidentModifiedSequences.find(peptide.Protein.NonVariantProtein);
					if (proteinToConfidentModifiedSequences_iterator != proteinToConfidentModifiedSequences.end())
					{
						myPepList = proteinToConfidentModifiedSequences_iterator->second;
						myPepList->Add(peptide);
					}
					else
					{
						myPepList = proteinToConfidentModifiedSequences_iterator->second;
						proteinToConfidentModifiedSequences.emplace(peptide->Protein.NonVariantProtein, std::vector<PeptideWithSetModifications*> {peptide});
					}
				}
			}

			// mods included in pruned database will only be confidently localized mods (peptide's FullSequence != null)
			for (auto nonVariantProtein : getParameters()->getProteinList().Select([&] (std::any p)
			{
				p::NonVariantProtein;
			}).Distinct())
			{
				if (!nonVariantProtein::IsDecoy)
				{
					TValue psms;
					std::unordered_map<Protein*, std::vector<PeptideWithSetModifications*>>::const_iterator proteinToConfidentModifiedSequences_iterator = proteinToConfidentModifiedSequences.find(nonVariantProtein);
					psms = proteinToConfidentModifiedSequences_iterator->second;
					std::unordered_set<(int, Modification, SequenceVariation)*> modsObservedOnThisProtein; // sequence variant is null if mod is not on a variant
					for (auto psm : (psms != nullptr) ? psms : std::vector<PeptideWithSetModifications*>())
					{
						for (auto idxModKV : psm->AllModsOneIsNterminus)
						{
							int proteinIdx = GetOneBasedIndexInProtein(idxModKV->Key, psm);
							SequenceVariation *relevantVariant = psm->Protein.AppliedSequenceVariations.FirstOrDefault([&] (std::any sv)
							{
								VariantApplication::IsSequenceVariantModification(sv, proteinIdx);
							});
							SequenceVariation *unappliedVariant = relevantVariant == nullptr ? nullptr : psm->Protein.SequenceVariations.FirstOrDefault([&] (std::any sv)
							{
								return sv::Description != nullptr && sv::Description->Equals(relevantVariant->Description);
							});
							modsObservedOnThisProtein.insert((VariantApplication::RestoreModificationIndex(psm->Protein, proteinIdx), idxModKV->Value, unappliedVariant));
						}
					}

					std::unordered_map<(SequenceVariation, int)*, std::vector<Modification*>> modsToWrite;

					//Add if observed (regardless if in database)
					for (auto observedMod : modsObservedOnThisProtein)
					{
						auto tempMod = observedMod->Item2;

						if (std::find(modificationsToWriteIfObserved.begin(), modificationsToWriteIfObserved.end(), tempMod) != modificationsToWriteIfObserved.end())
						{
							auto svIdxKey = (observedMod->Item3, observedMod->Item1);
							if (modsToWrite.find(svIdxKey) == modsToWrite.end())
							{
								modsToWrite.emplace(svIdxKey, std::vector<Modification*> {observedMod->Item2});
							}
							else
							{
								modsToWrite[svIdxKey].push_back(observedMod->Item2);
							}
						}
					}

					// Add modification if in database (two cases: always or if observed)
					for (auto modkv : nonVariantProtein::OneBasedPossibleLocalizedModifications)
					{
						for (auto mod : modkv->Value)
						{
							//Add if always In Database or if was observed and in database and not set to not include
							if (std::find(modificationsToWriteIfInDatabase.begin(), modificationsToWriteIfInDatabase.end(), mod) != modificationsToWriteIfInDatabase.end() || (std::find(modificationsToWriteIfBoth.begin(), modificationsToWriteIfBoth.end(), mod) != modificationsToWriteIfBoth.end() && std::find(modsObservedOnThisProtein.begin(), modsObservedOnThisProtein.end(), (modkv->Key, mod, nullptr)) != modsObservedOnThisProtein.end())))
							{
								if (modsToWrite.find((nullptr, modkv->Key)) == modsToWrite.end())
								{
									modsToWrite.emplace((nullptr, modkv->Key), std::vector<Modification*> {mod});
								}
								else
								{
									modsToWrite[(nullptr, modkv->Key)].push_back(mod);
								}
							}
						}
					}

					// Add variant modification if in database (two cases: always or if observed)
					for (SequenceVariation *sv : nonVariantProtein::SequenceVariations)
					{
						for (auto modkv : sv->OneBasedModifications)
						{
							for (auto mod : modkv->Value)
							{
								//Add if always In Database or if was observed and in database and not set to not include
								if (std::find(modificationsToWriteIfInDatabase.begin(), modificationsToWriteIfInDatabase.end(), mod) != modificationsToWriteIfInDatabase.end() || (std::find(modificationsToWriteIfBoth.begin(), modificationsToWriteIfBoth.end(), mod) != modificationsToWriteIfBoth.end() && std::find(modsObservedOnThisProtein.begin(), modsObservedOnThisProtein.end(), (modkv->Key, mod, sv)) != modsObservedOnThisProtein.end())))
								{
									if (modsToWrite.find((sv, modkv->Key)) == modsToWrite.end())
									{
										modsToWrite.emplace((sv, modkv->Key), std::vector<Modification*> {mod});
									}
									else
									{
										modsToWrite[(sv, modkv->Key)].push_back(mod);
									}
								}
							}
						}
					}

					if (proteinToConfidentBaseSequences.find(nonVariantProtein::NonVariantProtein) != proteinToConfidentBaseSequences.end())
					{
						// adds confidently localized and identified mods
						nonVariantProtein::OneBasedPossibleLocalizedModifications->Clear();
						for (auto kvp : modsToWrite.Where([&] (std::any kv)
						{
							return kv::Key->Item1 == nullptr;
						}))
						{
							nonVariantProtein::OneBasedPossibleLocalizedModifications->Add(kvp::Key->Item2, kvp->Value);
						}
						for (auto sv : nonVariantProtein::SequenceVariations)
						{
							sv->OneBasedModifications->Clear();
							for (auto kvp : modsToWrite.Where([&] (std::any kv)
							{
								return kv::Key->Item1 != nullptr && kv::Key->Item1->Equals(sv);
							}))
							{
								sv->OneBasedModifications->Add(kvp::Key->Item2, kvp->Value);
							}
						}
					}
				}
			}

			//writes all proteins
			if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)
			{
				!b::IsContaminant;
			}))
			{
				std::string outputXMLdbFullName = FileSystem::combine(getParameters()->getOutputFolder(), std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)
				{
					!b::IsContaminant;
				})->Select([&] (std::any b)
				{
					Path::GetFileNameWithoutExtension(b::FilePath);
				})) + "pruned.xml");
				ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(), getParameters()->getProteinList().Select([&] (std::any p)
				{
					p::NonVariantProtein;
				}).Where([&] (std::any b)
				{
					return !b::IsDecoy && !b::IsContaminant;
				}).ToList(), outputXMLdbFullName);
				FinishedWritingFile(outputXMLdbFullName, std::vector<std::string> {getParameters()->getSearchTaskId()});
			}
			if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)
			{
				b::IsContaminant;
			}))
			{
				std::string outputXMLdbFullNameContaminants = FileSystem::combine(getParameters()->getOutputFolder(), std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)
				{
					b::IsContaminant;
				})->Select([&] (std::any b)
				{
					Path::GetFileNameWithoutExtension(b::FilePath);
				})) + "pruned.xml");
				ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(), getParameters()->getProteinList().Select([&] (std::any p)
				{
					p::NonVariantProtein;
				}).Where([&] (std::any b)
				{
					return !b::IsDecoy && b::IsContaminant;
				}).ToList(), outputXMLdbFullNameContaminants);
				FinishedWritingFile(outputXMLdbFullNameContaminants, std::vector<std::string> {getParameters()->getSearchTaskId()});
			}

			//writes only detected proteins
			if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)
			{
				!b::IsContaminant;
			}))
			{
				std::string outputXMLdbFullName = FileSystem::combine(getParameters()->getOutputFolder(), std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)
				{
					!b::IsContaminant;
				})->Select([&] (std::any b)
				{
					Path::GetFileNameWithoutExtension(b::FilePath);
				})) + "proteinPruned.xml");
				ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(), proteinToConfidentBaseSequences.Keys->Where([&] (std::any b)
				{
					return !b::IsDecoy && !b::IsContaminant;
				}).ToList(), outputXMLdbFullName);
				FinishedWritingFile(outputXMLdbFullName, std::vector<std::string> {getParameters()->getSearchTaskId()});
			}
			if (getParameters()->getDatabaseFilenameList().Any([&] (std::any b)
			{
				b::IsContaminant;
			}))
			{
				std::string outputXMLdbFullNameContaminants = FileSystem::combine(getParameters()->getOutputFolder(), std::string::Join("-", getParameters()->getDatabaseFilenameList().Where([&] (std::any b)
				{
					b::IsContaminant;
				})->Select([&] (std::any b)
				{
					Path::GetFileNameWithoutExtension(b::FilePath);
				})) + "proteinPruned.xml");
				ProteinDbWriter::WriteXmlDatabase(std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>>(), proteinToConfidentBaseSequences.Keys->Where([&] (std::any b)
				{
					return !b::IsDecoy && b::IsContaminant;
				}).ToList(), outputXMLdbFullNameContaminants);
				FinishedWritingFile(outputXMLdbFullNameContaminants, std::vector<std::string> {getParameters()->getSearchTaskId()});
			}
		}
	}

	int PostSearchAnalysisTask::GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications *peptideWithSetModifications)
	{
		if (oneIsNterminus == 1)
		{
			return peptideWithSetModifications->OneBasedStartResidueInProtein;
		}
		if (oneIsNterminus == peptideWithSetModifications->Length + 2)
		{
			return peptideWithSetModifications->OneBasedEndResidueInProtein;
		}
		return peptideWithSetModifications->OneBasedStartResidueInProtein + oneIsNterminus - 2;
	}

	void PostSearchAnalysisTask::WriteTree(BinTreeStructure *myTreeStructure, const std::string &writtenFile)
	{
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(writtenFile))
		{
			StreamWriter output = StreamWriter(writtenFile);
			output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tUnimodDiffs\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tprotNtermLocFrac\tpepNtermLocFrac\tpepCtermLocFrac\tprotCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot");
			for (Bin *bin : myTreeStructure->getFinalBins().OrderByDescending([&] (std::any b)
			{
				b->Count;
			}))
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				output.WriteLine(bin::MassShift.ToString("F4", CultureInfo::InvariantCulture) + "\t" + bin->Count.ToString(CultureInfo::InvariantCulture) + "\t" + bin::CountDecoy.ToString(CultureInfo::InvariantCulture) + "\t" + bin::CountTarget.ToString(CultureInfo::InvariantCulture) + "\t" + bin::LocalizeableTarget.ToString(CultureInfo::InvariantCulture) + "\t" + (bin::CountTarget - bin::LocalizeableTarget).ToString(CultureInfo::InvariantCulture) + "\t" + (bin->Count == 0 ? NAN : static_cast<double>(bin::CountDecoy) / bin->Count).ToString("F3", CultureInfo::InvariantCulture) + "\t" + (Normal::CDF(0, 1, bin::ComputeZ(0.01))).ToString("F3", CultureInfo::InvariantCulture) + "\t" + (Normal::CDF(0, 1, bin::ComputeZ(0.255))).ToString("F3", CultureInfo::InvariantCulture) + "\t" + (bin->CountTarget == 0 ? NAN : static_cast<double>(bin::LocalizeableTarget) / bin::CountTarget).ToString("F3", CultureInfo::InvariantCulture) + "\t" + bin::Mine + "\t" + bin::UnimodId + "\t" + bin::UnimodFormulas + "\t" + bin::UnimodDiffs + "\t" + bin::AA + "\t" + bin::Combos + "\t" + std::string::Join(",", bin::ModsInCommon::OrderByDescending([&] (std::any b)
				{
					b->Value;
				}).Where([&] (std::any b)
				{
					return b->Value > bin::CountTarget / 10.0;
				})->Select([&] (std::any b)
				{
					return b::Key + ":" + (static_cast<double>(b->Value) / bin::CountTarget).ToString("F3", CultureInfo::InvariantCulture);
				})) + "\t" + std::string::Join(",", bin::AAsInCommon::OrderByDescending([&] (std::any b)
				{
					b->Value;
				}).Where([&] (std::any b)
				{
					return b->Value > bin::CountTarget / 10.0;
				})->Select([&] (std::any b)
				{
					return b::Key + ":" + (static_cast<double>(b->Value) / bin::CountTarget).ToString("F3", CultureInfo::InvariantCulture);
				})) + "\t" + std::string::Join(",", bin::ResidueCount::OrderByDescending([&] (std::any b)
				{
					b->Value;
				})->Select([&] (std::any b)
				{
					return b::Key + ":" + b->Value;
				})) + "\t" + (bin->LocalizeableTarget == 0 ? NAN : static_cast<double>(bin::ProtNlocCount) / bin::LocalizeableTarget).ToString("F3", CultureInfo::InvariantCulture) + "\t" + (bin->LocalizeableTarget == 0 ? NAN : static_cast<double>(bin::PepNlocCount) / bin::LocalizeableTarget).ToString("F3", CultureInfo::InvariantCulture) + "\t" + (bin->LocalizeableTarget == 0 ? NAN : static_cast<double>(bin::PepClocCount) / bin::LocalizeableTarget).ToString("F3", CultureInfo::InvariantCulture) + "\t" + (bin->LocalizeableTarget == 0 ? NAN : static_cast<double>(bin::ProtClocCount) / bin::LocalizeableTarget).ToString("F3", CultureInfo::InvariantCulture) + "\t" + (bin::FracWithSingle).ToString("F3", CultureInfo::InvariantCulture) + "\t" + (static_cast<double>(bin::Overlapping) / bin::CountTarget).ToString("F3", CultureInfo::InvariantCulture) + "\t" + (bin::MedianLength).ToString("F3", CultureInfo::InvariantCulture) + "\t" + bin::UniprotID);
			}
		}
	}

	void PostSearchAnalysisTask::WritePsmsForPercolator(std::vector<PeptideSpectralMatch*> &psmList, const std::string &writtenFileForPercolator, double qValueCutoff)
	{
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(writtenFileForPercolator))
		{
			StreamWriter output = StreamWriter(writtenFileForPercolator);
			output.WriteLine("SpecId\tLabel\tScanNr\tF1\tF2\tPeptide\tProteins");
			output.WriteLine("DefaultDirection\t-\t-\t1\t1\t\t");
			for (int i = 0; i < psmList.size(); i++)
			{
				auto psm = psmList[i];

				if (psm->getFdrInfo()->getQValue() > qValueCutoff || psm->getFdrInfo()->getQValueNotch() > qValueCutoff)
				{
					continue;
				}

				output.Write(std::to_string(i));
				output.Write(L'\t' + std::to_string(psm->getIsDecoy() ? -1 : 1));
				output.Write(L'\t' + std::to_string(psm->getScanNumber()));

				// Features
				output.Write(StringHelper::toString(L'\t') + std::string::Join("\t", psm->getFeatures()));

				// HACKY: Ignores all ambiguity
				auto pwsm = psm->BestMatchingPeptides.First().Peptide;

//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				output.Write(L'\t' + (pwsm->PreviousAminoAcid + "." + pwsm->FullSequence + "." + pwsm->NextAminoAcid).ToString());
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				output.Write(L'\t' + (pwsm->Protein.Accession).ToString());
				output.WriteLine();
			}
		}
	}

	void PostSearchAnalysisTask::WriteProteinGroupsToTsv(std::vector<EngineLayer::ProteinGroup*> &proteinGroups, const std::string &filePath, std::vector<std::string> &nestedIds, double qValueCutoff)
	{
		if (proteinGroups.size() > 0 && proteinGroups.Any())
		{
//C# TO C++ CONVERTER NOTE: The following 'using' block is replaced by its C++ equivalent:
//ORIGINAL LINE: using (StreamWriter output = new StreamWriter(filePath))
			{
				StreamWriter output = StreamWriter(filePath);
				output.WriteLine(proteinGroups.front().GetTabSeparatedHeader());
				for (int i = 0; i < proteinGroups.size(); i++)
				{
					if ((!getParameters()->getSearchParameters()->getWriteDecoys() && proteinGroups[i]->getIsDecoy()) || (!getParameters()->getSearchParameters()->getWriteContaminants() && proteinGroups[i]->getIsContaminant()))
					{
						continue;
					}
					else if (proteinGroups[i]->getQValue() <= qValueCutoff)
					{
						output.WriteLine(proteinGroups[i]);
					}
				}
			}

			FinishedWritingFile(filePath, nestedIds);
		}
	}

	void PostSearchAnalysisTask::WritePeptideQuantificationResultsToTsv(FlashLfqResults *flashLFQResults, const std::string &outputFolder, const std::string &fileName, std::vector<std::string> &nestedIds)
	{
		auto fullSeqPath = FileSystem::combine(outputFolder, fileName + ".tsv");

		flashLFQResults->WriteResults(nullptr, fullSeqPath, nullptr);

		FinishedWritingFile(fullSeqPath, nestedIds);
	}

       void PostSearchAnalysisTask::WritePeakQuantificationResultsToTsv(FlashLfqResults *flashLFQResults,
                                                                         const std::string &outputFolder,
                                                                         const std::string &fileName,
                                                                         std::vector<std::string> &nestedIds)
	{
            auto peaksPath = FileSystem::combine(outputFolder, fileName + ".tsv");
            
            flashLFQResults->WriteResults(peaksPath, nullptr, nullptr);
            
            FinishedWritingFile(peaksPath, nestedIds);
	}
}
