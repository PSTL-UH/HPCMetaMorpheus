#include "ParameterTest.h"
#include "../EngineLayer/CommonParameters.h"
#include "../TaskLayer/FileSpecificParameters.h"
#include "../TaskLayer/MetaMorpheusTask.h"

using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace NUnit::Framework;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
using namespace TaskLayer;

namespace Test
{

	void ParameterTest::TestFileSpecifcParameterOverwrite()
	{
		CommonParameters *defaultParameters = new CommonParameters();
		AbsoluteTolerance tempVar(69);
		DigestionParams tempVar4(protease: L"Asp-N", maxMissedCleavages: 69, minPeptideLength: 1, maxPeptideLength: 69, maxModificationIsoforms: 69, initiatorMethionineBehavior: InitiatorMethionineBehavior::Retain, maxModsForPeptides: 69, searchModeType: CleavageSpecificity::Semi, fragmentationTerminus: FragmentationTerminus::C);
		AbsoluteTolerance tempVar3(69);
		AbsoluteTolerance tempVar2(69);
		CommonParameters *notDefaultParameters = new CommonParameters(, DissociationType::ETD, false, false, 69, 69, false, true, 69, 69, 69, 69, true, false, true, true, &tempVar, &tempVar2, &tempVar3, 69, &tempVar4, {(L"asdf", L"asdf")}, std::vector<(std::wstring, std::wstring)*> {(L"asdf", L"asdf")});

		//check that the defaults are not the same as the not defaults.
		//IF ONE OF THESE FAILS, PLEASE UPDATE THE "notDefaultParameters"
		Assert::AreNotEqual(defaultParameters->getDissociationType(), notDefaultParameters->getDissociationType());
		Assert::AreNotEqual(defaultParameters->getDoPrecursorDeconvolution(), notDefaultParameters->getDoPrecursorDeconvolution());
		Assert::AreNotEqual(defaultParameters->getUseProvidedPrecursorInfo(), notDefaultParameters->getUseProvidedPrecursorInfo());
		Assert::AreNotEqual(defaultParameters->getDeconvolutionIntensityRatio(), notDefaultParameters->getDeconvolutionIntensityRatio());
		Assert::AreNotEqual(defaultParameters->getDeconvolutionMaxAssumedChargeState(), notDefaultParameters->getDeconvolutionMaxAssumedChargeState());
		Assert::AreNotEqual(defaultParameters->getReportAllAmbiguity(), notDefaultParameters->getReportAllAmbiguity());
		Assert::AreNotEqual(defaultParameters->getAddCompIons(), notDefaultParameters->getAddCompIons());
		Assert::AreNotEqual(defaultParameters->getTotalPartitions(), notDefaultParameters->getTotalPartitions());
		Assert::AreNotEqual(defaultParameters->getScoreCutoff(), notDefaultParameters->getScoreCutoff());
		Assert::AreNotEqual(defaultParameters->getTopNpeaks(), notDefaultParameters->getTopNpeaks());
		Assert::AreNotEqual(defaultParameters->getMinRatio(), notDefaultParameters->getMinRatio());
		Assert::AreNotEqual(defaultParameters->getTrimMs1Peaks(), notDefaultParameters->getTrimMs1Peaks());
		Assert::AreNotEqual(defaultParameters->getTrimMsMsPeaks(), notDefaultParameters->getTrimMsMsPeaks());
		Assert::AreNotEqual(defaultParameters->getUseDeltaScore(), notDefaultParameters->getUseDeltaScore());
		Assert::AreNotEqual(defaultParameters->getCalculateEValue(), notDefaultParameters->getCalculateEValue());
		Assert::AreNotEqual(defaultParameters->getProductMassTolerance(), notDefaultParameters->getProductMassTolerance());
		Assert::AreNotEqual(defaultParameters->getPrecursorMassTolerance(), notDefaultParameters->getPrecursorMassTolerance());
		Assert::AreNotEqual(defaultParameters->getDeconvolutionMassTolerance(), notDefaultParameters->getDeconvolutionMassTolerance());
		Assert::AreNotEqual(defaultParameters->getMaxThreadsToUsePerFile(), notDefaultParameters->getMaxThreadsToUsePerFile());
		Assert::AreNotEqual(defaultParameters->getDigestionParams(), notDefaultParameters->getDigestionParams());
		Assert::AreNotEqual(defaultParameters->ListOfModsVariable, notDefaultParameters->ListOfModsVariable);
		Assert::AreNotEqual(defaultParameters->ListOfModsFixed, notDefaultParameters->ListOfModsFixed);

		FileSpecificParameters *emptyFileSpecificParameters = new FileSpecificParameters();
		CommonParameters *updatedParameters = MetaMorpheusTask::SetAllFileSpecificCommonParams(notDefaultParameters, emptyFileSpecificParameters);

		//CHECK THAT NOTHING CHANGED
		Assert::AreEqual(updatedParameters->getDissociationType(), notDefaultParameters->getDissociationType());
		Assert::AreEqual(updatedParameters->getDoPrecursorDeconvolution(), notDefaultParameters->getDoPrecursorDeconvolution());
		Assert::AreEqual(updatedParameters->getUseProvidedPrecursorInfo(), notDefaultParameters->getUseProvidedPrecursorInfo());
		Assert::AreEqual(updatedParameters->getDeconvolutionIntensityRatio(), notDefaultParameters->getDeconvolutionIntensityRatio());
		Assert::AreEqual(updatedParameters->getDeconvolutionMaxAssumedChargeState(), notDefaultParameters->getDeconvolutionMaxAssumedChargeState());
		Assert::AreEqual(updatedParameters->getReportAllAmbiguity(), notDefaultParameters->getReportAllAmbiguity());
		Assert::AreEqual(updatedParameters->getAddCompIons(), notDefaultParameters->getAddCompIons());
		Assert::AreEqual(updatedParameters->getTotalPartitions(), notDefaultParameters->getTotalPartitions());
		Assert::AreEqual(updatedParameters->getScoreCutoff(), notDefaultParameters->getScoreCutoff());
		Assert::AreEqual(updatedParameters->getTopNpeaks(), notDefaultParameters->getTopNpeaks());
		Assert::AreEqual(updatedParameters->getMinRatio(), notDefaultParameters->getMinRatio());
		Assert::AreEqual(updatedParameters->getTrimMs1Peaks(), notDefaultParameters->getTrimMs1Peaks());
		Assert::AreEqual(updatedParameters->getTrimMsMsPeaks(), notDefaultParameters->getTrimMsMsPeaks());
		Assert::AreEqual(updatedParameters->getUseDeltaScore(), notDefaultParameters->getUseDeltaScore());
		Assert::AreEqual(updatedParameters->getCalculateEValue(), notDefaultParameters->getCalculateEValue());
		Assert::AreEqual(updatedParameters->getProductMassTolerance(), notDefaultParameters->getProductMassTolerance());
		Assert::AreEqual(updatedParameters->getPrecursorMassTolerance(), notDefaultParameters->getPrecursorMassTolerance());
		Assert::AreEqual(updatedParameters->getDeconvolutionMassTolerance(), notDefaultParameters->getDeconvolutionMassTolerance());
		Assert::AreEqual(updatedParameters->getMaxThreadsToUsePerFile(), notDefaultParameters->getMaxThreadsToUsePerFile());
		Assert::AreEqual(updatedParameters->getDigestionParams(), notDefaultParameters->getDigestionParams());
		Assert::AreEqual(updatedParameters->ListOfModsVariable, notDefaultParameters->ListOfModsVariable);
		Assert::AreEqual(updatedParameters->ListOfModsFixed, notDefaultParameters->ListOfModsFixed);

		FileSpecificParameters *basicFileSpecificParameters = new FileSpecificParameters();
		PpmTolerance tempVar5(10);
		basicFileSpecificParameters->setPrecursorMassTolerance(&tempVar5);
		PpmTolerance tempVar6(30);
		basicFileSpecificParameters->setProductMassTolerance(&tempVar6);
		Protease tempVar7(L"Arg-C", CleavageSpecificity::Full, nullptr, nullptr, new std::vector<DigestionMotif*> {new DigestionMotif(L"K", nullptr, 1, L"")});
		basicFileSpecificParameters->setProtease(&tempVar7);
		basicFileSpecificParameters->setMinPeptideLength(std::make_optional(1));
		basicFileSpecificParameters->setMaxPeptideLength(std::make_optional(50));
		basicFileSpecificParameters->setMaxMissedCleavages(std::make_optional(2));
		basicFileSpecificParameters->setMaxModsForPeptide(std::make_optional(1));
		basicFileSpecificParameters->setDissociationType(DissociationType::CID);
		updatedParameters = MetaMorpheusTask::SetAllFileSpecificCommonParams(notDefaultParameters, basicFileSpecificParameters);
		//CHECK THAT SOMETHINGS CHANGED AND OTHERS DIDN'T
		Assert::AreEqual(updatedParameters->getDissociationType(), basicFileSpecificParameters->getDissociationType());
		Assert::AreEqual(updatedParameters->getProductMassTolerance(), basicFileSpecificParameters->getProductMassTolerance());
		Assert::AreEqual(updatedParameters->getPrecursorMassTolerance(), basicFileSpecificParameters->getPrecursorMassTolerance());
		Assert::AreEqual(updatedParameters->getDigestionParams()->MaxModsForPeptide, basicFileSpecificParameters->getMaxModsForPeptide());
		Assert::AreEqual(updatedParameters->getDigestionParams()->MaxMissedCleavages, basicFileSpecificParameters->getMaxMissedCleavages());
		Assert::AreEqual(updatedParameters->getDigestionParams()->MinPeptideLength, basicFileSpecificParameters->getMinPeptideLength());
		Assert::AreEqual(updatedParameters->getDigestionParams()->MaxPeptideLength, basicFileSpecificParameters->getMaxPeptideLength());
		Assert::AreEqual(updatedParameters->getDigestionParams()->Protease, basicFileSpecificParameters->getProtease());

		Assert::AreEqual(updatedParameters->getDoPrecursorDeconvolution(), notDefaultParameters->getDoPrecursorDeconvolution());
		Assert::AreEqual(updatedParameters->getUseProvidedPrecursorInfo(), notDefaultParameters->getUseProvidedPrecursorInfo());
		Assert::AreEqual(updatedParameters->getDeconvolutionIntensityRatio(), notDefaultParameters->getDeconvolutionIntensityRatio());
		Assert::AreEqual(updatedParameters->getDeconvolutionMaxAssumedChargeState(), notDefaultParameters->getDeconvolutionMaxAssumedChargeState());
		Assert::AreEqual(updatedParameters->getReportAllAmbiguity(), notDefaultParameters->getReportAllAmbiguity());
		Assert::AreEqual(updatedParameters->getAddCompIons(), notDefaultParameters->getAddCompIons());
		Assert::AreEqual(updatedParameters->getTotalPartitions(), notDefaultParameters->getTotalPartitions());
		Assert::AreEqual(updatedParameters->getScoreCutoff(), notDefaultParameters->getScoreCutoff());
		Assert::AreEqual(updatedParameters->getTopNpeaks(), notDefaultParameters->getTopNpeaks());
		Assert::AreEqual(updatedParameters->getMinRatio(), notDefaultParameters->getMinRatio());
		Assert::AreEqual(updatedParameters->getTrimMs1Peaks(), notDefaultParameters->getTrimMs1Peaks());
		Assert::AreEqual(updatedParameters->getTrimMsMsPeaks(), notDefaultParameters->getTrimMsMsPeaks());
		Assert::AreEqual(updatedParameters->getUseDeltaScore(), notDefaultParameters->getUseDeltaScore());
		Assert::AreEqual(updatedParameters->getCalculateEValue(), notDefaultParameters->getCalculateEValue());
		Assert::AreEqual(updatedParameters->getDeconvolutionMassTolerance(), notDefaultParameters->getDeconvolutionMassTolerance());
		Assert::AreEqual(updatedParameters->getMaxThreadsToUsePerFile(), notDefaultParameters->getMaxThreadsToUsePerFile());
		Assert::AreEqual(updatedParameters->getDigestionParams()->InitiatorMethionineBehavior, notDefaultParameters->getDigestionParams()->InitiatorMethionineBehavior);
		Assert::AreEqual(updatedParameters->ListOfModsVariable, notDefaultParameters->ListOfModsVariable);
		Assert::AreEqual(updatedParameters->ListOfModsFixed, notDefaultParameters->ListOfModsFixed);

//C# TO C++ CONVERTER TODO TASK: A 'delete basicFileSpecificParameters' statement was not added since basicFileSpecificParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete emptyFileSpecificParameters' statement was not added since emptyFileSpecificParameters was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete notDefaultParameters' statement was not added since notDefaultParameters was passed to a method or constructor. Handle memory management manually.
		delete defaultParameters;
	}
}
