#pragma once

#include "../MetaMorpheusTask.h"
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include "stringhelper.h"
#include "stringbuilder.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
//namespace TaskLayer { class CalibrationParameters; }
#include "CalibrationParameters.h"

//namespace TaskLayer { class DbForTask; }
#include "../DbForTask.h"

namespace TaskLayer { class FileSpecificParameters; }
namespace TaskLayer { class MyTaskResults; }

//namespace EngineLayer { class CommonParameters; }
#include "../EngineLayer/CommonParameters.h"
using namespace EngineLayer;

//using namespace EngineLayer::Calibration;
//using namespace EngineLayer::ClassicSearch;

#include "../EngineLayer/FdrAnalysis/FdrAnalysisEngine.h"
using namespace EngineLayer::FdrAnalysis;
//using namespace IO::MzML;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "MzLibUtil.h"
using namespace MzLibUtil;

//using namespace Nett;
#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace UsefulProteomicsDatabases;

namespace TaskLayer
{
	class CalibrationTask : public MetaMorpheusTask
	{
	private:
		TaskLayer::CalibrationParameters *privateCalibrationParameters;

	public:
		CalibrationTask();

		TaskLayer::CalibrationParameters *getCalibrationParameters() const;
		void setCalibrationParameters(TaskLayer::CalibrationParameters *value);

	protected:
		MyTaskResults *RunSpecific(const std::string &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::string> &currentRawFileList, const std::string &taskId, std::vector<FileSpecificParameters*> &fileSettingsList) override;

	private:
		int NumRequiredPsms = 20;
		int NumRequiredMs1Datapoints = 50;
		int NumRequiredMs2Datapoints = 100;
	public:
		static const std::string CalibSuffix;

	private:
		bool ImprovGlobal(double prevPrecTol, double prevProdTol, int prevPsmCount, int thisRoundPsmCount, double thisRoundPrecTol, double thisRoundProdTol);

		DataPointAquisitionResults *GetDataAcquisitionResults(MsDataFile *myMsDataFile, const std::string &currentDataFile, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<Protein*> &proteinList, const std::string &taskId, EngineLayer::CommonParameters *combinedParameters, Tolerance *initPrecTol, Tolerance *initProdTol);

		static void WriteNewExperimentalDesignFile(const std::string &assumedPathToExperDesign, const std::string &outputFolder, std::vector<std::string> &spectraFilesAfterCalibration);
	};
}
