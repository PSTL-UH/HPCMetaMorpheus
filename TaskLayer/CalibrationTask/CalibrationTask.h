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
namespace TaskLayer { class CalibrationParameters; }
namespace TaskLayer { class DbForTask; }
namespace TaskLayer { class FileSpecificParameters; }
namespace TaskLayer { class MyTaskResults; }
namespace EngineLayer { class CommonParameters; }

using namespace EngineLayer;
using namespace EngineLayer::Calibration;
using namespace EngineLayer::ClassicSearch;
using namespace EngineLayer::FdrAnalysis;
using namespace IO::MzML;
using namespace MassSpectrometry;
using namespace MzLibUtil;
using namespace Nett;
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
		MyTaskResults *RunSpecific(const std::wstring &OutputFolder, std::vector<DbForTask*> &dbFilenameList, std::vector<std::wstring> &currentRawFileList, const std::wstring &taskId, std::vector<FileSpecificParameters*> &fileSettingsList) override;

	private:
		int NumRequiredPsms = 20;
		int NumRequiredMs1Datapoints = 50;
		int NumRequiredMs2Datapoints = 100;
	public:
		static const std::wstring CalibSuffix;

	private:
		bool ImprovGlobal(double prevPrecTol, double prevProdTol, int prevPsmCount, int thisRoundPsmCount, double thisRoundPrecTol, double thisRoundProdTol);

		DataPointAquisitionResults *GetDataAcquisitionResults(MsDataFile *myMsDataFile, const std::wstring &currentDataFile, std::vector<Modification*> &variableModifications, std::vector<Modification*> &fixedModifications, std::vector<Protein*> &proteinList, const std::wstring &taskId, EngineLayer::CommonParameters *combinedParameters, Tolerance *initPrecTol, Tolerance *initProdTol);

		static void WriteNewExperimentalDesignFile(const std::wstring &assumedPathToExperDesign, const std::wstring &outputFolder, std::vector<std::wstring> &spectraFilesAfterCalibration);
	};
}
