#pragma once

#include <string>
#include <vector>
#include <cmath>
#include "exceptionhelper.h"

using namespace Chemistry;
using namespace MassSpectrometry;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace Test
{
	class TestDataFile : public MsDataFile
	{
	public:
		TestDataFile();

		TestDataFile(double closeMassDifference);

		TestDataFile(bool emptyScan);

		TestDataFile(const std::wstring &slightlyLargerDataFile);

		TestDataFile(double precursor, std::vector<double> &products);

		TestDataFile(std::vector<PeptideWithSetModifications*> &pepWithSetModss, bool additionalMasses = false);

		TestDataFile(PeptideWithSetModifications *pepWithSetMods);

		TestDataFile(PeptideWithSetModifications *pepWithSetMods, const std::wstring &v);

		TestDataFile(PeptideWithSetModifications *pepWithSetMods, int charge, double intensity, double rt);

		TestDataFile(int MS3 = 5);

		TestDataFile(std::vector<double> &ms2Mz, std::vector<double> &ms2Intensities, double precursorMass, int precursorZ, double rt = 1.0);

		std::wstring getFilePath() const;

		std::wstring getName() const;

		void ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities);

		MsDataScan *GetOneBasedScan(int scanNumber) override;

		std::vector<MsDataScan*> GetMS1Scans() override;
	};
}
