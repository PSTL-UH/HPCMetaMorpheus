#pragma once

#include "TomlReadWrite.h"

#include "../EngineLayer/CommonParameters.h"
using namespace EngineLayer;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include"Proteomics/ProteolyticDigestion/ProteinDigestion.h"
using namespace Proteomics::ProteolyticDigestion;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

namespace TaskLayer
{
	class FileSpecificParameters
	{
	private:
		Tolerance *privatePrecursorMassTolerance;
		Tolerance *privateProductMassTolerance;
		Protease *privateProtease;
		std::optional<int> privateMinPeptideLength;
		std::optional<int> privateMaxPeptideLength;
		std::optional<int> privateMaxMissedCleavages;
		std::optional<int> privateMaxModsForPeptide;
		std::optional<DissociationType*> privateDissociationType;

	public:
		
		FileSpecificParameters(toml::Table tomlTable);
		FileSpecificParameters(FileSpecificParameters *filep);
		FileSpecificParameters();

		Tolerance *getPrecursorMassTolerance() const;
		void setPrecursorMassTolerance(Tolerance *value);
		Tolerance *getProductMassTolerance() const;
		void setProductMassTolerance(Tolerance *value);
		Protease *getProtease() const;
		void setProtease(Protease *value);
		std::optional<int> getMinPeptideLength() const;
		void setMinPeptideLength(std::optional<int> value);
		std::optional<int> getMaxPeptideLength() const;
		void setMaxPeptideLength(std::optional<int> value);
		std::optional<int> getMaxMissedCleavages() const;
		void setMaxMissedCleavages(std::optional<int> value);
		std::optional<int> getMaxModsForPeptide() const;
		void setMaxModsForPeptide(std::optional<int> value);
		std::optional<DissociationType*> getDissociationType() const;
		void setDissociationType(std::optional<DissociationType*> value);

		toml::Table getFileSpecificParametersTomlTable();

		// This method is to make sure developers keep consistent naming between CommonParameters and FileSpecificParameters.
		// It's supposed to immediately crash MetaMorpheus if you rename a Common Parameter and don't rename it here.
		// The reason this method exists is to make sure toml settings are written and parsed consistently between the tasks
		// and the file-specific settings.
                //
                // EDGAR: not needed for now in the C++ version
		// static void ValidateFileSpecificVariableNames();

                // EDGAR: replaced in the C++ version with a copy-constructor.
		//FileSpecificParameters *Clone();
	};
}
