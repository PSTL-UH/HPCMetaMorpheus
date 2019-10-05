#pragma once

#include "TomlReadWrite.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Nett;
using namespace Proteomics::ProteolyticDigestion;
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
#ifdef ORIG
		FileSpecificParameters(TomlTable *tomlTable);
#endif
		//New constructor using tinytoml
		FileSpecificParameters(toml::Table tomlTable);
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

		// This method is to make sure developers keep consistent naming between CommonParameters and FileSpecificParameters.
		// It's supposed to immediately crash MetaMorpheus if you rename a Common Parameter and don't rename it here.
		// The reason this method exists is to make sure toml settings are written and parsed consistently between the tasks
		// and the file-specific settings.
		static void ValidateFileSpecificVariableNames();

		FileSpecificParameters *Clone();
	};
}
