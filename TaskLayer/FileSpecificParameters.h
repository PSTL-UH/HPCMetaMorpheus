#pragma once

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
		Nullable<int> privateMinPeptideLength;
		Nullable<int> privateMaxPeptideLength;
		Nullable<int> privateMaxMissedCleavages;
		Nullable<int> privateMaxModsForPeptide;
		Nullable<DissociationType*> privateDissociationType;

	public:
		FileSpecificParameters(TomlTable *tomlTable);

		FileSpecificParameters();

		Tolerance *getPrecursorMassTolerance() const;
		void setPrecursorMassTolerance(Tolerance *value);
		Tolerance *getProductMassTolerance() const;
		void setProductMassTolerance(Tolerance *value);
		Protease *getProtease() const;
		void setProtease(Protease *value);
		Nullable<int> getMinPeptideLength() const;
		void setMinPeptideLength(Nullable<int> value);
		Nullable<int> getMaxPeptideLength() const;
		void setMaxPeptideLength(Nullable<int> value);
		Nullable<int> getMaxMissedCleavages() const;
		void setMaxMissedCleavages(Nullable<int> value);
		Nullable<int> getMaxModsForPeptide() const;
		void setMaxModsForPeptide(Nullable<int> value);
		Nullable<DissociationType*> getDissociationType() const;
		void setDissociationType(Nullable<DissociationType*> value);

		// This method is to make sure developers keep consistent naming between CommonParameters and FileSpecificParameters.
		// It's supposed to immediately crash MetaMorpheus if you rename a Common Parameter and don't rename it here.
		// The reason this method exists is to make sure toml settings are written and parsed consistently between the tasks
		// and the file-specific settings.
		static void ValidateFileSpecificVariableNames();

		FileSpecificParameters *Clone();
	};
}
