#include "FileSpecificParameters.h"
#include "../EngineLayer/MetaMorpheusException.h"
#include "../EngineLayer/CommonParameters.h"

using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Nett;
using namespace Proteomics::ProteolyticDigestion;
using namespace MassSpectrometry;

namespace TaskLayer
{

	FileSpecificParameters::FileSpecificParameters(TomlTable *tomlTable)
	{
		for (auto keyValuePair : tomlTable)
		{
			switch (keyValuePair->Key)
			{
				// we're using the name of the variable here and not a fixed string
				// in case the variable name changes at some point
				case "PrecursorMassTolerance":
					setPrecursorMassTolerance(keyValuePair->Value->Get<Tolerance*>());
					break;
				case "ProductMassTolerance":
					setProductMassTolerance(keyValuePair->Value->Get<Tolerance*>());
					break;
				case "Protease":
					setProtease(keyValuePair->Value->Get<getProtease()*>());
					break;
				case "MinPeptideLength":
					setMinPeptideLength(keyValuePair->Value->Get<int>());
					break;
				case "MaxPeptideLength":
					setMaxPeptideLength(keyValuePair->Value->Get<int>());
					break;
				case "MaxMissedCleavages":
					setMaxMissedCleavages(keyValuePair->Value->Get<int>());
					break;
				case "MaxModsForPeptide":
					setMaxModsForPeptide(keyValuePair->Value->Get<int>());
					break;


				case "DissociationType":
					setDissociationType(keyValuePair->Value->Get<MassSpectrometry::DissociationType*>());
					break;


				default:
					throw MetaMorpheusException("Unrecognized parameter \"" + keyValuePair->Key + "\" in file-specific parameters toml");
			}
		}
	}

	FileSpecificParameters::FileSpecificParameters()
	{
		// everything initialized to null
	}

	Tolerance *FileSpecificParameters::getPrecursorMassTolerance() const
	{
		return privatePrecursorMassTolerance;
	}

	void FileSpecificParameters::setPrecursorMassTolerance(Tolerance *value)
	{
		privatePrecursorMassTolerance = value;
	}

	Tolerance *FileSpecificParameters::getProductMassTolerance() const
	{
		return privateProductMassTolerance;
	}

	void FileSpecificParameters::setProductMassTolerance(Tolerance *value)
	{
		privateProductMassTolerance = value;
	}

	Protease *FileSpecificParameters::getProtease() const
	{
		return privateProtease;
	}

	void FileSpecificParameters::setProtease(Protease *value)
	{
		privateProtease = value;
	}

	std::optional<int> FileSpecificParameters::getMinPeptideLength() const
	{
		return privateMinPeptideLength;
	}

	void FileSpecificParameters::setMinPeptideLength(std::optional<int> value)
	{
		privateMinPeptideLength = value;
	}

	std::optional<int> FileSpecificParameters::getMaxPeptideLength() const
	{
		return privateMaxPeptideLength;
	}

	void FileSpecificParameters::setMaxPeptideLength(std::optional<int> value)
	{
		privateMaxPeptideLength = value;
	}

	std::optional<int> FileSpecificParameters::getMaxMissedCleavages() const
	{
		return privateMaxMissedCleavages;
	}

	void FileSpecificParameters::setMaxMissedCleavages(std::optional<int> value)
	{
		privateMaxMissedCleavages = value;
	}

	std::optional<int> FileSpecificParameters::getMaxModsForPeptide() const
	{
		return privateMaxModsForPeptide;
	}

	void FileSpecificParameters::setMaxModsForPeptide(std::optional<int> value)
	{
		privateMaxModsForPeptide = value;
	}

	std::optional<DissociationType*> FileSpecificParameters::getDissociationType() const
	{
		return privateDissociationType;
	}

	void FileSpecificParameters::setDissociationType(std::optional<DissociationType*> value)
	{
		privateDissociationType = value;
	}

	void FileSpecificParameters::ValidateFileSpecificVariableNames()
	{
		CommonParameters *temp = new CommonParameters();

		if ("PrecursorMassTolerance" != "PrecursorMassTolerance")
		{
			delete temp;
			throw MetaMorpheusException("Precursor tol variable name is inconsistent");
		}
		if ("ProductMassTolerance" != "ProductMassTolerance")
		{
			delete temp;
			throw MetaMorpheusException("Product tol variable name is inconsistent");
		}
		if ("Protease" != "Protease")
		{
			delete temp;
			throw MetaMorpheusException("Protease variable name is inconsistent");
		}
		if ("MinPeptideLength" != "MinPeptideLength")
		{
			delete temp;
			throw MetaMorpheusException("Min peptide length variable name is inconsistent");
		}
		if ("MaxPeptideLength" != "MaxPeptideLength")
		{
			delete temp;
			throw MetaMorpheusException("Max peptide length variable name is inconsistent");
		}
		if ("MaxMissedCleavages" != "MaxMissedCleavages")
		{
			delete temp;
			throw MetaMorpheusException("Max missed cleavages variable name is inconsistent");
		}
		if ("MaxModsForPeptide" != "MaxModsForPeptide")
		{
			delete temp;
			throw MetaMorpheusException("Max mods per peptide variable name is inconsistent");
		}

		delete temp;
	}

	FileSpecificParameters *FileSpecificParameters::Clone()
	{
		return static_cast<FileSpecificParameters*>(this->MemberwiseClone());
	}
}
