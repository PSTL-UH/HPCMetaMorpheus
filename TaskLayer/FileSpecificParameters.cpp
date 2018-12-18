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
				case L"PrecursorMassTolerance":
					setPrecursorMassTolerance(keyValuePair->Value->Get<Tolerance*>());
					break;
				case L"ProductMassTolerance":
					setProductMassTolerance(keyValuePair->Value->Get<Tolerance*>());
					break;
				case L"Protease":
					setProtease(keyValuePair->Value->Get<getProtease()*>());
					break;
				case L"MinPeptideLength":
					setMinPeptideLength(keyValuePair->Value->Get<int>());
					break;
				case L"MaxPeptideLength":
					setMaxPeptideLength(keyValuePair->Value->Get<int>());
					break;
				case L"MaxMissedCleavages":
					setMaxMissedCleavages(keyValuePair->Value->Get<int>());
					break;
				case L"MaxModsForPeptide":
					setMaxModsForPeptide(keyValuePair->Value->Get<int>());
					break;


				case L"DissociationType":
					setDissociationType(keyValuePair->Value->Get<MassSpectrometry::DissociationType*>());
					break;


				default:
					throw MetaMorpheusException(L"Unrecognized parameter \"" + keyValuePair->Key + L"\" in file-specific parameters toml");
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

	Nullable<int> FileSpecificParameters::getMinPeptideLength() const
	{
		return privateMinPeptideLength;
	}

	void FileSpecificParameters::setMinPeptideLength(Nullable<int> value)
	{
		privateMinPeptideLength = value;
	}

	Nullable<int> FileSpecificParameters::getMaxPeptideLength() const
	{
		return privateMaxPeptideLength;
	}

	void FileSpecificParameters::setMaxPeptideLength(Nullable<int> value)
	{
		privateMaxPeptideLength = value;
	}

	Nullable<int> FileSpecificParameters::getMaxMissedCleavages() const
	{
		return privateMaxMissedCleavages;
	}

	void FileSpecificParameters::setMaxMissedCleavages(Nullable<int> value)
	{
		privateMaxMissedCleavages = value;
	}

	Nullable<int> FileSpecificParameters::getMaxModsForPeptide() const
	{
		return privateMaxModsForPeptide;
	}

	void FileSpecificParameters::setMaxModsForPeptide(Nullable<int> value)
	{
		privateMaxModsForPeptide = value;
	}

	Nullable<DissociationType*> FileSpecificParameters::getDissociationType() const
	{
		return privateDissociationType;
	}

	void FileSpecificParameters::setDissociationType(Nullable<DissociationType*> value)
	{
		privateDissociationType = value;
	}

	void FileSpecificParameters::ValidateFileSpecificVariableNames()
	{
		CommonParameters *temp = new CommonParameters();

		if (L"PrecursorMassTolerance" != L"PrecursorMassTolerance")
		{
			delete temp;
			throw MetaMorpheusException(L"Precursor tol variable name is inconsistent");
		}
		if (L"ProductMassTolerance" != L"ProductMassTolerance")
		{
			delete temp;
			throw MetaMorpheusException(L"Product tol variable name is inconsistent");
		}
		if (L"Protease" != L"Protease")
		{
			delete temp;
			throw MetaMorpheusException(L"Protease variable name is inconsistent");
		}
		if (L"MinPeptideLength" != L"MinPeptideLength")
		{
			delete temp;
			throw MetaMorpheusException(L"Min peptide length variable name is inconsistent");
		}
		if (L"MaxPeptideLength" != L"MaxPeptideLength")
		{
			delete temp;
			throw MetaMorpheusException(L"Max peptide length variable name is inconsistent");
		}
		if (L"MaxMissedCleavages" != L"MaxMissedCleavages")
		{
			delete temp;
			throw MetaMorpheusException(L"Max missed cleavages variable name is inconsistent");
		}
		if (L"MaxModsForPeptide" != L"MaxModsForPeptide")
		{
			delete temp;
			throw MetaMorpheusException(L"Max mods per peptide variable name is inconsistent");
		}

		delete temp;
	}

	FileSpecificParameters *FileSpecificParameters::Clone()
	{
		return static_cast<FileSpecificParameters*>(this->MemberwiseClone());
	}
}
