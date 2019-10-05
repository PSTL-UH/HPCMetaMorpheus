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

#ifdef ORIG
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
#endif

//Dr. Gabriel,  I got a small version of this working without using the different Set functions like setPrecursorMassTolerance.
//It looks to me like they are reading the file specific parameters from a toml table (that has already been read), and then setting various values
//based on the parameters in the table.  I changed how the table is passed to it so it is in line with the tinytoml types we are using
//and I changed the switch statement to if/else statements because the toml::Table keys are strings.  Im not sure if this is the best 
//way to handle this.  I also have not yet found where the fileSpecificParameters are being read.
	FileSpecificParameters::FileSpecificParameters(toml::Table tomlTable)
	{

		for (auto const& keyValuePair : tomlTable)
		{
			//From Nick:  Im not sure if the changes I made are in line with the following comment
			//from the original constructor above:
			//
			// we're using the name of the variable here and not a fixed string
			// in case the variable name changes at some point
			if (keyValuePair.first == "PrecursorMassTolerance") {
				
				//Im not quite sure yet how to deal with the ->Get<Tolerance*>() parts of these
				//of the functions calls.  The keyValuePair.second is the value associated with the 
				//key found in the if else statements.  Here the PrecursorMassTolerance value has already been found
				//Im not sure why we need to perform a Get() operation.  All of the if else statements have something 
				//similar and Im not sure if it is an operation we need to perform or not.
				setPrecursorMassTolerance(keyValuePair.second->Get<Tolerance*>());
				break;
			}
			else if (keyValuePair.first == "ProductMassTolerance") {
				setProductMassTolerance(keyValuePair.second->Get<Tolerance*>());
				break;
			}
			else if (keyValuePair.first == "Protease") {
				setProtease(keyValuePair.second->Get<getProtease()*>());
				break;
			}
			else if (keyValuePair.first == "MinPeptideLength") {
				setMinPeptideLength(keyValuePair.second->Get<int>());
				break;
			}
			else if (keyValuePair.first == "MaxPeptideLength") {
				setMaxPeptideLength(keyValuePair.second->Get<int>());
				break;
			}
			else if (keyValuePair.first == "MaxMissedCleavages") {
				setMaxMissedCleavages(keyValuePair.second->Get<int>());
				break;
			}
			else if (keyValuePair.first == "MaxModsForPeptide") {
				setMaxModsForPeptide(keyValuePair.second->Get<int>());
				break;
			}
			else if (keyValuePair.first == "DissociationType") {
				setDissociationType(keyValuePair.second->Get<MassSpectrometry::DissociationType*>());
				break;
			}

			//Should this be an Else statement?  Not sure how to handle exception
			// default:
			// 	throw MetaMorpheusException("Unrecognized parameter \"" + keyValuePair->Key + "\" in file-specific parameters toml");
			else {
				std::cout << "Unrecognized parameter " << keyValuePair.first << " in file-specific parameters toml" << std::endl;
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
