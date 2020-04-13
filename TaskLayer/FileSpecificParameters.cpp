#include "FileSpecificParameters.h"
#include "../EngineLayer/MetaMorpheusException.h"
#include "../EngineLayer/CommonParameters.h"

using namespace EngineLayer;
using namespace MzLibUtil;
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

	FileSpecificParameters::FileSpecificParameters(toml::Table tomlTable)
	{

		for (auto const& keyValuePair : tomlTable)
		{
			// we're using the name of the variable here and not a fixed string
			// in case the variable name changes at some point
			if (keyValuePair.first == "PrecursorMassTolerance") {
				Tolerance* t = Tolerance::ParseToleranceString(keyValuePair.second.as<std::string>());
				setPrecursorMassTolerance(t);
			}
			else if (keyValuePair.first == "ProductMassTolerance") {
				Tolerance* t = Tolerance::ParseToleranceString(keyValuePair.second.as<std::string>());
				setProductMassTolerance(t);
			}
			else if (keyValuePair.first == "Protease") {
				std::string name = keyValuePair.second.as<std::string>();
				std::vector<DigestionMotif*> dm;
				
				//create Protease instance
				//Where do we find these parameters, I havent not found a toml file that specifies any other than "name"?
				//Protease(const std::string &name, Proteomics::ProteolyticDigestion::CleavageSpecificity cleavageSpecificity, const std::string &psiMSAccessionNumber, const std::string &psiMSName, std::vector<DigestionMotif*> &motifList)
				Protease p = Protease(name, Proteomics::ProteolyticDigestion::CleavageSpecificity::Unknown, "", "", dm);
				setProtease(&p);
			}
			else if (keyValuePair.first == "MinPeptideLength") {
				int MinPeptideLength = keyValuePair.second.as<int>();
				setMinPeptideLength(std::make_optional(MinPeptideLength));
			}
			else if (keyValuePair.first == "MaxPeptideLength") {
				int MaxPeptideLength = keyValuePair.second.as<int>();
				setMaxPeptideLength(std::make_optional(MaxPeptideLength));
			}
			else if (keyValuePair.first == "MaxMissedCleavages") {
				int MaxMissedCleavages = keyValuePair.second.as<int>();
				setMaxMissedCleavages(std::make_optional(MaxMissedCleavages));
			}
			else if (keyValuePair.first == "MaxModsForPeptide") {
				int MaxModsForPeptide = keyValuePair.second.as<int>();
				setMaxModsForPeptide(std::make_optional(MaxModsForPeptide));
			}
			else if (keyValuePair.first == "DissociationType") {
				// setDissociationType(keyValuePair.second->Get<MassSpectrometry::DissociationType*>());
				std::string type = keyValuePair.second.as<std::string>();
				MassSpectrometry::DissociationType d = MassSpectrometry::GetDissocationType::GetDissocationTypeFromString(type);
				setDissociationType(std::make_optional(&d));
			}

			//Should this be an Else statement?  Not sure how to handle exception.  This is commented out because 
			//both of the toml config files I was testing with had many more key-value pairs than those accounted for above,
			//which led this to printing many Unrecognized parameter lines.
			// else {
			// 	std::cout << "Unrecognized parameter " << keyValuePair.first << " in file-specific parameters toml" << std::endl;
			// }
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

	toml::Table FileSpecificParameters::getFileSpecificParametersTomlTable()
	{
		toml::Table FileParameters;

		Tolerance *precursorMassTolerance = FileSpecificParameters::getPrecursorMassTolerance();
		Tolerance *productMassTolerance = FileSpecificParameters::getProductMassTolerance();
		Protease *protease = FileSpecificParameters::getProtease();
		std::optional<int> minPeptideLength = FileSpecificParameters::getMinPeptideLength();
		std::optional<int> maxPeptideLength = FileSpecificParameters::getMaxPeptideLength();
		std::optional<int> maxMissedCleavages = FileSpecificParameters::getMaxMissedCleavages();
		std::optional<int> maxModsforPeptide = FileSpecificParameters::getMaxModsForPeptide();
		std::optional<DissociationType*> dissociationType = FileSpecificParameters::getDissociationType();

		if (precursorMassTolerance != nullptr){
			FileParameters["PrecursorMassTolerance"] = std::to_string(precursorMassTolerance->getValue());
		}
		if (productMassTolerance != nullptr){
			FileParameters["ProductMassTolerance"] = std::to_string(productMassTolerance->getValue());
		}
		if (protease != nullptr){
			FileParameters["Protease"] = protease->getName();
			FileParameters["CleavageSpecificity"] = Proteomics::ProteolyticDigestion::CleavageSpecificityExtension::GetCleavageSpecificityAsString(protease->getCleavageSpecificity());
			FileParameters["PsiMsAccessionNumber"] = protease->getPsiMsAccessionNumber();
			FileParameters["PsiMsName"] = protease->getPsiMsName();
			
			//how to handle DigestionMotif?  Protease has a vector of DigestionMotifs.
			//Each DigestionMotif has 3 public string attributes:  InducingCleavage, PreventingCleavage, ExcludeFromWildcard
			//And one int attribute:  CutIndex
		}
		if (minPeptideLength.has_value()){
			FileParameters["MinPeptideLength"] = std::to_string(minPeptideLength.value());
		}
		if (maxPeptideLength.has_value()){
			FileParameters["MaxPeptideLength"] = std::to_string(maxPeptideLength.value());
		}
		if (maxMissedCleavages.has_value()){
			FileParameters["MaxMissedCleavages"] = std::to_string(maxMissedCleavages.value());
		}
		if (maxModsforPeptide.has_value()){
			FileParameters["MaxModsforPeptide"] = std::to_string(maxModsforPeptide.value());
		}
		if (dissociationType.has_value()){
			FileParameters["DissociationType"] = MassSpectrometry::GetDissocationType::GetDissocationTypeAsString(*dissociationType.value());
		}

		return FileParameters;
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
