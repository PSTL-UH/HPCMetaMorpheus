#include "GlobalVariables.h"
#include "IGlobalSettings.h"

using namespace MassSpectrometry;
using namespace Nett;
using namespace Proteomics;

namespace EngineLayer
{

std::wstring GlobalVariables::privateDataDir;
bool GlobalVariables::privateStopLoops = false;
std::wstring GlobalVariables::privateElementsLocation;
std::wstring GlobalVariables::privateMetaMorpheusVersion;
IGlobalSettings *GlobalVariables::privateGlobalSettings;
std::vector<Modification*> GlobalVariables::privateUnimodDeserialized;
std::vector<Modification*> GlobalVariables::privateUniprotDeseralized;
UsefulProteomicsDatabases::Generated::obo *GlobalVariables::privatePsiModDeserialized;
std::unordered_map<std::wstring, Modification*> GlobalVariables::privateAllModsKnownDictionary;
std::unordered_map<std::wstring, DissociationType*> GlobalVariables::privateAllSupportedDissociationTypes;
std::wstring GlobalVariables::privateExperimentalDesignFileName;
std::vector<Modification*> GlobalVariables::_AllModsKnown;
std::unordered_set<std::wstring> GlobalVariables::_AllModTypesKnown;

	GlobalVariables::StaticConstructor::StaticConstructor()
	{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		MetaMorpheusVersion = GlobalVariables::typeid->Assembly->GetName()->Version->ToString();
        
		if (getMetaMorpheusVersion() == L"1.0.0.0")
		{
		#if defined(DEBUG)
			MetaMorpheusVersion = L"Not a release version. DEBUG.";
		#else
			MetaMorpheusVersion = L"Not a release version.";
		#endif
		}
		else
		{
			// as of 0.0.277, AppVeyor appends the build number
			// this is intentional; it's to avoid conflicting AppVeyor build numbers
			// trim the build number off the version number for displaying/checking versions, etc
			auto foundIndexes = std::vector<int>();
			for (int i = 0; i < getMetaMorpheusVersion().length(); i++)
			{
				if (getMetaMorpheusVersion()[i] == L'.')
				{
					foundIndexes.push_back(i);
				}
			}
			MetaMorpheusVersion = getMetaMorpheusVersion().substr(0, foundIndexes.back());
		}
        
		{
			auto pathToProgramFiles = Environment::GetFolderPath(Environment::SpecialFolder::ProgramFiles);
			if (!StringHelper::isEmptyOrWhiteSpace(pathToProgramFiles) && AppDomain::CurrentDomain->BaseDirectory.find(pathToProgramFiles) != std::wstring::npos && !AppDomain::CurrentDomain->BaseDirectory.find(L"Jenkins") != std::wstring::npos)
			{
				DataDir = FileSystem::combine(Environment::GetFolderPath(Environment::SpecialFolder::LocalApplicationData), L"MetaMorpheus");
			}
			else
			{
				DataDir = AppDomain::CurrentDomain->BaseDirectory;
			}
		}
        
		ElementsLocation = FileSystem::combine(getDataDir(), LR"(Data)", LR"(elements.dat)");
		UsefulProteomicsDatabases::Loaders::LoadElements(getElementsLocation());
        
		ExperimentalDesignFileName = L"ExperimentalDesign.tsv";
        
		UnimodDeserialized = UsefulProteomicsDatabases::Loaders::LoadUnimod(FileSystem::combine(getDataDir(), LR"(Data)", LR"(unimod.xml)")).ToList();
		PsiModDeserialized = UsefulProteomicsDatabases::Loaders::LoadPsiMod(FileSystem::combine(getDataDir(), LR"(Data)", LR"(PSI-MOD.obo.xml)"));
		auto formalChargesDictionary = UsefulProteomicsDatabases::Loaders::GetFormalChargesDictionary(getPsiModDeserialized());
		UniprotDeseralized = UsefulProteomicsDatabases::Loaders::LoadUniprot(FileSystem::combine(getDataDir(), LR"(Data)", LR"(ptmlist.txt)"), formalChargesDictionary).ToList();
        
		for (auto modFile : Directory::GetFiles(FileSystem::combine(getDataDir(), LR"(Mods)")))
		{
			std::any errorMods;
			AddMods(UsefulProteomicsDatabases::PtmListLoader::ReadModsFromFile(modFile, errorMods), false);
		}
        
		AddMods(getUniprotDeseralized().OfType<Modification*>(), false);
		AddMods(getUnimodDeserialized().OfType<Modification*>(), false);
        
		// populate dictionaries of known mods/proteins for deserialization
		setAllModsKnownDictionary(std::unordered_map<std::wstring, Modification*>());
		for (auto mod : getAllModsKnown())
		{
			if (getAllModsKnownDictionary().find(mod->IdWithMotif) == getAllModsKnownDictionary().end())
			{
				getAllModsKnownDictionary().emplace(mod->IdWithMotif, mod);
			}
			// no error thrown if multiple mods with this ID are present - just pick one
		}
        
		GlobalSettings = Toml::ReadFile<getGlobalSettings()*>(FileSystem::combine(getDataDir(), LR"(settings.toml)"));
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		setAllSupportedDissociationTypes(std::unordered_map<std::wstring, DissociationType*>
		{
			{DissociationType::CID.ToString(), DissociationType::CID},
			{DissociationType::ECD.ToString(), DissociationType::ECD},
			{DissociationType::ETD.ToString(), DissociationType::ETD},
			{DissociationType::HCD.ToString(), DissociationType::HCD},
			{DissociationType::EThcD.ToString(), DissociationType::EThcD}
		});
	}

GlobalVariables::StaticConstructor GlobalVariables::staticConstructor;
std::vector<std::wstring> GlobalVariables::ErrorsReadingMods;

	std::wstring GlobalVariables::getDataDir()
	{
		return privateDataDir;
	}

	bool GlobalVariables::getStopLoops()
	{
		return privateStopLoops;
	}

	void GlobalVariables::setStopLoops(bool value)
	{
		privateStopLoops = value;
	}

	std::wstring GlobalVariables::getElementsLocation()
	{
		return privateElementsLocation;
	}

	std::wstring GlobalVariables::getMetaMorpheusVersion()
	{
		return privateMetaMorpheusVersion;
	}

	IGlobalSettings *GlobalVariables::getGlobalSettings()
	{
		return privateGlobalSettings;
	}

	std::vector<Modification*> GlobalVariables::getUnimodDeserialized()
	{
		return privateUnimodDeserialized;
	}

	std::vector<Modification*> GlobalVariables::getUniprotDeseralized()
	{
		return privateUniprotDeseralized;
	}

	UsefulProteomicsDatabases::Generated::obo *GlobalVariables::getPsiModDeserialized()
	{
		return privatePsiModDeserialized;
	}

	std::vector<Modification*> GlobalVariables::getAllModsKnown()
	{
		return _AllModsKnown.AsEnumerable();
	}

	std::vector<std::wstring> GlobalVariables::getAllModTypesKnown()
	{
		return _AllModTypesKnown.AsEnumerable();
	}

	std::unordered_map<std::wstring, Modification*> GlobalVariables::getAllModsKnownDictionary()
	{
		return privateAllModsKnownDictionary;
	}

	void GlobalVariables::setAllModsKnownDictionary(const std::unordered_map<std::wstring, Modification*> &value)
	{
		privateAllModsKnownDictionary = value;
	}

	std::unordered_map<std::wstring, DissociationType*> GlobalVariables::getAllSupportedDissociationTypes()
	{
		return privateAllSupportedDissociationTypes;
	}

	void GlobalVariables::setAllSupportedDissociationTypes(const std::unordered_map<std::wstring, DissociationType*> &value)
	{
		privateAllSupportedDissociationTypes = value;
	}

	std::wstring GlobalVariables::getExperimentalDesignFileName()
	{
		return privateExperimentalDesignFileName;
	}

	void GlobalVariables::AddMods(std::vector<Modification*> &modifications, bool modsAreFromTheTopOfProteinXml)
	{
		for (auto mod : modifications)
		{
			if (mod->ModificationType.empty() || mod->IdWithMotif.empty())
			{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				ErrorsReadingMods.push_back(mod->ToString() + L"\r\n" + L" has null or empty modification type");
				continue;
			}
			if (getAllModsKnown().Any([&] (std::any b)
			{
				return b::IdWithMotif->Equals(mod->IdWithMotif) && b::ModificationType->Equals(mod->ModificationType) && !b->Equals(mod);
			}))
			{
				if (modsAreFromTheTopOfProteinXml)
				{
					_AllModsKnown.RemoveAll([&] (std::any p)
					{
						return p::IdWithMotif->Equals(mod->IdWithMotif) && p::ModificationType->Equals(mod->ModificationType) && !p->Equals(mod);
					});
					_AllModsKnown.push_back(mod);
					_AllModTypesKnown.insert(mod->ModificationType);
				}
				else
				{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					ErrorsReadingMods.push_back(std::wstring(L"Modification id and type are equal, but some fields are not! ") + L"The following mod was not read in: " + L"\r\n" + mod->ToString());
				}
				continue;
			}
			else if (getAllModsKnown().Any([&] (std::any b)
			{
				return b::IdWithMotif->Equals(mod->IdWithMotif) && b::ModificationType->Equals(mod->ModificationType);
			}))
			{
				// same ID, same mod type, and same mod properties; continue and don't output an error message
				// this could result from reading in an XML database with mods annotated at the top
				// that are already loaded in MetaMorpheus
				continue;
			}
			else if (getAllModsKnown().Any([&] (std::any m)
			{
				return m->IdWithMotif == mod->IdWithMotif;
			}))
			{
				// same ID but different mod types. This can happen if the user names a mod the same as a UniProt mod
				// this is problematic because if a mod is annotated in the database, all we have to go on is an ID ("description" tag).
				// so we don't know which mod to use, causing unnecessary ambiguity
				if (modsAreFromTheTopOfProteinXml)
				{
					_AllModsKnown.RemoveAll([&] (std::any p)
					{
						return p::IdWithMotif->Equals(mod->IdWithMotif) && !p->Equals(mod);
					});
					_AllModsKnown.push_back(mod);
					_AllModTypesKnown.insert(mod->ModificationType);
				}
				else if (!mod->ModificationType->Equals(L"Unimod"))
				{
					ErrorsReadingMods.push_back(L"Duplicate mod IDs! Skipping " + mod->ModificationType + L":" + mod->IdWithMotif);
				}
				continue;
			}
			else
			{
				// no errors! add the mod
				_AllModsKnown.push_back(mod);
				_AllModTypesKnown.insert(mod->ModificationType);
			}
		}
	}

	std::wstring GlobalVariables::CheckLengthOfOutput(const std::wstring &psmString)
	{
		if (psmString.length() > 32000 && getGlobalSettings()->getWriteExcelCompatibleTSVs())
		{
			return L"Output too long for Excel";
		}
		else
		{
			return psmString;
		}
	}
}
