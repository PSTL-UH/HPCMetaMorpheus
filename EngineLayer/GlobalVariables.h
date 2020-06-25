#pragma once

#include "GlobalSettings.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include "stringhelper.h"
#include "tangible_filesystem.h"

#include "IGlobalSettings.h"

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"

//for toml files
#include "TomlReadWrite.h"

namespace EngineLayer
{
	class GlobalVariables final
	{
	private:
		static std::string privateDataDir;
		static bool privateStopLoops;
		static std::string privateElementsLocation;
		static std::string privateMetaMorpheusVersion;
		static GlobalSettings privateGlobalSettings;
		static std::vector<Modification*> privateUnimodDeserialized;
		static std::vector<Modification*> privateUniprotDeseralized;

		static std::unique_ptr<obo> privatePsiModDeserialized;
		static std::unordered_map<std::string, Modification*> privateAllModsKnownDictionary;
		static std::unordered_map<std::string, DissociationType> privateAllSupportedDissociationTypes;
		static std::string privateExperimentalDesignFileName;

		static std::vector<Modification*> _AllModsKnown;
		static std::unordered_set<std::string> _AllModTypesKnown;
                
                
	public:
                static void GlobalVariables_init();

		static std::vector<std::string> ErrorsReadingMods;

		// File locations
		static std::string getDataDir();
		static bool getStopLoops();
		static void setStopLoops(bool value);
		static std::string getElementsLocation();
		static std::string getMetaMorpheusVersion();
		static GlobalSettings getGlobalSettings();
		static std::vector<Modification*> getUnimodDeserialized();
		static std::vector<Modification*> getUniprotDeseralized();

		//static std::unique_ptr<obo> getPsiModDeserialized();
		static std::vector<Modification*> getAllModsKnown();

		static std::unordered_set<std::string> getAllModTypesKnown();
		static std::unordered_map<std::string, Modification*> getAllModsKnownDictionary();
		static void setAllModsKnownDictionary(const std::unordered_map<std::string, Modification*> &value);
		static std::unordered_map<std::string, DissociationType> getAllSupportedDissociationTypes();
		static void setAllSupportedDissociationTypes(const std::unordered_map<std::string, DissociationType> &value);

		static std::string getExperimentalDesignFileName();

		static void AddMods(std::vector<Modification*> &modifications, bool modsAreFromTheTopOfProteinXml);

		static std::string CheckLengthOfOutput(const std::string &psmString);
	};
}
