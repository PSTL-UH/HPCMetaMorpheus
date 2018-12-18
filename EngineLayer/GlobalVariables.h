#pragma once

#include "GlobalSettings.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "stringhelper.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class IGlobalSettings; }

using namespace MassSpectrometry;
using namespace Nett;
using namespace Proteomics;

namespace EngineLayer
{
	class GlobalVariables final
	{
	private:
		static std::wstring privateDataDir;
		static bool privateStopLoops;
		static std::wstring privateElementsLocation;
		static std::wstring privateMetaMorpheusVersion;
		static IGlobalSettings *privateGlobalSettings;
		static std::vector<Modification*> privateUnimodDeserialized;
		static std::vector<Modification*> privateUniprotDeseralized;
		static UsefulProteomicsDatabases::Generated::obo *privatePsiModDeserialized;
		static std::unordered_map<std::wstring, Modification*> privateAllModsKnownDictionary;
		static std::unordered_map<std::wstring, DissociationType*> privateAllSupportedDissociationTypes;
		static std::wstring privateExperimentalDesignFileName;

		static std::vector<Modification*> _AllModsKnown;
		static std::unordered_set<std::wstring> _AllModTypesKnown;

	private:
		class StaticConstructor
		{
		public:
			StaticConstructor();
		};

	private:
		static GlobalVariables::StaticConstructor staticConstructor;


	public:
		static std::vector<std::wstring> ErrorsReadingMods;
		// File locations
		static std::wstring getDataDir();
		static bool getStopLoops();
		static void setStopLoops(bool value);
		static std::wstring getElementsLocation();
		static std::wstring getMetaMorpheusVersion();
		static IGlobalSettings *getGlobalSettings();
		static std::vector<Modification*> getUnimodDeserialized();
		static std::vector<Modification*> getUniprotDeseralized();
		static UsefulProteomicsDatabases::Generated::obo *getPsiModDeserialized();
		static std::vector<Modification*> getAllModsKnown();
		static std::vector<std::wstring> getAllModTypesKnown();
		static std::unordered_map<std::wstring, Modification*> getAllModsKnownDictionary();
		static void setAllModsKnownDictionary(const std::unordered_map<std::wstring, Modification*> &value);
		static std::unordered_map<std::wstring, DissociationType*> getAllSupportedDissociationTypes();
		static void setAllSupportedDissociationTypes(const std::unordered_map<std::wstring, DissociationType*> &value);

		static std::wstring getExperimentalDesignFileName();

		static void AddMods(std::vector<Modification*> &modifications, bool modsAreFromTheTopOfProteinXml);

		static std::wstring CheckLengthOfOutput(const std::wstring &psmString);
	};
}
