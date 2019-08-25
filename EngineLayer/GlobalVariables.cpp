#include "GlobalVariables.h"
#include "IGlobalSettings.h"

#include <experimental/filesystem>

using namespace MassSpectrometry;
//using namespace Nett;
using namespace Proteomics;

namespace EngineLayer
{

    std::string GlobalVariables::privateDataDir;
    bool GlobalVariables::privateStopLoops = false;
    std::string GlobalVariables::privateElementsLocation;
    std::string GlobalVariables::privateMetaMorpheusVersion;
    IGlobalSettings *GlobalVariables::privateGlobalSettings;
    std::vector<Modification*> GlobalVariables::privateUnimodDeserialized;
    std::vector<Modification*> GlobalVariables::privateUniprotDeseralized;
    UsefulProteomicsDatabases::Generated::obo *GlobalVariables::privatePsiModDeserialized;
    std::unordered_map<std::string, Modification*> GlobalVariables::privateAllModsKnownDictionary;
    std::unordered_map<std::string, DissociationType> GlobalVariables::privateAllSupportedDissociationTypes;
    std::string GlobalVariables::privateExperimentalDesignFileName;
    std::vector<Modification*> GlobalVariables::_AllModsKnown;
    std::unordered_set<std::string> GlobalVariables::_AllModTypesKnown;
    
    GlobalVariables::StaticConstructor::StaticConstructor()
    {
        privateMetaMorpheusVersion = "0.0.295";
        
        if (getMetaMorpheusVersion() == "1.0.0.0")
        {
#if defined(DEBUG)
            privateMetaMorpheusVersion += "Not a release version. DEBUG.";
#else
            privateMetaMorpheusVersion += "Not a release version.";
#endif
        }
        else
        {
            // as of 0.0.277, AppVeyor appends the build number
            // this is intentional; it's to avoid conflicting AppVeyor build numbers
            // trim the build number off the version number for displaying/checking versions, etc
            auto foundIndexes = std::vector<int>();
            for (int i = 0; i < (int)getMetaMorpheusVersion().length(); i++)
            {
                if (getMetaMorpheusVersion()[i] == L'.')
                {
                    foundIndexes.push_back(i);
                }
            }
            privateMetaMorpheusVersion = getMetaMorpheusVersion().substr(0, foundIndexes.back());
        }
        
#ifdef ORIG
        // Edgar: directory structure for C++/Linux not as sophisticated for now. Can revisit later if necessary.
        auto pathToProgramFiles = Environment::GetFolderPath(Environment::SpecialFolder::ProgramFiles);
        if (!StringHelper::isEmptyOrWhiteSpace(pathToProgramFiles)                                &&
            AppDomain::CurrentDomain->BaseDirectory.find(pathToProgramFiles) != std::string::npos &&
            !AppDomain::CurrentDomain->BaseDirectory.find("Jenkins") != std::string::npos)
        {
            privateDataDir = FileSystem::combine(Environment::GetFolderPath(Environment::SpecialFolder::LocalApplicationData), "MetaMorpheus");
        }
        else
        {
            privateDataDir = AppDomain::CurrentDomain->BaseDirectory;
        }
#endif
        privateDataDir = std::experimental::filesystem::current_path();
        
        std::string datadir = FileSystem::combine(getDataDir(), "Data" );
        privateElementsLocation = FileSystem::combine(datadir, "elements.dat");
        UsefulProteomicsDatabases::Loaders::LoadElements(getElementsLocation());
        
        privateExperimentalDesignFileName = "ExperimentalDesign.tsv";

        privateUnimodDeserialized = UsefulProteomicsDatabases::Loaders::LoadUnimod(FileSystem::combine(datadir, "unimod.xml")); //.ToList();
        privatePsiModDeserialized = UsefulProteomicsDatabases::Loaders::LoadPsiMod(FileSystem::combine(datadir, "PSI-MOD.obo.xml"));
        auto formalChargesDictionary = UsefulProteomicsDatabases::Loaders::GetFormalChargesDictionary(getPsiModDeserialized());
        privateUniprotDeseralized = UsefulProteomicsDatabases::Loaders::LoadUniprot(FileSystem::combine(datadir, "ptmlist.txt"), formalChargesDictionary); //.ToList();
        
        //for (auto modFile : Directory::GetFiles(FileSystem::combine(getDataDir(), "Mods")))
        for (auto  modFile : std::experimental::filesystem::directory_iterator(FileSystem::combine(getDataDir(), "Mods")))
        {
            std::vector<std::tuple<Modification *, std::string>> errorMods;
            std::vector<Modification *> mods;
            mods = UsefulProteomicsDatabases::PtmListLoader::ReadModsFromFile(modFile.path().string(), errorMods);
            AddMods( mods, false);
        }

#ifdef ORIG
        AddMods(getUniprotDeserialized().OfType<Modification*>(), false);
        AddMods(getUnimodDeserialized().OfType<Modification*>(), false);
#endif
        AddMods(privateUniprotDeseralized, false);
        AddMods(privateUnimodDeserialized, false);
        
        // populate dictionaries of known mods/proteins for deserialization
        setAllModsKnownDictionary(std::unordered_map<std::string, Modification*>());
        for (auto mod : getAllModsKnown())
        {
            if (getAllModsKnownDictionary().find(mod->getIdWithMotif()) == getAllModsKnownDictionary().end())
            {
                getAllModsKnownDictionary().emplace(mod->getIdWithMotif(), mod);
            }
            // no error thrown if multiple mods with this ID are present - just pick one
        }
        
#ifdef ORIG
        // Todo for Nick: Replace the Nett functionality here.
        GlobalSettings = Toml::ReadFile<getGlobalSettings()*>(FileSystem::combine(getDataDir(), "settings.toml"));
#endif

#ifdef ORIG
        setAllSupportedDissociationTypes(std::unordered_map<std::string, DissociationType> {
                {DissociationType::CID.ToString(), DissociationType::CID},
                {DissociationType::ECD.ToString(), DissociationType::ECD},
                {DissociationType::ETD.ToString(), DissociationType::ETD},
                {DissociationType::HCD.ToString(), DissociationType::HCD},
                {DissociationType::EThcD.ToString(), DissociationType::EThcD}
            });
#endif
        setAllSupportedDissociationTypes(std::unordered_map<std::string, DissociationType> {
                {"CID", DissociationType::CID},
                {"ECD", DissociationType::ECD},
                {"ETD", DissociationType::ETD},
                {"HCD", DissociationType::HCD},
                {"EThcD", DissociationType::EThcD}
            });
    }
    
    //GlobalVariables::StaticConstructor GlobalVariables::staticConstructor;
    //std::vector<std::string> GlobalVariables::ErrorsReadingMods;
    
    std::string GlobalVariables::getDataDir()
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
    
    std::string GlobalVariables::getElementsLocation()
    {
        return privateElementsLocation;
    }
    
    std::string GlobalVariables::getMetaMorpheusVersion()
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
        //return _AllModsKnown.AsEnumerable();
        return _AllModsKnown;
    }
    
    //std::vector<std::string> GlobalVariables::getAllModTypesKnown()
    std::unordered_set<std::string> GlobalVariables::getAllModTypesKnown()
    {
        //return _AllModTypesKnown.AsEnumerable();
        return _AllModTypesKnown;
    }
    
    std::unordered_map<std::string, Modification*> GlobalVariables::getAllModsKnownDictionary()
    {
        return privateAllModsKnownDictionary;
    }
    
    void GlobalVariables::setAllModsKnownDictionary(const std::unordered_map<std::string, Modification*> &value)
    {
        privateAllModsKnownDictionary = value;
    }
    
    std::unordered_map<std::string, DissociationType> GlobalVariables::getAllSupportedDissociationTypes()
    {
        return privateAllSupportedDissociationTypes;
    }
    
    void GlobalVariables::setAllSupportedDissociationTypes(const std::unordered_map<std::string, DissociationType> &value)
    {
        privateAllSupportedDissociationTypes = value;
    }
    
    std::string GlobalVariables::getExperimentalDesignFileName()
    {
        return privateExperimentalDesignFileName;
    }
    
    void GlobalVariables::AddMods(std::vector<Modification*> &modifications, bool modsAreFromTheTopOfProteinXml)
    {
        for (auto mod : modifications)
        {
            if (mod->getModificationType().empty() || mod->getIdWithMotif().empty())
            {
                ErrorsReadingMods.push_back(mod->ToString() + "\r\n" + " has null or empty modification type");
                continue;
            }
#ifdef ORIG
            //if (getAllModsKnown().Any([&] (std::any b) {
            //            return b::IdWithMotif->Equals(mod->IdWithMotif) &&
            //                b::ModificationType->Equals(mod->ModificationType) && !b->Equals(mod);
            //}))
#endif
            bool found = false;
            bool found2 = false;
            bool found3 = false;
            for  ( auto b: getAllModsKnown() ) {
                if ( b->getIdWithMotif() == mod->getIdWithMotif() ) {
                    found3 = true;
                    if ( b->getModificationType() == mod->getModificationType() ) {
                        found2 = true;
                        if ( !b->Equals(mod) ) {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if ( found ) {
                if (modsAreFromTheTopOfProteinXml)
                {
#ifdef ORIG
                    _AllModsKnown.RemoveAll([&] (std::any p) {
                            return p::IdWithMotif->Equals(mod->IdWithMotif) &&
                                p::ModificationType->Equals(mod->ModificationType) && !p->Equals(mod);
                        });
#endif
                    for ( auto p = _AllModsKnown.begin(); p != _AllModsKnown.end(); p++ ) {
                        if ( (*p)->getIdWithMotif() == mod->getIdWithMotif()              &&
                             (*p)->getModificationType() == mod->getModificationType()  &&
                             !(*p)->Equals(mod) ) {
                            _AllModsKnown.erase(p);
                        }
                    }
                    _AllModsKnown.push_back(mod);
                    _AllModTypesKnown.insert(mod->getModificationType());
                }
                else
                {
                    ErrorsReadingMods.push_back(std::string("Modification id and type are equal, but some fields are not! ") +
                                                "The following mod was not read in: " + "\r\n" + mod->ToString());
                }
                continue;
            }
#ifdef ORIG
            //else if (getAllModsKnown().Any([&] (std::any b) {
            //            return b::IdWithMotif->Equals(mod->IdWithMotif) &&
            //                b::ModificationType->Equals(mod->ModificationType);
            //        }))
#endif
            else if ( found2 ) {
                // same ID, same mod type, and same mod properties; continue and don't output an error message
                // this could result from reading in an XML database with mods annotated at the top
                // that are already loaded in MetaMorpheus
                continue;
            }
#ifdef ORIG
            //else if (getAllModsKnown().Any([&] (std::any m) {
            //            return m->IdWithMotif == mod->getIdWithMotif();
            //        }))
#endif            
            else if ( found3 ) {
                // same ID but different mod types. This can happen if the user names a mod the same as a UniProt mod
                // this is problematic because if a mod is annotated in the database, all we have to go on is an ID ("description" tag).
                // so we don't know which mod to use, causing unnecessary ambiguity
                if (modsAreFromTheTopOfProteinXml)
                {
#ifdef ORIG
                    _AllModsKnown.RemoveAll([&] (std::any p) {
                            return p::IdWithMotif->Equals(mod->getIdWithMotif()) && !p->Equals(mod);
                        });
#endif
                    for ( auto p = _AllModsKnown.begin(); p != _AllModsKnown.end(); p++ ) {
                        if ( (*p)->getIdWithMotif() == mod->getIdWithMotif()   &&
                             !(*p)->Equals(mod) ) {
                            _AllModsKnown.erase(p);
                        }
                    }
                    _AllModsKnown.push_back(mod);
                    _AllModTypesKnown.insert(mod->getModificationType());
                }
                else if (mod->getModificationType() != ("Unimod"))
                {
                    ErrorsReadingMods.push_back("Duplicate mod IDs! Skipping " + mod->getModificationType() + ":" + mod->getIdWithMotif());
                }
                continue;
            }
            else
            {
                // no errors! add the mod
                _AllModsKnown.push_back(mod);
                _AllModTypesKnown.insert(mod->getModificationType());
            }
        }
    }
    
    std::string GlobalVariables::CheckLengthOfOutput(const std::string &psmString)
    {
        if (psmString.length() > 32000 && getGlobalSettings()->getWriteExcelCompatibleTSVs())
        {
            return "Output too long for Excel";
        }
        else
        {
            return psmString;
        }
    }
}
