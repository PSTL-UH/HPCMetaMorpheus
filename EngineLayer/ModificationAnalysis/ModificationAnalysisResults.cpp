#include "ModificationAnalysisResults.h"
#include "ModificationAnalysisEngine.h"

using namespace Chemistry;
namespace EngineLayer
{
    namespace ModificationAnalysis
    {
        
        ModificationAnalysisResults::ModificationAnalysisResults(ModificationAnalysisEngine *modificationAnalysisEngine) : MetaMorpheusEngineResults(modificationAnalysisEngine)
        {
        }
        
        std::unordered_map<std::string, int> ModificationAnalysisResults::getCountOfAmbiguousButLocalizedModsSeen() const
        {
            return privateCountOfAmbiguousButLocalizedModsSeen;
        }
        
        void ModificationAnalysisResults::setCountOfAmbiguousButLocalizedModsSeen(const std::unordered_map<std::string, int> &value)
        {
            privateCountOfAmbiguousButLocalizedModsSeen = value;
        }
        
        std::unordered_map<std::string, int> ModificationAnalysisResults::getCountOfModsSeenAndLocalized() const
        {
            return privateCountOfModsSeenAndLocalized;
        }
        
        void ModificationAnalysisResults::setCountOfModsSeenAndLocalized(const std::unordered_map<std::string, int> &value)
        {
            privateCountOfModsSeenAndLocalized = value;
        }
        
        std::unordered_map<std::string, int> ModificationAnalysisResults::getCountOfEachModSeenOnProteins() const
        {
            return privateCountOfEachModSeenOnProteins;
        }
        
        void ModificationAnalysisResults::setCountOfEachModSeenOnProteins(const std::unordered_map<std::string, int> &value)
        {
            privateCountOfEachModSeenOnProteins = value;
        }
        
        std::unordered_map<std::string, int> ModificationAnalysisResults::getCountOfUnlocalizedMods() const
        {
            return privateCountOfUnlocalizedMods;
        }
        
        void ModificationAnalysisResults::setCountOfUnlocalizedMods(const std::unordered_map<std::string, int> &value)
        {
            privateCountOfUnlocalizedMods = value;
        }
        
        std::unordered_map<ChemicalFormula*, int> ModificationAnalysisResults::getCountOfUnlocalizedFormulas() const
        {
            return privateCountOfUnlocalizedFormulas;
        }
        
        void ModificationAnalysisResults::setCountOfUnlocalizedFormulas(const std::unordered_map<ChemicalFormula*, int> &value)
        {
            privateCountOfUnlocalizedFormulas = value;
        }
        
        std::string ModificationAnalysisResults::ToString()
        {
            auto sb = new StringBuilder();

            sb->appendLine(MetaMorpheusEngineResults::ToString());
            sb->appendLine("Localized mods seen below q-value 0.01:");
#ifdef ORIG
            sb->appendLine(std::string::Join("\r\n", getCountOfModsSeenAndLocalized().OrderBy([&] (std::any b) {
                            -b->Value;
                        })->Select([&] (std::any b)   {
				return "\t" + b::Key + "\t" + b->Value;
                            })));
#endif
            std::vector<std::tuple<std::string, int>> tmp;
            for ( auto p = getCountOfModsSeenAndLocalized().begin(); p != getCountOfModsSeenAndLocalized().end(); p++ ) {
                tmp.push_back(std::make_tuple(std::get<0>(*p), std::get<1>(*p)));
            }
            std::sort(tmp.begin(), tmp.end(), [&] (std::tuple<std::string, int> l, std::tuple<std::string, int> r) {
                    return -std::get<1>(l) < -std::get<1>(r);
                });
            std::string stmp = "\t";
            std::string stmp2 = "\r\n";
            std::string sp;
            for ( int p = 0; p < (int)tmp.size(); p++ ) {               
                sp += stmp + std::get<0>(tmp[p]) + stmp + std::to_string(std::get<1>(tmp[p]));
                if ( p != ((int)tmp.size() - 1) ) 
                    sp += stmp2;
            }
            sb->appendLine (sp);
            
            sb->appendLine("(Approx) Additional localized but protein ambiguous mods seen below q-value 0.01:");
#ifdef ORIG
            sb->appendLine(std::string::Join("\r\n", getCountOfAmbiguousButLocalizedModsSeen().OrderBy([&] (std::any b){
                            -b->Value;                            
			})->Select([&] (std::any b)  {
				return "\t" + b::Key + "\t" + b->Value;
                            })));
#endif
            tmp.clear();
            for ( auto p = getCountOfAmbiguousButLocalizedModsSeen().begin(); p != getCountOfAmbiguousButLocalizedModsSeen().end(); p++ ) {
                tmp.push_back(std::make_tuple(std::get<0>(*p), std::get<1>(*p)));
            }
            std::sort(tmp.begin(), tmp.end(), [&] (std::tuple<std::string, int> l, std::tuple<std::string, int> r) {
                    return -std::get<1>(l) < -std::get<1>(r);
                });
            sp="";
            for ( int p = 0; p < (int)tmp.size(); p++ ) {               
                sp += stmp + std::get<0>(tmp[p]) + stmp + std::to_string(std::get<1>(tmp[p]));
                if ( p != ((int)tmp.size() - 1) ) 
                    sp += stmp2;
            }
            sb->appendLine (sp);
            
            sb->appendLine("(Approx) Additional unlocalized mods seen below q-value 0.01:");
#ifdef ORIG
            sb->appendLine(std::string::Join("\r\n", getCountOfUnlocalizedMods().OrderBy([&] (std::any b) {
                            -b->Value;
			})->Select([&] (std::any b){
				return "\t" + b::Key + "\t" + b->Value;
                            })));
#endif
            tmp.clear();
            for ( auto p = getCountOfUnlocalizedMods().begin(); p != getCountOfUnlocalizedMods().end(); p++ ) {
                tmp.push_back(std::make_tuple(std::get<0>(*p), std::get<1>(*p)));
            }
            std::sort(tmp.begin(), tmp.end(), [&] (std::tuple<std::string, int> l, std::tuple<std::string, int> r) {
                    return -std::get<1>(l) < -std::get<1>(r);
                });
            sp="";
            for ( int p = 0; p < (int)tmp.size(); p++ ) {               
                sp += stmp + std::get<0>(tmp[p]) + stmp + std::to_string(std::get<1>(tmp[p]));
                if ( p != ((int)tmp.size() - 1) ) 
                    sp += stmp2;
            }
            sb->appendLine (sp);

            sb->appendLine("(Approx) Additional unlocalized modification formulas seen below q-value 0.01:");
#ifdef ORIG
            sb->appendLine(std::string::Join("\r\n", getCountOfUnlocalizedFormulas().OrderBy([&] (std::any b){
                            -b->Value;
			})->Select([&] (std::any b){
				return "\t" + b::Key->Formula + "\t" + b->Value;
                            })));
#endif
            std::vector<std::tuple<ChemicalFormula*, int>> tmp2;            
            for ( auto p = getCountOfUnlocalizedFormulas().begin(); p != getCountOfUnlocalizedFormulas().end(); p++ ) {
                tmp2.push_back(std::make_tuple(std::get<0>(*p), std::get<1>(*p)));
            }
            std::sort(tmp2.begin(), tmp2.end(), [&] (std::tuple<ChemicalFormula*, int> l, std::tuple<ChemicalFormula*, int> r) {
                    return -std::get<1>(l) < -std::get<1>(r);
                });
            sp="";
            for ( int p = 0; p < (int)tmp2.size(); p++ ) {               
                sp += stmp + std::get<0>(tmp2[p])->getFormula() + stmp + std::to_string(std::get<1>(tmp2[p]));
                if ( p != ((int)tmp2.size() - 1) ) 
                    sp += stmp2;
            }
            sb->appendLine (sp);


            sb->appendLine();
            sb->appendLine("All mods in database limited to peptides observed in the results:");
#ifdef ORIG
            sb->appendLine(std::string::Join("\r\n", getCountOfEachModSeenOnProteins().OrderBy([&] (std::any b)	{
                            -b->Value;
			})->Select([&] (std::any b){
				return "\t" + b::Key + "\t" + b->Value;
                            })));
#endif
            tmp.clear();
            for ( auto p = getCountOfEachModSeenOnProteins().begin(); p != getCountOfEachModSeenOnProteins().end(); p++ ) {
                tmp.push_back(std::make_tuple(std::get<0>(*p), std::get<1>(*p)));
            }
            std::sort(tmp.begin(), tmp.end(), [&] (std::tuple<std::string, int> l, std::tuple<std::string, int> r) {
                    return -std::get<1>(l) < -std::get<1>(r);
                });
            sp="";
            for ( int p = 0; p < (int)tmp.size(); p++ ) {               
                sp += stmp + std::get<0>(tmp[p]) + stmp + std::to_string(std::get<1>(tmp[p]));
                if ( p != ((int)tmp.size() - 1) ) 
                    sp += stmp2;
            }
            sb->appendLine (sp);
            
            std::string s = sb->toString();
            delete sb;
            return s;
        }
    }
}
