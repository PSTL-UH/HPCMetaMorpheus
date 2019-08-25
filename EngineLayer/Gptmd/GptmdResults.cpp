#include "GptmdResults.h"
#include "../MetaMorpheusEngine.h"

using namespace Proteomics;
namespace EngineLayer
{
	namespace Gptmd
	{

		GptmdResults::GptmdResults(MetaMorpheusEngine *s, std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>> &mods, int modsAdded) : MetaMorpheusEngineResults(s), ModsAdded(modsAdded)
		{
			setMods(mods);
		}

		std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>> GptmdResults::getMods() const
		{
			return privateMods;
		}

		void GptmdResults::setMods(const std::unordered_map<std::string, std::unordered_set<std::tuple<int, Modification*>>> &value)
		{
			privateMods = value;
		}

		std::string GptmdResults::ToString()
		{
			auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->appendLine(MetaMorpheusEngineResults::ToString());
			sb->appendLine("Modifications trying to add: " + std::to_string(ModsAdded));
			sb->appendLine("Proteins trying to expand: " + std::to_string(getMods().size()));
			sb->appendLine("Mods types and counts:");
			sb->append(std::string::Join("\r\n", getMods().SelectMany([&] (std::any b)
			{
				b->Value;
			}).GroupBy([&] (std::any b)
			{
				b::Item2;
			}).OrderBy([&] (std::any b)
			{
				-b->Count();
			})->Select([&] (std::any b)
			{
			delete sb;
				return "\t" + b::Key->IdWithMotif + "\t" + std::to_string(b->Count());
			})));

			delete sb;
			return sb->toString();
		}
	}
}
