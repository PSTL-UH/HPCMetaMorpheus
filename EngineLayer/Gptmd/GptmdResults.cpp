#include "GptmdResults.h"
#include "../MetaMorpheusEngine.h"

using namespace Proteomics;
namespace EngineLayer
{
	namespace Gptmd
	{

		GptmdResults::GptmdResults(MetaMorpheusEngine *s, std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> &mods, int modsAdded) : MetaMorpheusEngineResults(s), ModsAdded(modsAdded)
		{
			setMods(mods);
		}

		std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> GptmdResults::getMods() const
		{
			return privateMods;
		}

		void GptmdResults::setMods(const std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> &value)
		{
			privateMods = value;
		}

		std::wstring GptmdResults::ToString()
		{
			auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->appendLine(MetaMorpheusEngineResults::ToString());
			sb->appendLine(L"Modifications trying to add: " + std::to_wstring(ModsAdded));
			sb->appendLine(L"Proteins trying to expand: " + std::to_wstring(getMods().size()));
			sb->appendLine(L"Mods types and counts:");
			sb->append(std::wstring::Join(L"\r\n", getMods().SelectMany([&] (std::any b)
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
				return L"\t" + b::Key->IdWithMotif + L"\t" + std::to_wstring(b->Count());
			})));

			delete sb;
			return sb->toString();
		}
	}
}
