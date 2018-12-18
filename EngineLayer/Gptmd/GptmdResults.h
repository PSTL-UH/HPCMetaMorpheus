#pragma once

#include "../MetaMorpheusEngineResults.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include "stringbuilder.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MetaMorpheusEngine; }

using namespace Proteomics;

namespace EngineLayer
{
	namespace Gptmd
	{
		class GptmdResults : public MetaMorpheusEngineResults
		{
		private:
			std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> privateMods;

			const int ModsAdded;

		public:
			GptmdResults(MetaMorpheusEngine *s, std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> &mods, int modsAdded);

				std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> getMods() const;
				void setMods(const std::unordered_map<std::wstring, std::unordered_set<std::tuple<int, Modification*>>> &value);

			std::wstring ToString() override;
		};
	}
}
