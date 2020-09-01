#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../PeptideSpectralMatch.h"
#include "../CommonParameters.h"
#include "../MetaMorpheusEngineResults.h"

using namespace Chemistry;

namespace EngineLayer
{
    namespace ModificationAnalysis
    {

        typedef std::tuple<std::string, std::string, int> ModTuple;

        struct ModTuple_hash: public std::unary_function<ModTuple, std::size_t>{
            std::size_t operator() (const ModTuple& k ) const
            {
                size_t h1= std::hash<std::string>{}(std::get<0>(k));
                size_t h2= std::hash<std::string>{}(std::get<1>(k));
                size_t h3= std::hash<int>{}(std::get<2>(k));
                return h1 ^ (h2 << 1) ^ (h3 << 2);
            }
        };

        struct ModTuple_equal: public std::binary_function<ModTuple, ModTuple, bool>{
            bool operator() (const ModTuple& lhs, const ModTuple& rhs) const
            {
                return std::get<0>(lhs) == std::get<0>(rhs) &&
                    std::get<1>(lhs) == std::get<1>(rhs)    &&
                    std::get<2>(lhs) == std::get<2>(rhs);
            }
        };

        typedef std::unordered_set<ModTuple, ModTuple_hash, ModTuple_equal> ModTuple_set;

        class ModificationAnalysisEngine : public MetaMorpheusEngine
        {
        private:
            const std::vector<PeptideSpectralMatch*> NewPsms;
            
        public:
            ModificationAnalysisEngine(std::vector<PeptideSpectralMatch*> &newPsms, CommonParameters *commonParameters,
                                       std::vector<std::string> &nestedIds, int verbosityLevel=0);
            
        protected:
            MetaMorpheusEngineResults *RunSpecific() override;
        };
    }
}
