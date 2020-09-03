#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <tuple>
#include <utility>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { namespace CrosslinkSearch { class Crosslinker; } }

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    namespace CrosslinkSearch
    {

        // Having a std::tuple as the key of a std::unordered_map is not directly allowed
        // in C++.
        // See https://stackoverflow.com/questions/11408934/using-a-stdtuple-as-key-for-stdunordered-map
        // for how to handle this.
        
        typedef std::tuple<int, int> XLTuple;

        struct XLTuple_hash: public std::unary_function<XLTuple, std::size_t>{
            std::size_t operator() (const XLTuple& k ) const
            {
                size_t h1= std::hash<int>{}(std::get<0>(k));
                size_t h2= std::hash<int>{}(std::get<1>(k));
                return h1 ^ (h2 << 1);
            }
        };

        struct XLTuple_equal: public std::binary_function<XLTuple, XLTuple, bool>{
            bool operator() (const XLTuple& lhs, const XLTuple& rhs) const
            {
                return std::get<0>(lhs) == std::get<0>(rhs) &&
                    std::get<1>(lhs) == std::get<1>(rhs);
            }
        };
            
        typedef std::unordered_map<XLTuple, std::vector<Product *>, XLTuple_hash, XLTuple_equal> XLumap;
        
        class CrosslinkedPeptide
        {
        public:
            static std::vector<std::tuple<int, std::vector<Product*>>> XlGetTheoreticalFragments(DissociationType dissociationType, Crosslinker *crosslinker,
                                                                                                 std::vector<int> &possibleCrosslinkerPositions, double otherPeptideMass,
                                                                                                 PeptideWithSetModifications *peptide);
            
            static XLumap XlLoopGetTheoreticalFragments(DissociationType dissociationType, Modification *loopMass, std::vector<int> &modPos, PeptideWithSetModifications *peptide);
        };
    }
}
