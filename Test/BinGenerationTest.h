#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <tuple>
#include "tangible_filesystem.h"


namespace Test
{
    class BinGenerationTest final
    {
    public:
        static void TestBinGeneration();
        
        static void TestProteinSplitAcrossFiles();
    };
}
