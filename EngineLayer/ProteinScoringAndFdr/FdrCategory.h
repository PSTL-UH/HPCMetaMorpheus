#pragma once

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
    enum class FdrCategory
    {
        //Cleavage Specificity
        FullySpecific = 0,
        SemiSpecific = 1,
        NonSpecific = 2
            
        //New category here
   };

    static std::string FdrCategoryToString( FdrCategory t ) {
        std::string s ;
        if ( t ==  FdrCategory::FullySpecific ) {
            s = "FullySpecific";
        }
        else if ( t ==  FdrCategory::SemiSpecific ) {
            s = "SemiSpecific";
        }
        else if ( t ==  FdrCategory::NonSpecific ) {
            s = "NonSpecific";
        }

        return s;
    }
}
