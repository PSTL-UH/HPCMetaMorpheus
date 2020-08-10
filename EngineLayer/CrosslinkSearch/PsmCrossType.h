#pragma once

#include <string>

namespace EngineLayer
{
	namespace CrosslinkSearch
	{
		enum class PsmCrossType
		{
			Single,
			Cross,
			DeadEnd,
			Loop,
			Inter,
			Intra,
			DeadEndH2O,
			DeadEndNH2,
			DeadEndTris
		};

                static std::string PsmCrossTypeToString( PsmCrossType &t ) {
                    std::string s;
                    if ( t == PsmCrossType::Single ) {
                        s = "Single";
                    }
                    else if ( t ==  PsmCrossType::Cross) {
                        s = "Cross";
                    }
                    else if ( t == PsmCrossType::DeadEnd) {
                        s = "DeadEnd";
                    }
                    else if ( t == PsmCrossType::Loop) {
                        s = "Loop";
                    }
                    else if ( t == PsmCrossType::Inter) {
                        s = "Inter";
                    }
                    else if ( t == PsmCrossType::Intra ) {
                        s = "Intra";
                    }
                    else if ( t == PsmCrossType::DeadEndH2O ) {
                        s = "DeadEndH2O";
                    }
                    else if ( t == PsmCrossType::DeadEndNH2) {
                        s = "DeadEndNH2";
                    }
                    else if ( t == PsmCrossType::DeadEndTris) {
                        s = "DeadEndTris";
                    }
                    return s;
                }

        }
}
