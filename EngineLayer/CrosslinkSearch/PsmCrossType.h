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

                static PsmCrossType PsmCrossTypeFromString( std::string &s ) {
                    PsmCrossType t;
                    if ( s == "Single" ) {
                        t = PsmCrossType::Single;
                    }
                    else if ( s == "Cross" ) {
                        t =  PsmCrossType::Cross;
                    }
                    else if ( s == "DeadEnd" ) {
                        t = PsmCrossType::DeadEnd;
                    }
                    else if ( s == "Loop" ) {
                        t = PsmCrossType::Loop;
                    }
                    else if ( s == "Inter" ) {
                        t = PsmCrossType::Inter;
                    }
                    else if ( s == "Intra" ) {
                        t = PsmCrossType::Intra ;
                    }
                    else if ( s == "DeadEndH2O" ) {
                        t = PsmCrossType::DeadEndH2O ;
                    }
                    else if ( s == "DeadEndNH2" ) {
                        t = PsmCrossType::DeadEndNH2;
                    }
                    else if ( s == "DeadEndTris" ) {
                        t = PsmCrossType::DeadEndTris;
                    }
                    return t;
                }

                
        }
}
