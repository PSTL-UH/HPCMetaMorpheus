#pragma once

#include <string>

namespace TaskLayer
{
	enum class MassDiffAcceptorType
	{
		Exact,
		OneMM,
		TwoMM,
		ThreeMM,
		ModOpen,
		Open,
		Custom
	};

        static MassDiffAcceptorType  MassDiffAcceptorTypeFromString(std::string s) {
            MassDiffAcceptorType t;

            if ( s == "Exact" ) {
                t = MassDiffAcceptorType::Exact;
            }
            else if ( s == "OneMM" ) {
                t = MassDiffAcceptorType::OneMM;
            }
            else if ( s == "TwoMM" ) {
                t = MassDiffAcceptorType::TwoMM;                
            }
            else if ( s ==  "ThreMM" ) {
                t =  MassDiffAcceptorType::ThreeMM;
            }
            else if ( s == "ModOpen" ) {
                t = MassDiffAcceptorType::ModOpen;
            }
            else if ( s == "Open" ) {
                t = MassDiffAcceptorType::Open;
            }
            else if ( s ==  "Custom" ) {
                t = MassDiffAcceptorType::Custom;
            }
            
            return t;
        }
        
        static std::string  MassDiffAcceptorTypeToString( MassDiffAcceptorType t ) {
            std::string s;

            if ( t == MassDiffAcceptorType::Exact ) {
                s = "Exact";
            }
            else if ( t == MassDiffAcceptorType::OneMM ) {
                s = "OneMM";
            }
            else if ( t == MassDiffAcceptorType::TwoMM ) {
                s = "TwoMM";
            }
            else if ( t == MassDiffAcceptorType::ThreeMM ) {
                s = "ThreMM";
            }
            else if ( t == MassDiffAcceptorType::ModOpen ) {
                s = "ModOpen";
            }
            else if ( t == MassDiffAcceptorType::Open ) {
                s = "Open";
            }
            else if ( t == MassDiffAcceptorType::Custom ) {
                s = "Custom";
            }
            
            return s;
        }
        
}
