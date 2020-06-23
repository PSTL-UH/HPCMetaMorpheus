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
