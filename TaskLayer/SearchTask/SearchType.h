#pragma once

#include <string>

namespace TaskLayer
{
	enum class SearchType
	{
		Classic,
		Modern,
		NonSpecific
	};

        static std::string SearchTypeToString( SearchType t ) {
            std::string s;

            if ( t == SearchType::Classic ) {
                s = "Classic";
            }
            else if ( t == SearchType::Modern ) {
                s = "Modern";
            }
            else if ( t == SearchType::NonSpecific ) {
                s = "NonSpecific";
            }
            
            return s;
        }
}
