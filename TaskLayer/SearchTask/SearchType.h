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

        static SearchType SearchTypeFromString ( std::string s ) {
            SearchType t;

            if  ( s == "Classic" ) {
                t = SearchType::Classic;
            }
            else if ( s == "Modern" ) {
                t = SearchType::Modern;
            }
            else if ( s == "NonSpecific" ) {
                t = SearchType::NonSpecific;
            }
            return t;
        }
        
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
