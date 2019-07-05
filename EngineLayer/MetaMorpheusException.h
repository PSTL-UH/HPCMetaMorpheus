#pragma once

#include <string>
#include <stdexcept>


namespace EngineLayer
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [Serializable] public class MetaMorpheusException : Exception
	class MetaMorpheusException : public std::runtime_error
	{
	public:
		MetaMorpheusException(const std::string &message);
	};
}
