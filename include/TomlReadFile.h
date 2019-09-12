#pragma once
#include "toml.h"

class TomlReadFile
{
public:
	
	toml::Value* tomlReadFile(std::string FilePath, std::string ValueName);
};

