#include <iostream>
#include <experimental/filesystem>
#include "TomlReadFile.h"


toml::Value* TomlReadFile::tomlReadFile(std::string FilePath, std::string ValueName)
{
	//Open and parse toml file located at 'FilePath'
	std::ifstream ifs(FilePath);
	toml::ParseResult parsed_toml = toml::parse(ifs);

	// check if parsed toml is valid
	if (!parsed_toml.valid()) {
		std::cout << parsed_toml.errorReason << std::endl;
		//return -1;
	}

	// parsed_toml.value is the parsed value.  The parsed value is used when 
	//looking up the 'ValueName'
	const toml::Value& toml_value = parsed_toml.value;
	const toml::Value* found_toml_value = toml_value.find(ValueName);
	toml::Value* value_to_return = const_cast<toml::Value *>(found_toml_value);

	return value_to_return;
}
