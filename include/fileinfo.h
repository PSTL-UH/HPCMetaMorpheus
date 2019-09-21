#pragma once

#include <string>
#include <ctime>

class FileInfo
{
public:
    std::string DirectoryName;
    std::string Extension;
    std::string Name;
    std::string Fullname;

    bool Exists;
    bool IsReadOnly;

    long int Length;
    time_t CreationTime;
    time_t LastAccessTime;
    time_t LastWriteTime;
    
};
