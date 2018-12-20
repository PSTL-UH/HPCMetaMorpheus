#pragma once

//----------------------------------------------------------------------------------------
//	Copyright © 2004 - 2018 Tangible Software Solutions, Inc.
//	This class can be used by anyone provided that the copyright notice remains intact.
//
//	This class is used to replace some static .NET file and folder method calls
//	with std::filesystem method calls.
//----------------------------------------------------------------------------------------
#include <string>
#include <experimental/filesystem>

class FileSystem
{
public:
	static void createDirectory(const std::wstring &path)
	{
		std::filesystem::create_directory(pathFromString(path));
	}

	static bool fileExists(const std::wstring &path)
	{
		return std::filesystem::is_regular_file(pathFromString(path));
	}

	static bool directoryExists(const std::wstring &path)
	{
		return std::filesystem::is_directory(pathFromString(path));
	}

	static std::wstring combine(const std::wstring &path1, const std::wstring &path2)
	{
		return (pathFromString(path1) / pathFromString(path2)).generic_wstring();
	}

	static bool isPathRooted(const std::wstring &path)
	{
		return pathFromString(path).has_root_path();
	}

	static std::wstring getFullPath(const std::wstring &path)
	{
		return std::filesystem::absolute(pathFromString(path)).generic_wstring();
	}

	static std::wstring getFileName(const std::wstring &path)
	{
		return std::filesystem::path(pathFromString(path)).filename().generic_wstring();
	}

	static std::wstring getDirectoryName(const std::wstring &path)
	{
		return std::filesystem::path(pathFromString(path)).parent_path().generic_wstring();
	}

	static std::wstring getCurrentDirectory()
	{
		return std::filesystem::current_path().generic_wstring();
	}

	static void copyFile(const std::wstring &path1, const std::wstring &path2)
	{
		std::filesystem::copy_file(pathFromString(path1), pathFromString(path2));
	}

	static void renamePath(const std::wstring &path1, const std::wstring &path2)
	{
		std::filesystem::rename(pathFromString(path1), pathFromString(path2));
	}

	static wchar_t preferredSeparator()
	{
		return std::filesystem::path::preferred_separator;
	}

private:
	static std::filesystem::path pathFromString(const std::wstring &path)
	{
		return std::filesystem::path(&path[0]);
	}
};
