#pragma once

#include <string>
#include <iostream>
#include <any>
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class StringEventArgs; }

// Copyright 2016 Stefan Solntsev
using namespace EngineLayer;
using namespace NUnit::Framework;
using namespace TaskLayer;
using namespace UsefulProteomicsDatabases;

namespace Test
{
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [SetUpFixture] public class MySetUpClass
	class MySetUpClass
	{
	public:
		static std::wstring outputFolder;

	private:
		static const std::wstring elementsLocation;

	public:
//C# TO C++ CONVERTER NOTE: The following .NET attribute has no direct equivalent in native C++:
//ORIGINAL LINE: [OneTimeSetUp] public static void Setup()
		static void Setup();

	private:
		static void SuccessfullyFinishedAllTasks(std::any sender, StringEventArgs *rootOutputFolderPath);

		static void WarnStatusHandler(std::any sender, StringEventArgs *e);
	};
}
