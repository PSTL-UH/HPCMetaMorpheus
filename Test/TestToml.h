#pragma once

#include <string>
#include <vector>
#include "tangible_filesystem.h"


#include "../EngineLayer/EngineLayer.h"
using namespace EngineLayer;

#include "MzLibUtil.h"
using namespace MzLibUtil;

#include "../TaskLayer/TaskLayer.h"
using namespace TaskLayer;

namespace Test
{
	class TestToml final
	{
	public:
		static void TestTomlFunction();

		static void TestTomlForSpecficFiles();
	};
}
