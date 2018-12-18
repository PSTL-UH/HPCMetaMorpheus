#include "MyEngineTest.h"

using namespace EngineLayer;
using namespace NUnit::Framework;

namespace Test
{

	void MyEngineTest::TestMyEngine()
	{
		MetaMorpheusEngine *level0engine = new TestEngine(0);

		level0engine = new TestEngine(0);
		level0engine->Run();

		delete level0engine;
	}

	MyEngineTest::TestEngine::TestEngine(int level) : MetaMorpheusEngine()
	{
	}

	MetaMorpheusEngineResults *MyEngineTest::TestEngine::RunSpecific()
	{
		return new TestResults(this);
	}

	MyEngineTest::TestEngine::TestResults::TestResults(MetaMorpheusEngine *e) : MetaMorpheusEngineResults(e)
	{
	}

	std::wstring MyEngineTest::TestEngine::TestResults::ToString()
	{
		auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		sb->appendLine(MetaMorpheusEngineResults::ToString());
		sb->append(L"String for the TestResults results class");

		delete sb;
		return sb->toString();
	}
}
