#pragma once

#include <vector>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace TaskLayer { class DbForTask; }


namespace TaskLayer
{
	class XmlForTaskListEventArgs : public EventArgs
	{
	public:
		std::vector<DbForTask*> NewDatabases;

		XmlForTaskListEventArgs(std::vector<DbForTask*> &newDatabases);
	};
}
