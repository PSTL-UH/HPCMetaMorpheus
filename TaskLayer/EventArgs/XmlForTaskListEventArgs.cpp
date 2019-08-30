#include "XmlForTaskListEventArgs.h"
#include "../DbForTask.h"


namespace TaskLayer
{

	XmlForTaskListEventArgs::XmlForTaskListEventArgs(std::vector<DbForTask*> &newDatabases)
	{
		NewDatabases = newDatabases;
	}
}
