#pragma once

#include <string>
#include <vector>

#include "EventArgs.h"

namespace EngineLayer
{
	class StringListEventArgs : public EventArgs
	{
	private:
            std::vector<std::string> privateStringList;

	public:
            StringListEventArgs(std::vector<std::string> &stringList);
                
            std::vector<std::string> getStringList() const;

            bool Equals( EventArgs *obj) const override;
    
            int GetHashCode() const override;
            
            std::string ToString() const override;
        };
}
