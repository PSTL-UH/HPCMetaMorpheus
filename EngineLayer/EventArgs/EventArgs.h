#pragma once

namespace EngineLayer
{
    class EventArgs
    {
    public:
        virtual bool Equals( EventArgs *obj) const = 0;
        
        virtual int GetHashCode() const = 0;

        virtual std::string ToString()  const = 0;
    };
}
        
