#pragma once

// This file is an extended version of code provided by The Tangible Software Solution C# to
// C++ convertor
//----------------------------------------------------------------------------------------
//	Copyright © 2004 - 2018 Tangible Software Solutions, Inc.
//	This class can be used by anyone provided that the copyright notice remains intact.
//
//----------------------------------------------------------------------------------------
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>

template<typename T>
class EventHandler
{
private:
    std::unordered_map<std::string, std::function<void(T)>> namedListeners;

public:
    void addListener(const std::string &methodName, std::function<void(T)> namedEventHandlerMethod)
    {
        if (namedListeners.find(methodName) == namedListeners.end())
            namedListeners[methodName] = namedEventHandlerMethod;
    }
    
    void removeListener(const std::string &methodName)
    {
        if (namedListeners.find(methodName) != namedListeners.end())
            namedListeners.erase(methodName);
    }
    
    void Invoke (T args) {
        for ( auto p= namedListeners.begin(); p!= namedListeners.end(); p++ ) {
            p->second(args);
        }
    }
};
