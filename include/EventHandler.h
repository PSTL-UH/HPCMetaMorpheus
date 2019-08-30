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

#include "EventArgs.h"

template<typename T>
class EventHandler
{
private:
    std::unordered_map<std::string, T> namedListeners;

public:
    void addListener(const std::string &methodName, T namedEventHandlerMethod)
    {
        if (namedListeners.find(methodName) == namedListeners.end())
            namedListeners[methodName] = namedEventHandlerMethod;
    }
    
    void removeListener(const std::string &methodName)
    {
        if (namedListeners.find(methodName) != namedListeners.end())
            namedListeners.erase(methodName);
    }
    
private:
    std::vector<T> anonymousListeners;
public:
    void addListener(T unnamedEventHandlerMethod)
    {
        anonymousListeners.push_back(unnamedEventHandlerMethod);
    }
    
    std::vector<T> listeners()
    {
        std::vector<T> allListeners;
        for (auto listener : namedListeners)
        {
            allListeners.push_back(listener.second);
        }
        allListeners.insert(allListeners.end(), anonymousListeners.begin(), anonymousListeners.end());
        return allListeners;
    }

public:
    void Invoke (void *vent, EventArgs *args);
};
