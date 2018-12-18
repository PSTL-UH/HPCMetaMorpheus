#pragma once

#include <exception>

class OutOfMemoryException : public std::exception
{
private:
    std::string msg;

public:
    OutOfMemoryException(const std::string& message = "") : msg(message)
    {
    }

    const char * what() const throw()
    {
        return msg.c_str();
    }
};

class NotImplementedException : public std::exception
{
private:
    std::string msg;

public:
    NotImplementedException(const std::string& message = "") : msg(message)
    {
    }

    const char * what() const throw()
    {
        return msg.c_str();
    }
};

class MzLibException : public std::exception
{
private:
    std::string msg;

public:
    MzLibException(const std::string& message = "") : msg(message)
    {
    }

    const char * what() const throw()
    {
        return msg.c_str();
    }
};
