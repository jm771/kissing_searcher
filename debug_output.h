#pragma once

#include "types.h"
#include <format>

template <size_t Dim>
void PrintVect(Vector<Dim> const & vec, std::ostream & out)
{
    out << '[';
    bool first = true;
    for (size_t i = 0; i < Dim; i++)
    {
        auto const & coord = vec.mValues[i];
        if (first) { first = false;} else { out << ", ";}
        out << coord;
    }
    out << ']';
}

template <size_t Dim>
void PrintVects(std::vector<Vector<Dim>> const & vects)
{
    for (auto const & vec : vects)
    {
        PrintVect(vec, std::cout);
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void PrintLookups(NeighboursLookup const & lookup)
{
    for (size_t i = 0; i < lookup.size(); i++)
    {
        std::cout << i << ": [";
        for (auto const & val : lookup[i])
        {
            std::cout << val << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << std::endl;
}


// #define DEBUG_LOGGING

#ifdef DEBUG_LOGGING
#define DEBUG_LOG(format, ...) \
do {\
    printf(format, ##__VA_ARGS__); \
} while (false)
# else 
#define DEBUG_LOG(fmt, ...)
#endif

#ifdef DEBUG_LOGGING
#define DEBUG_LOG_VECTS(vects) PrintVects(vects)
#else
#define DEBUG_LOG_VECTS(vects)
#endif

#ifdef DEBUG_LOGGING
#define DEBUG_LOG_LOOKUP(lookup) PrintLookups(lookup)
#else
#define DEBUG_LOG_LOOKUP(lookup)
#endif

struct Exception : std::exception {
    Exception(std::string const & what) : mWhat(what) {

    }

    const char* what() const noexcept override
    {
        return mWhat.c_str();
    }
    

    std::string mWhat; 
};

#define ASSERT_MSG(x, msg, ...) do { \
     if (!(x)) { throw Exception(std::format("Assertion at {}:{}\nmsg: {}", __FILE__, __LINE__, std::format(msg, ##__VA_ARGS__))); } \
} while(false)

#define ASSERT(x) ASSERT_MSG(x, "");