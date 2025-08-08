#pragma once

#include <iostream>
#include <stdint.h>
#include <vector>
#include <cstring>
#include <cmath>
#include <array>

using PointId = size_t;
using NeighboursLookup = std::vector<std::vector<PointId>>;
using PointType = double;

template <size_t Dim>
struct Vector
{
    std::array<PointType, Dim> mValues;

    void Zero()
    {
        std::memset(mValues.data(), 0, sizeof(PointType)* Dim);
    }

    void Add(Vector<Dim> const & other)
    {
        for (size_t i = 0; i < Dim; i++)
        {
            mValues[i] += other.mValues[i];
        }
    }
};

#define FLOAT_POINTS
//#define INT_POINTS

#ifdef INT_POINTS
static constexpr PointType ScaledOne = 1L << 30;
// Can be relatively confident we won't overflow if we keep permutations relatively small before rescaling
static constexpr PointType ScaledOneSquared = 1L << 60;

PointType Divide(PointType a, PointType b)
{
    return a + b - 1 / b;
}

#elif defined(FLOAT_POINTS)
static constexpr PointType ScaledOne = 1;
static constexpr PointType ScaledOneSquared = 1;

PointType Divide(PointType a, PointType b)
{
    return a / b;
}

#else
#error must define one of FLOAT_POINTS or INT_POINTS
#endif


