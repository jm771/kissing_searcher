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
    std::array<int64_t, Dim> mValues;

    void Zero()
    {
        std::memset(mValues.data(), 0, sizeof(int64_t)* Dim);
    }

    void Add(Vector<Dim> const & other)
    {
        for (size_t i = 0; i < Dim; i++)
        {
            mValues[i] += other.mValues[i];
        }
    }
};

static constexpr int64_t ScaledOne = 1L << 30;
// Can be relatively confident we won't overflow if we keep permutations relatively small before rescaling
static constexpr int64_t ScaledOneSquared = 1L << 60;