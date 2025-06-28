#pragma once

#include <iostream>
#include <stdint.h>
#include <vector>
#include <cstring>
#include <cmath>
#include <array>

using PointId = size_t;
using NeighboursLookup = std::vector<std::vector<PointId>>;

template <size_t Dim>
struct Vector
{
    std::array<int64_t, Dim> mValues;

    void Zero()
    {
        std::memset(mValues.data(), 0, sizeof(int64_t)* Dim);
    }
};