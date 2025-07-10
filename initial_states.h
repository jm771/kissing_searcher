#pragma once

#include "types.h"
#include "vectors.h"
#include <random>

// template <size_t Dim> 
// Vector<Dim> UniformOnBall(std::random)
// {
//     std::default_random_generator(1234);

// }

template <size_t Dim, typename Rand>
Vector<Dim> RandPointOnBall(int64_t radius, Rand & rand)
{
    std::normal_distribution<double> gauss(0, 1.0 / 16.0);
    Vector<Dim> ret;
    for (size_t i = 0; i < Dim; i++)
    {
        ret.mValues[i] = static_cast<int64_t>(gauss(rand) * radius);
    }

    Normalize(ret, ScaledOne);
    return ret;
}

template <size_t Dim, typename Rand>
std::vector<Vector<Dim>> Initialize(size_t nBalls, int64_t radius, Rand & rand)
{
    std::vector<Vector<Dim>> ret;
    for (size_t i = 0; i < nBalls; i++)
    {
        ret.push_back(RandPointOnBall<Dim>(radius, rand));
    }
    
    return ret;
}


std::vector<Vector<2>> Initialize2D()
{
    return {
        {0, ScaledOne} ,
        // Extra one
        { {ScaledOne/100, ScaledOne} },
        // Extra two - unclear if this will converge now
        { {-ScaledOne/100, ScaledOne} },
        {ScaledOne, 0},
        {0, -ScaledOne},
        { {-ScaledOne, 0} }
    };
}
