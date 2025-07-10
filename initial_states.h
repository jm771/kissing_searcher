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

    Normalize(ret, radius);
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

template <typename Rand>
std::vector<Vector<4>> Initialize4D(Rand & rand)
{
    std::vector<Vector<4>> ret;
    std::vector<int64_t> signs{-ScaledOne, ScaledOne};
    for (int i = 0; i < 4; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            for (auto s1 : signs)
            {
                for (auto s2 : signs)
                {
                    auto & n = ret.emplace_back();
                    n.Zero();
                    n.mValues[i] = s1;
                    n.mValues[j] = s2;
                    n.Add(RandPointOnBall<4>(ScaledOne >> 10, rand));
                }
            }
    
        }
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
