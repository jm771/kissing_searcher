#pragma once

#include "types.h"
#include "vectors.h"
#include <random>
#include <numbers>

// template <size_t Dim> 
// Vector<Dim> UniformOnBall(std::random)
// {
//     std::default_random_generator(1234);

// }

template <size_t Dim, typename Rand>
Vector<Dim> RandPoint(PointType stddev, Rand & rand)
{
    std::normal_distribution<double> gauss(0, 1.0);
    Vector<Dim> ret;
    for (size_t i = 0; i < Dim; i++)
    {
        ret.mValues[i] = static_cast<PointType>(gauss(rand) * stddev);
    }

    return ret;
}

template <size_t Dim, typename Rand>
Vector<Dim> RandPointOnBall(PointType radius, Rand & rand)
{
    auto ret = RandPoint<Dim>(radius / 16, rand);

    Normalize(ret, radius);
    return ret;
}

template <size_t Dim, typename Rand>
std::vector<Vector<Dim>> Initialize(size_t nBalls, PointType radius, Rand & rand)
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
    std::vector<PointType> signs{-ScaledOne, ScaledOne};
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
                    n.Add(RandPointOnBall<4>(ScaledOne / 2, rand));
                }
            }
    
        }
    }

    return ret;
}

double GetCoord(size_t i, size_t stratInLayer)
{
    return ((i * ScaledOne * 2) / (stratInLayer - 1)) - ScaledOne;
}

template <size_t Dim>
std::vector<Vector<Dim>> StratifyPointsCubic(size_t nToStratefy)
{
    const size_t stratInLayer = std::ceil(std::pow(nToStratefy, 1.0/Dim));
    std::vector<Vector<Dim>> ret{};
    if constexpr (Dim == 1) {
        for (size_t i = 0; i < stratInLayer; i++)
        {
            ret.emplace_back(Vector<1>{GetCoord(i, stratInLayer)});
        }
    }
    else {
        for (size_t i = 0; i < stratInLayer; i++)
        {
            auto inner = StratifyPointsCubic<Dim - 1>((nToStratefy + i) / stratInLayer);
            for (auto const & el : inner) {
                auto & newEl = ret.emplace_back();
                std::memcpy(newEl.mValues.data(), el.mValues.data(), sizeof(PointType) * (Dim - 1));
                newEl.mValues[Dim - 1] = GetCoord(i, stratInLayer);
            }
        }
    }

    return ret;
}

// template <size_t Dim, typename Rand>
// std::vector<Vector<Dim>> InitializeStratified(size_t nBalls)
// {
//     if constexpr (Dim == 2) {
//         std::vector<Vector<2>> ret();
//         for (size_t i = 0; i < nBalls; i++)
//         {
//             auto rads = std::pi_v * 2 * i / nBalls;
//             ret.emplace_back({std::sin(rads) * ScaledOne / 2, d::cos(rads) * ScaledOne / 2});
//         }
//         return ret;
//     }
//     else
//     {
//         std::vector<Vector<Dim>> ret();

//     for (size_t i = 0; i < nBalls; i++)
//     {
//         auto sign = ((i % 2) * 2 - 1;
//         auto axis = (i / 2) % Dim;
//         auto
        
//     }

//     }
// }

// template <size_t Dim, typename Rand>
// std::vector<Vector<Dim>> InitializeStratefied(size_t nBalls, Rand & rand)
// {
    

//     return ret;
// }




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
