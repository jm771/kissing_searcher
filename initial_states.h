#include "types.h"


std::vector<Vector<Dim>> Initialize(size_t nBalls)
{
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(1.0, 10.0)
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
