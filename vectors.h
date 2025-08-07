#pragma once

#include "types.h"

template <size_t Dim>
int64_t Dot (Vector<Dim> const & a, Vector<Dim> const & b)
{
    int64_t ret = 0;
    for (size_t i = 0; i < Dim; i++)
    {
        ret += a.mValues[i] * b.mValues[i];
    }

    return ret;
}

template <size_t Dim>
double DotSafe (Vector<Dim> const & a, Vector<Dim> const & b)
{
    double ret = 0;
    for (size_t i = 0; i < Dim; i++)
    {
        ret += static_cast<double>(a.mValues[i]) * b.mValues[i];
    }

    return ret;
}

template <size_t Dim>
Vector<Dim> Diff(Vector<Dim> const & a, Vector<Dim> const & b)
{
    Vector<Dim> result;

    for (size_t i = 0; i < Dim; i++)
    {
        result.mValues[i] = a.mValues[i] - b.mValues[i];
    }

    return result;
}

template <size_t Dim>
Vector<Dim> Add(Vector<Dim> const & a, Vector<Dim> const & b)
{
    Vector<Dim> result;

    for (size_t i = 0; i < Dim; i++)
    {
        result.mValues[i] = a.mValues[i] + b.mValues[i];
    }

    return result;
}

template <size_t Dim>
double Dist(Vector<Dim> const & a, Vector<Dim> const & b)
{
    auto diffVec = Diff(a, b);
    return std::sqrt(DotSafe(diffVec, diffVec));
}

template <size_t Dim>
void Normalize(Vector<Dim> & point, int64_t mag)
{
    auto squareMag = Dot(point, point);
    auto divisor = std::sqrt(squareMag);
    for (auto & coord : point.mValues)
    {
        // Round up I guess? Todo think more
        coord = ((coord * mag) + divisor - 1)/ divisor;
    }
}

template <size_t Dim>
void Normalize(std::vector<Vector<Dim>> & points, int64_t mag)
{
    for (auto & point : points)
    {
        Normalize(point, mag);
    }
}

static consteval double ScaledBound(double val)
{
    ASSERT(val < 2);
    return ScaledOneSquared * val * val;
}
