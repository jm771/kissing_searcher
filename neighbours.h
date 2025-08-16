#pragma once

#include "types.h"
#include "debug_output.h"
#include "file_output.h"
#include "initial_states.h"
#include "vectors.h"
#include <stdint.h>
#include <random>

template <size_t Dim>
bool CloserThanSafe(Vector<Dim> const & a, Vector<Dim> const & b, double bound)
{
    auto diffVec = Diff(a, b);
    return DotSafe(diffVec, diffVec) <= bound;
}


template <size_t Dim>
NeighboursLookup ConstructPointNeighbours(std::vector<Vector<Dim>> const & points)
{
    static constexpr PointType margin = 1.2;

    std::vector<std::vector<PointId>> ret;
    for (PointId pointId = 0; pointId < points.size(); pointId++)
    {
        auto const & point = points[pointId];
        auto & neighbours = ret.emplace_back();
        for (PointId maybeNeighbourId = pointId+1; maybeNeighbourId < points.size(); maybeNeighbourId++)
        {
            if (CloserThanSafe(point, points[maybeNeighbourId], margin)) {
                neighbours.push_back(maybeNeighbourId);
            }
        }
    }

    return ret;
}

template <size_t Dim>
NeighboursLookup ConstructPointNeighboursBidi(std::vector<Vector<Dim>> const & points, double margin)
{
    std::vector<std::vector<PointId>> ret;
    for (PointId pointId = 0; pointId < points.size(); pointId++)
    {
        auto const & point = points[pointId];
        auto & neighbours = ret.emplace_back();
        for (PointId maybeNeighbourId = 0; maybeNeighbourId < points.size(); maybeNeighbourId++)
        {
            if (maybeNeighbourId == pointId)
            {
                continue;
            }
            
            if (CloserThanSafe(point, points[maybeNeighbourId], margin)) {
                neighbours.push_back(maybeNeighbourId);
            }
        }
    }

    return ret;
}