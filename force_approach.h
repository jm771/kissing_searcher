
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
NeighboursLookup ConstructPointNeighbours(std::vector<Vector<Dim>> const & points, double margin)
{
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

template <size_t Dim>
void CalcRoundOfDiffs(std::vector<Vector<Dim>> const & points, NeighboursLookup const & neighbours, std::vector<Vector<Dim>> & rets)
{
    for (auto & val : rets)
    {
        val.Zero();
    }

    for (PointId pointId = 0; pointId < points.size(); pointId++)
    {
        auto const & point = points[pointId];
        auto & ret = rets[pointId];

        for (PointId neighbourId : neighbours[pointId])
        {
            auto const & neighbour = points[neighbourId];

            auto diff = Diff(point, neighbour);
            auto distsq = Dot(diff, diff);

            if (distsq < ScaledOneSquared)
            {
                auto & other = rets[neighbourId];

                int64_t dist = std::sqrt(distsq);

                static constexpr auto MinScale = 1 << 10;

                auto scaleFactor = (ScaledOne - dist + MinScale);

                DEBUG_LOG("%lu\n", scaleFactor);

                for (size_t i = 0; i < Dim; i++)
                {
                    // Scale quadratically to encourage even distn of points when over saturated
                    static constexpr int64_t MinQuad = ScaledOne >> 5;

                    auto dval = (((diff.mValues[i] * scaleFactor) + dist - 1) / dist);
                    if (scaleFactor > MinQuad)
                    {
                        dval *= scaleFactor;
                        dval /= ScaledOne;
                        dval += MinQuad;
                    }

                    ret.mValues[i] += dval;
                    other.mValues[i] -= dval;
                }
            }
        }
    }
}




template <size_t Dim>
bool AllZero(std::vector<Vector<Dim>> const & diffs)
{
    for (auto const & vec : diffs)
    {
        for (auto const & coeff : vec.mValues)
        {
            if (coeff != 0)
            {
                return false;
            }
        }
    }

    return true;
}

template <size_t Dim, typename Rand> 
void ApplyUnstick(std::vector<Vector<Dim>> & diffVectors, Rand & rand)
{
    // (void)diffVectors;
    // (void)rand;
    for (auto & vect : diffVectors)
    {
        int64_t size = std::sqrt(Dot(vect, vect));
        vect = RandPointOnBall<Dim>(10*size, rand);
    }
}

template <size_t Dim> 
void ApplyRecenter(std::vector<Vector<Dim>> & state)
{
    // Should probably not actually make this static...
    static Vector<Dim> working;
    working.Zero();

    for (auto const & vect : state)
    {
        working.Add(vect);
    }

    for (auto & vect : state)
    {
        for(size_t i = 0; i < Dim; i++)
        {
            vect.mValues[i] -= working.mValues[i] / state.size();
        }
    }    
}


template <size_t Dim> 
void ApplyDiffs(std::vector<Vector<Dim>> & state, std::vector<Vector<Dim>> & diffVects)
{
    // If we applied all pushes with a factor of 1/2 (as we apply the push to both balls)
    // we'd pop all the balls out to not be overlapping in a single iteration (Unless it's overlapping multiple balls). 
    // Let's slow this process a bit as it feels like it could be unstable if it shoots stuff straight into other balls
    static constexpr auto DiffScaleDivisor = 8 * 2;

    for (size_t i = 0; i < state.size(); i++)
    {
        for (size_t j = 0; j < Dim; j++)
        {
            state[i].mValues[j] += ((diffVects[i].mValues[j] + DiffScaleDivisor) / DiffScaleDivisor);
        }
    }
}

template <size_t Dim> 
void ApplyDiffsOrth(std::vector<Vector<Dim>> & state, std::vector<Vector<Dim>> & diffVects)
{
    // If we applied all pushes with a factor of 1/2 (as we apply the push to both balls)
    // we'd pop all the balls out to not be overlapping in a single iteration (Unless it's overlapping multiple balls). 
    // Let's slow this process a bit as it feels like it could be unstable if it shoots stuff straight into other balls
    static constexpr auto DiffScaleDivisor = 1 * 2;

    // We'll scale this down at the end to finish

    for (size_t i = 0; i < state.size(); i++)
    {
        auto & diffVect = diffVects[i];
        auto & stateVect = state[i];
        auto diffSize = std::sqrt(Dot(diffVect, diffVect));
        auto orthComp = Dot(diffVect, stateVect) / ScaledOne;
        
        for (size_t j = 0; j < Dim; j++)
        {
            diffVect.mValues[j] -= (stateVect.mValues[j] * orthComp) / ScaledOne;
        }

        // Adding 1 to denom to avoid issues with dividing by tiny numbers
        auto orthDiffSize = std::sqrt(Dot(diffVect, diffVect)) + 1;

        for (size_t j = 0; j < Dim; j++)
        {
            stateVect.mValues[j] += ((diffVect.mValues[j] + DiffScaleDivisor) * diffSize / DiffScaleDivisor / orthDiffSize);
        }
    }
}

template <size_t Dim, typename Rand> 
bool RunRoutine(std::vector<Vector<Dim>> & initialState, Rand & rand)
{
    FileOutput frameOutput("viewer/frames.json");
    auto & state = initialState;
    static constexpr size_t OuterEpochs = 100 * 1000;
    static constexpr size_t UnstickCadence = 1000;
    static constexpr size_t RecenterCadence = 1000;
    static constexpr size_t InnerIterationLoops = 1;


    std::vector<Vector<Dim>> diffVect(state.size());
    for (auto & vec : diffVect)
    {
        vec.Zero();
    }

    for (size_t outerEpoch = 0; outerEpoch < OuterEpochs; outerEpoch++)
    {
        DEBUG_LOG("[OuterEpoch=%lu] Before Normalization:", outerEpoch);
        DEBUG_LOG_VECTS(state);

        // static constexpr int64_t MIN_SCALE = (7 * ScaledOne) / 8;
        // static constexpr int64_t END_SCALE_EPOCH = (7 * OuterEpochs) / 8;

        // int64_t scale = outerEpoch > END_SCALE_EPOCH ? ScaledOne : ((END_SCALE_EPOCH - outerEpoch) * MIN_SCALE + outerEpoch * ScaledOne) / END_SCALE_EPOCH;

        int64_t scale = ScaledOne;

        Normalize(state, scale);
        DEBUG_LOG("[OuterEpoch=%lu] After Normalization:", outerEpoch);
        DEBUG_LOG_VECTS(state);

        if (outerEpoch % RecenterCadence == 0)
        {
            ApplyRecenter(diffVect);
            Normalize(state, scale);
        }

        auto neighbourLookup = ConstructPointNeighbours(state, ScaledBound(1.2));
        // DEBUG_LOG_LOOKUP(neighbourLookup);

        for (size_t innerEpoch = 0; innerEpoch < InnerIterationLoops; innerEpoch++)
        {
            frameOutput.WriteRow(state);
            CalcRoundOfDiffs(state, neighbourLookup, diffVect);
            DEBUG_LOG("[OuterEpoch=%lu][InnerEpoch=%lu]Diffs pre-scaling:", outerEpoch, innerEpoch);
            DEBUG_LOG_VECTS(diffVect);

            if (innerEpoch == 0) {
                if (AllZero(diffVect)) {
                    printf("[OuterEpoch=%lu][InnerEpoch=%lu] No diffs in iteration, exiting, final values:\n", outerEpoch, innerEpoch);
                    PrintVects(state);
                    return true;
                }
            }

            if (outerEpoch % UnstickCadence == 0)
            {
                ApplyUnstick(diffVect, rand);
            }
            
	    ApplyDiffs(state, diffVect);
        }
    }

    
    printf("No Soltuion found after %lu iterations:\n", OuterEpochs);
    PrintVects(state);

    return false;
}