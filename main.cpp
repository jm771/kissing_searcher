
#include "types.h"
#include "debug_output.h"
#include "file_output.h"
#include "initial_states.h"
#include "vectors.h"
#include <stdint.h>
#include <random>

// static constexpr int64_t GoogleScaledOne = 1e13;


// // TODO - commit to fixed point and do fixed point upper / lower ops with the __128 stuff

// static constexpr int64_t RescaleGoogleValue(int64_t googleValue)
// {
//     // 1e13 < 2^44 - so that leaves us 2^19 beore 2^63
//     // So we'll shift by 19 before div and 11 after to make the most of our 64 bit precision
//     return ((googleValue << 19) / GoogleScaledOne) << 11;
// }

static consteval double ScaledBound(double val)
{
    ASSERT(val < 2);
    return ScaledOneSquared * val * val;
}



// template <size_t Dim>
// int64_t DistSquared (Vector<Dim> const & a, Vector<Dim> const & b)
// {
//     int64_t ret = 0;
//     for (size_t i = 0; i < Dim; i++)
//     {
//         int64_t diff = a.mValues[i] - b.mValues[i];
//         ret += diff * diff;
//     }
        
//     return ret;
// }



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
                    static constexpr int64_t MinQuad = ScaledOne >> 10;

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
    (void)diffVectors;
    (void)rand;
    // for (auto & vect : diffVectors)
    // {
    //     int64_t size = std::sqrt(Dot(vect, vect));
    //     vect = RandPointOnBall<Dim>(10*size, rand);
    // }
}

template <size_t Dim> 
void ApplyDiffs(std::vector<Vector<Dim>> & state, std::vector<Vector<Dim>> & diffVects)
{
    // If we applied all pushes with a factor of 1/2 (as we apply the push to both balls)
    // we'd pop all the balls out to not be overlapping in a single iteration (Unless it's overlapping multiple balls). 
    // Let's slow this process a bit as it feels like it could be unstable if it shoots stuff straight into other balls
    static constexpr auto DiffScaleDivisor = 1 * 2;

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

        auto orthDiffSize = std::sqrt(Dot(diffVect, diffVect));

        for (size_t j = 0; j < Dim; j++)
        {
            stateVect.mValues[j] += ((diffVect.mValues[j] + DiffScaleDivisor) * diffSize / DiffScaleDivisor / orthDiffSize);
        }
    }
}

template <size_t Dim, typename Rand> 
bool RunRoutine(std::vector<Vector<Dim>> & initialState, Rand & rand)
{
    FileOutput frameOutput("../kissing_search_viewer/frames.json");
    auto & state = initialState;
    static constexpr size_t OuterEpochs = 100 * 1000;
    static constexpr size_t UnstickCadance = 1000;
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
        Normalize(state);
        DEBUG_LOG("[OuterEpoch=%lu] After Normalization:", outerEpoch);
        DEBUG_LOG_VECTS(state);

        auto neighbourLookup = ConstructPointNeighbours(state, ScaledBound(1.2));
        DEBUG_LOG_LOOKUP(neighbourLookup);

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

            if (outerEpoch % UnstickCadance == 0)
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

template <size_t Dim> 
bool Validate(std::vector<Vector<Dim>> & result)
{
    int64_t maxSize = 0;
    DEBUG_LOG("Square Sizes:\n");
    for (auto const & vect : result)
    {
        auto size = Dot(vect, vect);
        DEBUG_LOG("%lu\n", size);
        maxSize = std::max(maxSize, size);
    }

    DEBUG_LOG("\n");

    int64_t minDist = std::numeric_limits<int64_t>::max();

    for (size_t i = 0; i < result.size(); i++)
    {
        for (size_t j = 0; j < result.size(); j++)
        {
            if (i != j)
            {
                auto diff = Diff(result[i], result[j]);
                auto dist = Dot(diff, diff);
                minDist = std::min(minDist, dist);

                // Crude check to avoid overflow
                ASSERT(dist < 5 * ScaledOneSquared);
                DEBUG_LOG("i=%lu j=%lu dist=%lu\n", i, j, dist);
            }
        }
    }

    std::cout << "Max size (squared): " << maxSize << std::endl;

    std::cout << "Min dist (squared): " << minDist << std::endl;

    bool valid = maxSize < minDist;

    auto message =  valid ? "Passed" : "Failed";

    std::cout << "Validation " << message << std::endl;

    return valid;
}

int main(int, char**){
    // auto init = Initialize2D();
    // static constexpr size_t DIMENSION = 2; static constexpr size_t targetBalls = 6;
    static constexpr size_t DIMENSION = 3; static constexpr size_t targetBalls = 12;
    // static constexpr size_t DIMENSION = 4; static constexpr size_t targetBalls = 24;


    std::mt19937 rand(12345);


    auto state = Initialize<DIMENSION>(targetBalls, ScaledOne, rand);
    auto success = RunRoutine<DIMENSION>(state, rand);

    if (success) {
        Validate(state);
    }
    return 0;
}
