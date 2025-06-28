
#include "types.h"
#include "debug_output.h"
#include "file_output.h"
#include <stdint.h>
#include <random>


// Maybe we should make this a power of 2;
static constexpr int64_t GoogleScaledOne = 1e13;
static constexpr int64_t ScaledOne = 1L << 30;
// Can be relatively confident we won't overflow if we keep permutations relatively small before rescaling
static constexpr int64_t ScaledOneSquared = 1L << 60;


// TODO - commit to fixed point and do fixed point upper / lower ops with the __128 stuff

static constexpr int64_t RescaleGoogleValue(int64_t googleValue)
{
    // 1e13 < 2^44 - so that leaves us 2^19 beore 2^63
    // So we'll shift by 19 before div and 11 after to make the most of our 64 bit precision
    return ((googleValue << 19) / GoogleScaledOne) << 11;
}

static consteval double ScaledBound(double val)
{
    ASSERT(val < 2);
    return ScaledOneSquared * val * val;
}

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
int64_t DistSquared (Vector<Dim> const & a, Vector<Dim> const & b)
{
    int64_t ret = 0;
    for (size_t i = 0; i < Dim; i++)
    {
        int64_t diff = a.mValues[i] - b.mValues[i];
        ret += diff * diff;
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
bool CloserThanSafe (Vector<Dim> const & a, Vector<Dim> const & b, double bound)
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
                auto scaleFactor = (ScaledOneSquared - distsq) / ScaledOne;

                std::cout << scaleFactor << std::endl;

                for (size_t i = 0; i < Dim; i++)
                {
                    auto const dval = (diff.mValues[i] * scaleFactor) / ScaledOne;
                    ret.mValues[i] += dval;
                    other.mValues[i] -= dval;
                }
            }
        }
    }
}

template <size_t Dim>
void Renormalize(std::vector<Vector<Dim>> & points)
{
    for (auto & point : points)
    {
        auto squareMag = Dot(point, point);
        auto divisor = std::sqrt(squareMag);
        for (auto & coord : point.mValues)
        {
            // Round up I guess? Todo think more
            coord = ((coord * ScaledOne) + divisor - 1)/ divisor;
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


template <size_t Dim> 
void RunRoutine(std::vector<Vector<Dim>> & initialState)
{
    FileOutput frameOutput("../kissing_search_viewer/frames.json");
    auto & state = initialState;
    static constexpr auto OuterEpochs = 1 * 1000;
    static constexpr auto InnerIterationLoops = 1;
    // If we applied all pushes with a factor of 1/2 (as we apply the push to both balls)
    // we'd pop all the balls out to not be overlapping in a single iteration (Unless it's overlapping multiple balls). 
    // Let's slow this process a bit as it feels like it could be unstable if it shoots stuff straight into other balls
    static constexpr auto DiffScaleDivisor = 32 * 2;

    std::vector<Vector<Dim>> diffVect(state.size());
    for (auto & vec : diffVect)
    {
        vec.Zero();
    }

    for (size_t outerEpoch = 0; outerEpoch < OuterEpochs; outerEpoch++)
    {
        DEBUG_LOG("[OuterEpoch=%lu] Before Normalization:", outerEpoch);
        DEBUG_LOG_VECTS(state);
        Renormalize(state);
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
                    return;
                }
            }
            

            for (size_t i = 0; i < state.size(); i++)
            {
                for (size_t j = 0; j < Dim; j++)
                {
                    state[i].mValues[j] += ((diffVect[i].mValues[j] + DiffScaleDivisor) / DiffScaleDivisor);
                }
            }
        }
    }

    
    printf("No Soltuion found after all iterations:");
    PrintVects(state);

    // static constexpr auto 

    // In each outer loop we:
    // Rescale the vectors to have a magnitude around ScaledOne
    // construct a cache of points that are close enough to eachother that we check for colision
    // 
    // 
}


    

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

int main(int, char**){
    auto init = Initialize2D();
    RunRoutine<2>(init);
    return 0;
}
