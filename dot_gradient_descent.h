#pragma once

#include "weight_boosting.h"

template <size_t Dim>
void ApplyDiff(Vector<Dim> const & point, Vector<Dim> const & neighbour, double cos_theta, double scale, Vector<Dim> & ret)
{
    // SubMult(ret, neighbour, scale);

    // Lets just assume mags are close enough to 1...
    auto neighbourCopy = neighbour;

    // Orthoganlize the vector before scaling
    // This does seem to help avoid degeneracy. Unclear if it helps much in non degenerate cases.
    SubMult(neighbourCopy, point, cos_theta);
    Normalize(neighbourCopy, ScaledOne);
    SubMult(ret, neighbourCopy, scale);
}

template <size_t Dim>
void CalcDotDiffs(std::vector<Vector<Dim>> const & points, NeighboursLookup const & neighbours, std::vector<Vector<Dim>> & rets)
{
    static constexpr double DELTA = 1e-5;
    static constexpr double QUAD_DELTA = 1;

    std::vector<PointType> mags(points.size());

    double maxForce = 0;

    for (size_t i = 0; i < points.size(); i++)
    {
        mags[i] = std::sqrt(Dot(points[i], points[i]));
        rets[i].Zero();
    }
    for (PointId pointId = 0; pointId < points.size(); pointId++)
    {
        auto const & point = points[pointId];


        for (PointId neighbourId : neighbours[pointId])
        {
            auto & neighbour = points[neighbourId];

            auto cos_theta = Dot(point, neighbour) / mags[pointId] / mags[neighbourId];

            // boost[pointId].RegisterCosTheta(cos_theta);
            // boost[neighbourId].RegisterCosTheta(cos_theta);

            ASSERT_MSG(cos_theta <= 1.0000000001, "Cos theta was {}", cos_theta);

            static constexpr PointType RAMP_IN = 5;

            // Maybe we ramp this up over time instead?

            // if (cos_theta > 0.5 - (DELTA * RAMP_IN)) // points too close
            {
                // Give it this tiny bit of ramp in to try to help stability
                auto THRESH = 0.5 - (DELTA * RAMP_IN);
                auto scale = std::min(DELTA, (cos_theta - THRESH) / RAMP_IN);
                // auto scale = DELTA;


                auto heaviside = (cos_theta > THRESH);
                // auto heaviside = (cos_theta > THRESH) - (cos_theta <= THRESH && cos_theta > 0.45);
                // auto heaviside = 1 / (1 + std::exp(-10000 *(cos_theta - THRESH)));
                

                scale *= heaviside;

                // scale *= (1 + std::min(boost[pointId].GetBoostValue(), boost[neighbourId].GetBoostValue()));

                // What if we start chasing these numbers downwards after a while? "Cooling" as it were
                // Would be nice to have a convergence checker

                // Does a good job of preventing degenercy up to like 1.5, 1.6
                // Seed 12359 converges to an optimum until about 1.2. 
                // auto sf = exp(50 * (cos_theta - 0.5));
                // maxForce = std::max(sf, maxForce);
                auto sf = 1 / std::max(0.01, (1-cos_theta));
                scale *= sf;

                ApplyDiff(point, neighbour, cos_theta, scale, rets[pointId]);
                ApplyDiff(neighbour, point, cos_theta, scale, rets[neighbourId]);      
            }      
        }

    for (size_t i = 0; i < points.size(); i++)
    {
        // Apply force to keep kissing dist - quadratic unlike the linear forces for pushing away

        auto magError = (mags[i] - ScaledOne);
        auto forceScale = std::min(magError * magError, ScaledOne);
        auto force = std::signbit(magError) ? forceScale * QUAD_DELTA : -forceScale * QUAD_DELTA;

        for (size_t j = 0; j < Dim; j++)
        {
            // rets[i].mValues[j] /= maxForce;
            rets[i].mValues[j] += force * points[i].mValues[j] / mags[i];
        }
    }

        // boost[pointId].EndLoop();
    }
}

template <size_t Dim>
double CalcScore(std::vector<Vector<Dim>> const & state, NeighboursLookup & neighbourLookup)
{
    double score = 0;
    for (PointId pointId = 0; pointId < state.size(); pointId++)
    {
        for (PointId neighbourId : neighbourLookup[pointId])
        {
            auto dotVal = Dot(state[pointId], state[neighbourId]) / ScaledOneSquared;
            if (dotVal > 0.500000001)
            {
                score += (dotVal - 0.5);
            }
        }
    }

    return score;
}

template <size_t Dim, typename OutputT> 
double RunGradientDescent(std::vector<Vector<Dim>> & initialState, OutputT & frameOutput, size_t OuterEpochs, size_t InnerIterationLoops)
{
    auto & state = initialState;
    frameOutput.WriteRow(state);

    std::vector<Vector<Dim>> diffVect(state.size());
    // std::vector<BoostState> boost(state.size());
    for (auto & vec : diffVect)
    {
        vec.Zero();
    }

    for (size_t outerEpoch = 0; outerEpoch < OuterEpochs; outerEpoch++)
    {
        // std::cout << outerEpoch << std::endl;
        auto neighbourLookup = ConstructPointNeighbours(state, ScaledBound(1.2));
        frameOutput.WriteRow(state);

        for (size_t innerEpoch = 0; innerEpoch < InnerIterationLoops; innerEpoch++)
        {
            CalcDotDiffs<Dim>(state, neighbourLookup, diffVect);        
    	    
            for (size_t i = 0; i < state.size(); i++)
            {
                Acc(state[i], diffVect[i]);
            }
        }
    }

    Normalize(state, ScaledOne);

    auto neighbourLookup = ConstructPointNeighbours(state, ScaledBound(1.2));
    return CalcScore(state, neighbourLookup);

}

template <size_t Dim, typename OutputT> 
double RunGradientDescent(std::vector<Vector<Dim>> & initialState, OutputT & frameOutput)
{
    static constexpr size_t OuterEpochs = 20 * 1000;
    static constexpr size_t InnerIterationLoops = 100;

    return RunGradientDescent(initialState, frameOutput, OuterEpochs, InnerIterationLoops);
}