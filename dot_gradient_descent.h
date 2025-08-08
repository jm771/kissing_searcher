
#pragma once




template <size_t Dim>
void CalcDotDiffs(std::vector<Vector<Dim>> const & points, NeighboursLookup const & neighbours, std::vector<Vector<Dim>> & rets, std::vector<double> & boost, bool enableBoost)
{
    static constexpr double DELTA = 1e-5;
    static constexpr double QUAD_DELTA = 1;
    static constexpr double BOOST_THRESHOLD = 0.53;
    static constexpr double BOOST_RAMP = 0.001;
    static constexpr double BOOST_DECAY = 0.90;

    std::vector<PointType> mags(points.size());

    for (size_t i = 0; i < points.size(); i++)
    {
        mags[i] = std::sqrt(Dot(points[i], points[i]));

        // Apply force to keep kissing dist - quadratic unlike the linear forces for pushing away

        auto magError = (mags[i] - ScaledOne);
        auto forceScale = std::min(magError * magError, ScaledOne);
        auto force = std::signbit(magError) ? forceScale * QUAD_DELTA : -forceScale * QUAD_DELTA;

        for (size_t j = 0; j < Dim; j++)
        {
            rets[i].mValues[j] = force * points[i].mValues[j] / mags[i];
        }
    }


    for (PointId pointId = 0; pointId < points.size(); pointId++)
    {
        auto const & point = points[pointId];
        auto & ret = rets[pointId];

        bool hasCloseNeighbour = false;

        for (PointId neighbourId : neighbours[pointId])
        {
            auto const & neighbour = points[neighbourId];

            auto cos_theta = Dot(point, neighbour) / mags[pointId] / mags[neighbourId];
            if (enableBoost)
            {
                if (cos_theta > BOOST_THRESHOLD)
                {
                    hasCloseNeighbour = true;
                }
            }

            ASSERT_MSG(cos_theta <= 1.0000000001, "Cos theta was {}", cos_theta);

            static constexpr PointType RAMP_IN = 5;
            if (cos_theta > 0.5 - (DELTA * RAMP_IN)) // points too close
            {
                // Give it this tiny bit of ramp in to try to help stability
                auto scale = std::min(DELTA, (cos_theta - (0.5 - (DELTA * RAMP_IN))) / RAMP_IN);

                if (enableBoost)
                {
                    scale *= (1 + std::min(boost[pointId], boost[neighbourId]));
                }

                // Let's try this quadratic too (cos_theta - 0.5) *
                SubMult(ret, neighbour,  scale / mags[neighbourId]);
            }      
        }

        if (enableBoost)
        {
            if (hasCloseNeighbour) {
                boost[pointId] += BOOST_RAMP;
            } else {
                boost[pointId] *= BOOST_DECAY;
            }
        }
    }
}

template <size_t Dim> 
bool RunGradientDescent(std::vector<Vector<Dim>> & initialState, FileOutput & frameOutput, size_t OuterEpochs, size_t InnerIterationLoops, bool enableBoost)
{
    
    auto & state = initialState;
    frameOutput.WriteRow(state);



    std::vector<Vector<Dim>> diffVect(state.size());
    std::vector<double> boost(state.size());
    for (auto & vec : diffVect)
    {
        vec.Zero();
    }

    for (size_t outerEpoch = 0; outerEpoch < OuterEpochs; outerEpoch++)
    {
        auto neighbourLookup = ConstructPointNeighboursBidi(state, ScaledBound(1.2));
        frameOutput.WriteRow(state);

        for (size_t innerEpoch = 0; innerEpoch < InnerIterationLoops; innerEpoch++)
        {
            CalcDotDiffs(state, neighbourLookup, diffVect, boost, enableBoost);        
    	    
            for (size_t i = 0; i < state.size(); i++)
            {
                Acc(state[i], diffVect[i]);
            }
        }
    }

    PrintVects(state);

    return false;
}

template <size_t Dim> 
bool RunGradientDescent(std::vector<Vector<Dim>> & initialState)
{
    FileOutput frameOutput("viewer/frames.json");
    static constexpr size_t OuterEpochs = 50 * 1000;
    static constexpr size_t InnerIterationLoops = 100;

    RunGradientDescent(initialState, frameOutput, OuterEpochs, InnerIterationLoops, false);
    // return RunGradientDescent(initialState, frameOutput, OuterEpochs, InnerIterationLoops, true);
    return false;


}