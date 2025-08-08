
#pragma once




template <size_t Dim>
void CalcDotDiffs(std::vector<Vector<Dim>> const & points, NeighboursLookup const & neighbours, std::vector<Vector<Dim>> & rets)
{
    static constexpr double DELTA = 1e-6;
    static constexpr double QUAD_DELTA = 1;

    std::vector<PointType> mags(points.size());

    for (size_t i = 0; i < points.size(); i++)
    {
        mags[i] = std::sqrt(Dot(points[i], points[i]));

        // Apply force to keep kissing dist - quadratic unlike the linear forces for pushing away
        auto force = -1 * (mags[i] - ScaledOne) * std::abs(mags[i] - ScaledOne) * QUAD_DELTA;

        for (size_t j = 0; j < Dim; j++)
        {
            rets[i].mValues[j] = force * points[i].mValues[j];
        }
    }


    for (PointId pointId = 0; pointId < points.size(); pointId++)
    {
        auto const & point = points[pointId];
        auto & ret = rets[pointId];

        for (PointId neighbourId : neighbours[pointId])
        {
            auto const & neighbour = points[neighbourId];

            auto cos_theta = Dot(point, neighbour) / mags[pointId] / mags[neighbourId];

            if (cos_theta > 0.5) // points too close
            {
                SubMult(ret, neighbour, DELTA / mags[neighbourId]);
            }      
        }
    }
}

template <size_t Dim> 
bool RunGradientDescent(std::vector<Vector<Dim>> & initialState)
{
    FileOutput frameOutput("viewer/frames.json");
    auto & state = initialState;
    static constexpr size_t OuterEpochs = 100 * 1000;
    static constexpr size_t InnerIterationLoops = 1;


    std::vector<Vector<Dim>> diffVect(state.size());
    for (auto & vec : diffVect)
    {
        vec.Zero();
    }

    for (size_t outerEpoch = 0; outerEpoch < OuterEpochs; outerEpoch++)
    {
        auto neighbourLookup = ConstructPointNeighbours(state, ScaledBound(1.2));

        frameOutput.WriteRow(state);

        for (size_t innerEpoch = 0; innerEpoch < InnerIterationLoops; innerEpoch++)
        {
            
            CalcDotDiffs(state, neighbourLookup, diffVect);        
    	    
            for (size_t i = 0; i < state.size(); i++)
            {
                Acc(state[i], diffVect[i]);
            }
        }
    }

    PrintVects(state);

    return false;
}