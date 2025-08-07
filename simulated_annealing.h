



#pragma once


#include "types.h"
#include "debug_output.h"
#include "file_output.h"
#include "initial_states.h"
#include "vectors.h"
#include <stdint.h>
#include <random>

template <size_t Dim> 
double Energy(std::vector<Vector<Dim>> const & state, NeighboursLookup const & lookup)
{
    double ret = 0;
    std::cout << ret << std::endl;
    for (auto & vect : state)
    {
        auto e = labs(ScaledOne - std::sqrt(Dot(vect, vect)));
        ret += e;
        std::cout << ret << " " << e << std::endl;
    } 

    for (size_t id1 = 0; id1 < lookup.size(); id1++)
    {
        for (auto id2 : lookup[id1])
        {
            auto e = ScaledOne - Dist(state[id1], state[id2]);
            if (e > 0)
            {
                ret += e / 2;
            }
            std::cout << ret << " " << e << std::endl;
        }
    }

    return ret;
}

template <size_t Dim> 
double EnergyContrib(std::vector<Vector<Dim>> const & state, NeighboursLookup const & lookup, size_t elIdx, Vector<Dim> const & el)
{
    double ret = 0;
    {
        auto const & vect = el;
        auto e = abs(ScaledOne - std::sqrt(Dot(vect, vect)));
        ret += e;
    }

    for (auto id2 : lookup[elIdx])
    {
        auto e = ScaledOne - Dist(el, state[id2]);
        if (e > 0)
        {
            ret += e;
        }
    }

    return ret;
}

template <typename Rand> 
bool AcceptTransition(double oldEnergy, double newEnergy, double temperature, Rand & rand)
{
    if (newEnergy < oldEnergy)
    {
        return true;
    }
    else
    {
        std::uniform_real_distribution<double> realDistn(0, 1);
        return realDistn(rand) < std::exp(-(newEnergy - oldEnergy) / temperature);
    }
}




template <size_t Dim, typename Rand>
void RunInnerLoop(std::vector<Vector<Dim>> & state, NeighboursLookup const & neighbourLookup, double temperature, Rand & rand)
{
    // Let's ballpark a number of iterations whilst it's safe to not update neigbour lookup
    // Lets say each ball moves N times - that makes it move in each dimension a gaussian
    // with variance sqrt(N) * MoveDist.
    // The square magnitude of our vector is a chi squared distribution with parameter Dim
    // This has variance 2*Dim * sqrt(N) * MoveDist - Let's crudely approx as normal and say
    // We want something stepping out of the 0.2 move bounary to be a 5 sigma event
    // So we have 
    // 5 sigma = 0.1
    // -> sigma = 0.5 
    // -> variance = root(2) / 2
    // -> Sqrt(N) = root(2) / (2 * 2 * Dim * MoveDist)
    // -> N = 1 / 2 * Dim^2 * MoveDist ^ 2

    constexpr auto MOVE_DIST = ScaledOne >> 5;
    static constexpr auto RecipMoveDist = ScaledOne / MOVE_DIST;
    static constexpr auto SafeIterations = RecipMoveDist * RecipMoveDist / Dim / Dim / 2;


    std::uniform_int_distribution<size_t> intDistn(0, state.size() - 1);

    for (size_t innerEpoch = 0; innerEpoch < SafeIterations; innerEpoch++)
        {
            auto randomEl = intDistn(rand);
            auto randomMove = RandPoint<Dim>(MOVE_DIST, rand);
            auto oldEnergy = EnergyContrib(state, neighbourLookup, randomEl, state[randomEl]);
            auto newPoint = Add(state[randomEl], randomMove);
            auto newEnergy = EnergyContrib(state, neighbourLookup, randomEl, newPoint);

            if (AcceptTransition(oldEnergy, newEnergy, temperature, rand))
            {
                state[randomEl] = newPoint;
            }
        }
}


template <size_t Dim, typename Rand> 
bool RunAnnealing(std::vector<Vector<Dim>> & initialState, Rand & rand)
{
    FileOutput frameOutput("viewer/frames.json");
    auto & state = initialState;
    static constexpr size_t OuterEpochs = 100;

    auto neighbourLookup = ConstructPointNeighboursBidi(state, ScaledBound(1.2));    
    auto systemEnergy = Energy(state, neighbourLookup);
    std::cout << "System started with energy of " << systemEnergy << std::endl;
    // Allow a 1/4 increase in temp with probability 1/e
    auto const initialTemperature = systemEnergy / 4;

    // Let's just run the simulation a bit to quickly drop out of the really high initial temp.
    // Really high initial temp probably leads to some large movements so 
    for (size_t i = 0; i < 100000; i++)
    {
        neighbourLookup = ConstructPointNeighboursBidi(state, ScaledBound(1.2));
        RunInnerLoop(state, neighbourLookup, initialTemperature, rand);
    }

    systemEnergy = Energy(state, neighbourLookup);
    std::cout << "After initial rounds, system energy is " << systemEnergy << std::endl;

    // std::uniform_int_distribution<size_t> intDistn(0, state.size() - 1);

    // for (size_t outerEpoch = 0; outerEpoch < OuterEpochs; outerEpoch++)
    // {
    //     DEBUG_LOG_VECTS(state);
    //     frameOutput.WriteRow(state);
    //     neighbourLookup = ConstructPointNeighboursBidi(state, ScaledBound(1.2));

        

        
    // }

    
    printf("Finished found after %lu iterations:\n", OuterEpochs);
    PrintVects(state);

    return false;
}