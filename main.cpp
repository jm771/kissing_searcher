
#include "force_approach.h"
#include "simulated_annealing.h"

// static constexpr PointType GoogleScaledOne = 1e13;


// // TODO - commit to fixed point and do fixed point upper / lower ops with the __128 stuff

// static constexpr PointType RescaleGoogleValue(PointType googleValue)
// {
//     // 1e13 < 2^44 - so that leaves us 2^19 beore 2^63
//     // So we'll shift by 19 before div and 11 after to make the most of our 64 bit precision
//     return ((googleValue << 19) / GoogleScaledOne) << 11;
// }


// template <size_t Dim>
// PointType DistSquared (Vector<Dim> const & a, Vector<Dim> const & b)
// {
//     PointType ret = 0;
//     for (size_t i = 0; i < Dim; i++)
//     {
//         PointType diff = a.mValues[i] - b.mValues[i];
//         ret += diff * diff;
//     }
        
//     return ret;
// }



template <size_t Dim> 
bool Validate(std::vector<Vector<Dim>> & result)
{
    PointType maxSize = 0;
    DEBUG_LOG("Square Sizes:\n");
    for (auto const & vect : result)
    {
        auto size = Dot(vect, vect);
        DEBUG_LOG("%lu\n", size);
        maxSize = std::max(maxSize, size);
    }

    DEBUG_LOG("\n");

    PointType minDist = std::numeric_limits<PointType>::max();

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
    // static constexpr sifze_t DIMENSION = 2; static constexpr size_t targetBalls = 6;
    // static constexpr size_t DIMENSION = 3; static constexpr size_t targetBalls = 12;
    static constexpr size_t DIMENSION = 4; static constexpr size_t targetBalls = 24;
    // static constexpr size_t DIMENSION = 5; static constexpr size_t targetBalls = 40;


    std::mt19937 rand(12345);


    auto state = Initialize<DIMENSION>(targetBalls, ScaledOne, rand);
    // auto state = Initialize4D(rand);
    Normalize(state, ScaledOne);
    auto success = RunAnnealing<DIMENSION>(state, rand);

    if (success) {
        Validate(state);
    }
    return 0;
}
