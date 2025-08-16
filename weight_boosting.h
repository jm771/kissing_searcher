#pragma once

class BoostState
{
    static constexpr double BOOST_THRESHOLD = 0.53;
    static constexpr double EMERGENCY = 0.85;
    static constexpr double BOOST_RAMP = 0.005;
    static constexpr double BOOST_DECAY = 0.90;


    public:
    void RegisterCosTheta(double cos_theta)
    {
        if (cos_theta > BOOST_THRESHOLD)
        {
            hasCloseNeighbour = true;
        }
    }

    void EndLoop()
    {
        if (hasCloseNeighbour) {
            boostValue += BOOST_RAMP;
        } else {
            boostValue *= BOOST_DECAY;
        }

        hasCloseNeighbour = false;
    }

    double GetBoostValue()
    {
        return boostValue;
    }

    private:
    bool hasCloseNeighbour{};
    double boostValue;
};