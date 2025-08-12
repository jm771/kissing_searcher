
#include "force_approach.h"
#include "simulated_annealing.h"
#include "dot_gradient_descent.h"

#include <atomic>
#include <mutex>
#include <thread>
#include <chrono>
#include <deque>

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

template<typename StoredT>
class ThreadSafeQueue
{
    public:
    ThreadSafeQueue(size_t nProducers) : mNRunningProducers(nProducers)
    {
    }

    void Push(StoredT && item)
    {
        std::scoped_lock lock{mMutex};
        mQueue.push_back(std::move(item));
    }

    std::optional<StoredT> PopWait()
    {
        while (true)
        {
            if (!mNRunningProducers)
            {
                return std::nullopt;
            }

            {
                std::scoped_lock lock {mMutex};
                if (mQueue.size())
                {
                    auto ret = mQueue.front();
                    mQueue.pop_front();
                    return ret;
                }
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }

    }

    bool Finished()
    {
        return ! mNRunningProducers;
    }

    void MarkFinishedProducer()
    {
        mNRunningProducers--;
    }


private:
    std::atomic<size_t> mNRunningProducers;
    std::mutex mMutex{};
    std::deque<StoredT> mQueue;
};

struct WorkResult
{
    size_t mSeed;
    double mScore;
};

// static constexpr size_t DIMENSION = 2; static constexpr size_t targetBalls = 6;
// static constexpr size_t DIMENSION = 3; static constexpr size_t targetBalls = 12;
// static constexpr size_t DIMENSION = 4; static constexpr size_t targetBalls = 24;
static constexpr size_t DIMENSION = 5; static constexpr size_t targetBalls = 40;

void workerThread(std::atomic<size_t> & inputQueue, ThreadSafeQueue<WorkResult> & resultQueue, size_t finishNumber)
{
    while(true)
    {
        size_t seed = inputQueue++;
        if (seed > finishNumber)
        {
            resultQueue.MarkFinishedProducer();
            return;
        }

        std::mt19937 rand(seed);
        auto state = Initialize<DIMENSION>(targetBalls, ScaledOne, rand);
        ASSERT(state.size() == targetBalls);
        Normalize(state, ScaledOne);
        auto score = RunGradientDescent<DIMENSION>(state);


        resultQueue.Push(WorkResult{seed, score});
    }

}



int main(int, char**){
     
    static constexpr size_t STARTING_SEED = 12345;
    static constexpr size_t STOPPING_SEED = 1234567;

    std::atomic<size_t> nextSeed{STARTING_SEED};

    // auto state = Initialize4D(rand);
    std::mt19937 rand(STARTING_SEED);
    auto state = Initialize5D(rand);

    static constexpr size_t nThreads = 15;

    ThreadSafeQueue<WorkResult> results{nThreads};

    std::vector<std::thread> threads;

    for (size_t i = 0; i < nThreads; i++)
    {
        threads.emplace_back([&nextSeed, &results]{ return workerThread(nextSeed, results, STOPPING_SEED);});
    }

    while (true)
    {
        auto entry = results.PopWait();
        if (!entry.has_value())
        {
            break;
        }

        std::cout << entry->mSeed << "," << entry->mScore << std::endl;
    }

    for (auto & thread : threads)
    {
        thread.join();
    }


    

    // auto state = StratifyPointsCubic<DIMENSION>(targetBalls);
    // DEBUG_LOG_VECTS(state);


    // if (success) {
    //     Validate(state);
    // }
    return 0;
}
