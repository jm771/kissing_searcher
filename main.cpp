
#include "force_approach.h"
#include "simulated_annealing.h"
#include "dot_gradient_descent.h"
#include "thread_safe_queue.h"

struct WorkResult
{
    size_t mSeed;
    double mStartScore;
    double mScore;
};

// static constexpr size_t DIMENSION = 2; static constexpr size_t targetBalls = 6;
// static constexpr size_t DIMENSION = 3; static constexpr size_t targetBalls = 12;
static constexpr size_t DIMENSION = 4; static constexpr size_t targetBalls = 24;
// static constexpr size_t DIMENSION = 5; static constexpr size_t targetBalls = 40;
// static constexpr size_t DIMENSION = 11; static constexpr size_t targetBalls = 593;

template <typename OutputT>
void workerThread(std::atomic<size_t> & inputQueue, ThreadSafeQueue<WorkResult> & resultQueue, OutputT & output, size_t finishNumber)
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

        // auto state = Initialize4D(rand);
        ASSERT(state.size() == targetBalls);
        Normalize(state, ScaledOne);

        auto neighbourLookup = ConstructPointNeighbours(state);
        auto startScore = CalcScore(state, neighbourLookup);

        auto score = RunGradientDescent<DIMENSION>(state, output);


        resultQueue.Push(WorkResult{seed, startScore, score});
    }
}



int main(int nargs, char** argv){
    ASSERT_MSG(nargs >= 2, "Missing arg - choose one of batch or analyse");
    std::string mode(argv[1]);
     
    size_t STARTING_SEED = 0;
    size_t STOPPING_SEED = 0;
    size_t nThreads = 0;

    FileOutput fileOutput{"viewer/frames.json"};
    NoOutput noOutput;

    if (mode == "batch")
    {
        STARTING_SEED = 12345;
        STOPPING_SEED = 1234567;
        // I'm gunna assume everything is hyperthreaded these days and that we don't want hyperthreads
        nThreads = std::thread::hardware_concurrency() /2 -1;
        ASSERT_MSG(nThreads > 0, "Could not determine thread count - pls hardcode");
        std::cerr << "Running on " << nThreads << " threads" << std::endl;
        nThreads = 7;
    }
    else if (mode == "analyse")
    {
        ASSERT_MSG(nargs >= 3, "use {} analyse <seed_number>", argv[0]);
        auto seed = std::stoll(argv[2]);
        STARTING_SEED = seed;
        STOPPING_SEED = seed;
        nThreads = 1;
    }
    else
    {
        ASSERT_MSG(false, "unkown mode");
    }

    std::atomic<size_t> nextSeed{STARTING_SEED};
    ThreadSafeQueue<WorkResult> results{nThreads};
    std::vector<std::thread> threads;

    for (size_t i = 0; i < nThreads; i++)
    {
        if (mode == "batch")
        {
            threads.emplace_back([&nextSeed, &results, STOPPING_SEED, &noOutput]{ return workerThread(nextSeed, results, noOutput, STOPPING_SEED);});
        }
        else
        {
            threads.emplace_back([&nextSeed, &results, STOPPING_SEED, &fileOutput]{ return workerThread(nextSeed, results, fileOutput, STOPPING_SEED);});
        }
    }

    while (true)
    {
        auto entry = results.PopWait();
        if (!entry.has_value())
        {
            break;
        }

        std::cout << "(" << entry->mSeed  << "," << entry->mStartScore << "," << entry->mScore << ")," << std::endl;
    }

    for (auto & thread : threads)
    {
        thread.join();
    }

    return 0;
}
