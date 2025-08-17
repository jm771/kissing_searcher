
#include <atomic>
#include <mutex>
#include <thread>
#include <chrono>
#include <deque>

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
            {
                std::scoped_lock lock {mMutex};
                if (mQueue.size())
                {
                    auto ret = mQueue.front();
                    mQueue.pop_front();
                    return ret;
                }
            }

            if (!mNRunningProducers)
            {
                return std::nullopt;
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