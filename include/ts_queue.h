//
// Created by Diaz, Diego on 31.3.2023.
//

#ifndef PARALLEL_PARSING_THREAD_SAFE_QUEUE_H
#define PARALLEL_PARSING_THREAD_SAFE_QUEUE_H

#include<mutex>
#include<queue>
#include <condition_variable>

template<typename T>
class ts_queue {
public:
    ts_queue() = default;
    ~ts_queue() {
        done();
    };

    void push(const T& item) {
        {
            std::unique_lock guard(m_queue_lock);
            m_queue.push(item);
        }
        m_condition.notify_one();
    }

    void push(T&& item) {
        {
            std::unique_lock guard(m_queue_lock);
            m_queue.push(std::move(item));
        }
        m_condition.notify_one();
    }

    bool pop(T& item) {
        std::unique_lock guard(m_queue_lock);
        m_condition.wait(guard, [&]() { return !m_queue.empty() || m_done; });
        if(m_done) return false;
        item = std::move(m_queue.front());
        m_queue.pop();
        return true;
    }

    std::size_t size() const {
        std::unique_lock guard(m_queue_lock);
        return m_queue.size();
    }

    bool empty() const {
        std::unique_lock guard(m_queue_lock);
        return m_queue.empty();
    }

    void done() {
        {
            std::unique_lock guard(m_queue_lock);
            m_done = true;
        }
        m_condition.notify_all();
    }

private:
    using queue_t = std::queue<T>;
    queue_t m_queue;
    mutable std::mutex m_queue_lock;
    std::condition_variable m_condition;
    bool m_done = false;
};
#endif //PARALLEL_PARSING_THREAD_SAFE_QUEUE_H
