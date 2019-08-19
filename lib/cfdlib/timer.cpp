#include <cfdlib/timer.hpp>

#if !defined(_MSC_VER)
#include <sys/time.h>
#else
#include <windows.h>
#endif

class Timer::Timer_impl
{
public:
    Timer_impl(){};

    void start()
    {
#if defined(_MSC_VER)
        QueryPerformanceCounter(&start_t);
#else
        gettimeofday(&start_t, 0);
#endif
    }

    void stop()
    {
#if defined(_MSC_VER)
        QueryPerformanceCounter(&stop_t);
#else
        gettimeofday(&stop_t, 0);
#endif
    }

    float milliseconds()
    {
#if defined(_MSC_VER)
        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        return 1000.0f * (float)(stop_t.QuadPart - start_t.QuadPart) / (float)freq.QuadPart;
#else
        return float(stop_t.tv_sec - start_t.tv_sec) * 1000 + float(stop_t.tv_usec - start_t.tv_usec) / 1000;

#endif
    }

    virtual ~Timer_impl(){};

private:
#if defined(_MSC_VER)
    LARGE_INTEGER start_t, stop_t;
#else
    timeval start_t, stop_t;
#endif
};

Timer::Timer()
    : timer(new Timer_impl)
{
}

float Timer::milliseconds()
{
    return timer->milliseconds();
}

void Timer::start()
{
    timer->start();
}

void Timer::stop()
{
    timer->stop();
}

Timer::~Timer()
{
    delete timer;
}
