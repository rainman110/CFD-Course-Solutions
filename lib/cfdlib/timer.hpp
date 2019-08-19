#ifndef TIMER_HPP
#define TIMER_HPP

/**
 * @brief Platform independent class for high
 * precision time measurements.
 */
class Timer {
public:
    Timer();
    void start();
    void stop();
    float milliseconds();
    virtual ~Timer();
private:
    class Timer_impl;
    Timer_impl * timer;
};

#endif // TIMER_HPP
