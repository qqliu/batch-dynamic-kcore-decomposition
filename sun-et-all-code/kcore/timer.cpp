#include "timer.h"

void Timer::reset() {
    start_ = std::chrono::steady_clock::now();
}

double Timer::elapsed() const {
    return std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now() - start_).count();
}

