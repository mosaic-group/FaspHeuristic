//
// Created by Krzysztof Gonciarz on 2019-03-02.
//

#ifndef APR2_TIMER_H
#define APR2_TIMER_H


#include <vector>
#include <chrono>
#include <ostream>
#include <string>
#include "tools/easylogging++.h"

template <bool SHOW_TIMINGS=false, bool PRINT_JUST_TIME=true, bool PRINT_INIT_MSG=false>
class Timer {
    std::vector<std::chrono::system_clock::time_point> iStartTimes;
    std::vector<std::string> names;
    bool iUseTimer;
    std::string iTimerName;

public:
    Timer(bool aUseTimer, const std::string &aTimerName = "NO_NAME") : iUseTimer(aUseTimer), iTimerName(aTimerName) {}
    ~Timer() {
        while(iStartTimes.size() > 0) {
            stop_timer();
        }
    }

    void start_timer(const std::string &timing_name) {
        if (iUseTimer) {
            std::ostringstream os;
            for (std::size_t i = 0; i < iStartTimes.size(); ++i) os << "    ";
            if (SHOW_TIMINGS && !PRINT_JUST_TIME && PRINT_INIT_MSG) LOG(TRACE) << os.str() << iTimerName << " [" << timing_name << "]";
            names.push_back(timing_name);
            std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();
            iStartTimes.push_back(startTime);
        }
    }

    double stop_timer() {
        if (iUseTimer) {
            std::chrono::system_clock::time_point endTime = std::chrono::system_clock::now();
            if (iStartTimes.size() == 0) {
                LOG(ERROR) << "ERROR: you have stopeed timer too many times..." << std::endl;
                return -1;
            }
            std::chrono::system_clock::time_point startTime = iStartTimes.back();
            iStartTimes.pop_back();
            std::chrono::duration<double> elapsedSeconds = endTime - startTime;
            auto name = names.back();
            names.pop_back();
            std::ostringstream os;
            if (SHOW_TIMINGS && !PRINT_JUST_TIME) for (std::size_t i = 0; i < iStartTimes.size(); ++i) os << "    ";
            double time = elapsedSeconds.count();
            if (SHOW_TIMINGS) {
                if (PRINT_JUST_TIME) {
                    LOG(TRACE) << time << std::endl;
                }
                else {
                    LOG(TRACE) << os.str() << iTimerName << " [" << name << "] = " << time;
                }
            }
            return time;
        }
        return -1;
    }
};


#endif //APR2_TIMER_H
