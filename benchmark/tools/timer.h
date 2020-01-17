//
// Created by Krzysztof Gonciarz on 2019-03-02.
//

#ifndef TIMER_H
#define TIMER_H


#include <vector>
#include <chrono>
#include <ostream>
#include <string>
#include <sstream>


/**
 * Very simple but usefull timer class. Usage:
 *
 * Timer<true> t("Timer Name");
 * t.start_timer("test1");
 * t.start_timer("test2");
 * t.stop_timer();
 * t.stop_timer();
 *
 * To have outputs like:
 * Timer Name [test1]
 *     Timer Name [test2]
 *     Timer Name [test2] = 1e-06s
 * Timer Name [test1] = 1.1e-05s
 *
 * If stop_timer() is not called same number of times as start_timer then it will be called automatically when
 * Timer object is destructed.
 *
 * @tparam SHOW_TIMINGS - true (default) to print timings
 * @tparam PRINT_JUST_TIME - false (default) to print just timeing values (without timer names etc.)
 * @tparam USE_TIMER - true (default), can be set to false to temporarily turn off this timer
 */
template <bool SHOW_TIMINGS=false, bool PRINT_JUST_TIME=false, bool USE_TIMER = true>
class Timer {
    std::vector<std::chrono::system_clock::time_point> iStartTimes;
    std::vector<std::string> names;
    std::string iTimerName;

public:
    Timer(const std::string &aTimerName = "TIMER") : iTimerName(aTimerName) {}
    ~Timer() {
        while(iStartTimes.size() > 0) {
            stop_timer();
        }
    }

    void start_timer(const std::string &timing_name) {
        if (USE_TIMER) {
            if (SHOW_TIMINGS && !PRINT_JUST_TIME) {
                std::ostringstream os;
                for (std::size_t i = 0; i < iStartTimes.size(); ++i) os << "    ";
                std::cout << os.str() << iTimerName << " [" << timing_name << "]" << "\n";
                names.push_back(timing_name);
            }
            std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();
            iStartTimes.push_back(startTime);
        }
    }

    double stop_timer() {
        if (USE_TIMER) {
            std::chrono::system_clock::time_point endTime = std::chrono::system_clock::now();
            if (iStartTimes.size() == 0) {
                std::cerr << "ERROR: you have stopeed timer too many times..." << std::endl;
                return -1;
            }
            std::chrono::system_clock::time_point startTime = iStartTimes.back();
            iStartTimes.pop_back();
            std::chrono::duration<double> elapsedSeconds = endTime - startTime;
            double time = elapsedSeconds.count();
            if (SHOW_TIMINGS) {
                if (PRINT_JUST_TIME) {
                    std::cout << time << "s\n";
                }
                else {
                    auto name = names.back();
                    names.pop_back();
                    std::ostringstream os;
                    for (std::size_t i = 0; i < iStartTimes.size(); ++i) os << "    ";
                    std::cout << os.str() << iTimerName << " [" << name << "] = " << time << "s\n";
                }
            }
            return time;
        }
        return -1;
    }
};


#endif
