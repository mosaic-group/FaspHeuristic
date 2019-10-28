//
// Created by gonciarz on 2019-03-04.
//

#ifndef TOOLS_H
#define TOOLS_H

#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <typeinfo>
#include <cxxabi.h>
#include <memory>
#include <algorithm>
#include <cmath>
#include <random>
#include <iomanip>

namespace Tools {

    [[maybe_unused]]
    static std::string demangle(const char* name) {
        int status = -4; // some arbitrary value to eliminate the compiler warning
        std::unique_ptr<char, void(*)(void*)> res {
                abi::__cxa_demangle(name, NULL, NULL, &status),
                std::free
        };

        return (status==0) ? res.get() : name ;
    }

    template <class T>
    static constexpr std::string_view type_name() {
        using namespace std;
        #ifdef __clang__
        string_view p = __PRETTY_FUNCTION__;
        return string_view(p.data() + 41, p.size() - 41 - 1);
        #elif defined(__GNUC__)
        string_view p = __PRETTY_FUNCTION__;
        #  if __cplusplus < 201402
        return string_view(p.data() + 36, p.size() - 36 - 1);
        #  else
        return string_view(p.data() + 49, p.find(';', 49) - 49);
        #  endif
        #elif defined(_MSC_VER)
        string_view p = __FUNCSIG__;
        return string_view(p.data() + 84, p.size() - 84 - 7);
        #endif
    }

    template <typename A>
    static void printType() {
        std::cout << type_name<A>() << std::endl;
    }

    template<typename A>
    static void printConstructorsAvailability() {
        std::ostringstream typeInfo{};
        typeInfo << "================= " << type_name<A>() << " =================";

        std::cout << typeInfo.str() << std::endl;
        std::cout << "constructible:      " << std::is_constructible<A>::value << " (default_constructible: " << std::is_default_constructible<A>::value << ")" << std::endl;
        std::cout << "destructible:       " << std::is_destructible<A>::value << std::endl;
        std::cout << "copy_assignable:    " << std::is_copy_assignable<A>::value << std::endl;
        std::cout << "copy_constructible: " << std::is_copy_constructible<A>::value << std::endl;
        std::cout << "move_assignable:    " << std::is_move_assignable<A>::value << std::endl;
        std::cout << "move_constructible: " << std::is_move_constructible<A>::value << std::endl;
        std::cout << std::string(typeInfo.str().length(), '=') << std::endl;
    }

    [[maybe_unused]]
    static bool endsWith(const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
    }

    [[maybe_unused]]
    static bool startsWith(const std::string& str, const std::string& prefix) {
        return str.size() >= prefix.size() && 0 == str.compare(0, prefix.size(), prefix);
    }

    [[maybe_unused]]
    static void replace(std::string &str, const std::string &from, const std::string &to) {
        if (from.empty()) return;

        size_t start_pos = 0;
        while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
            str.replace(start_pos, from.length(), to);
            start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
        }
    }

    [[maybe_unused]]
    static auto convertToStrWithLeadingZeros(int aNumber, uint16_t aNumOfLeadingZeros = 4) {
        std::stringstream ss;
        ss << std::setw(aNumOfLeadingZeros) << std::setfill('0') << aNumber;
        std::string s = ss.str();
        return s;
    }

    /**
     * Generates equally distributed and unique values in range [aMin, aMax]
     * @tparam T - type of generated elements
     * @param aMin - min value (first element of output)
     * @param aMax - max value (last element of output)
     * @param aNoOfSteps - requested number of steps (might be lower in output if elements repeats)
     * @return generated vector with values
     */
    template <typename T>
    auto linspace(T aMin, T aMax, int aNoOfSteps) {
        std::vector<T> result(aNoOfSteps);

        if (aNoOfSteps < 1) {
            throw std::runtime_error("Number of steps must be >= 1");
        }

        if (aNoOfSteps > 1) {
            double step = (static_cast<double>(aMax) - aMin) / (aNoOfSteps - 1);
            for (int i = 0; i < aNoOfSteps - 1; ++i) {
                result[i] = static_cast<T>(i * step + aMin);

            }
        }
        result[aNoOfSteps - 1] = aMax;

        // Make values in container unique
        auto it = std::unique(result.begin(), result.end());
        result.resize(distance(result.begin(), it));

        return result;
    }


    /**
 * Generates equally distributed and unique values in range [aMin, aMax]
 * @tparam T - type of generated elements
 * @param aMin - min value (first element of output)
 * @param aMax - max value (last element of output)
 * @param aNoOfSteps - requested number of steps (might be lower in output if elements repeats)
 * @return generated vector with values
 */
    template <typename T>
    auto logspace(T start, T stop, int num, double logBase = 3.0) {

        const double base = std::pow(logBase, 2.0/num);
        double value = 1.0;
        std::vector<double> retval; retval.reserve(num);
        std::generate_n(std::back_inserter(retval), num, [&](){double ret = value; value *= base; return ret;});

        std::vector<T> result(num);
        double b = retval.front();
        double e = retval.back();
        for (uint32_t i = 0; i < retval.size(); ++i) {
            auto v = retval[i];
            result[i] = static_cast<T>( floor( (stop-start) * (v-b) / (e-b) + start + 0.5) );
        }

        // Make values in container unique
        auto it = std::unique(result.begin(), result.end());
        result.resize(distance(result.begin(), it));

        return result;
    }

    auto randInt(int min, int max) {
        static std::mt19937 mt(std::random_device{}());
        return std::uniform_int_distribution<>(min, max)(mt);
    }
}
#endif
