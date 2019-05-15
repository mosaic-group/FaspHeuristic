//
// Created by gonciarz on 2019-03-04.
//

#ifndef TOOLS_H
#define TOOLS_H

#include <sstream>
#include <iostream>
#include <string>
#include <cstring>


namespace Tools {

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
}
#endif
