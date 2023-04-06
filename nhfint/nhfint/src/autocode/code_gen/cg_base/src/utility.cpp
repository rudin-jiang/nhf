#include "cg_base/utility.hpp"
#include <cassert>
#include <cstddef>
#include <string>

std::size_t shell_size(std::size_t sum) {
    assert(sum < 100);
    return (sum + 1) * (sum + 2) / 2;
}

std::size_t pair_size(std::size_t sumA, std::size_t sumB) {
    return shell_size(sumA) * shell_size(sumB);
}

std::size_t ijk_index(std::size_t i, std::size_t j, std::size_t k) {
    assert(i < 100);
    assert(j < 100);
    assert(k < 100);

    return (j + k) * (j + k + 1) / 2 + k;
}

std::string encode_number(std::size_t n) {
    assert(n <= 35);
    if (n < 10) return std::to_string(n);
    return std::string(1, 'A' + int(n) - 10);
}

std::string array_name(std::size_t na, std::size_t nb) {
    return "a" + encode_number(na) + "b" + encode_number(nb);
}

std::size_t number_length(std::size_t num) {
    if (num == 0) return 1;
    std::size_t cnt = 0;
    while (num) {
        num /= 10;
        ++cnt;
    }
    return cnt;
}

std::string format_number(std::size_t num, std::size_t len) {
    std::size_t numLen = number_length(num);
    assert(numLen <= len);

    std::size_t spaceLen = len - numLen;
    std::string space = std::string(spaceLen, ' ');
    return space + std::to_string(num);
}