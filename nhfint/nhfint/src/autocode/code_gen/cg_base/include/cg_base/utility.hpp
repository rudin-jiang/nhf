// some utility functions for generating code.

#include <string>
#include <vector>
#include <cstddef>

std::size_t shell_size(std::size_t sum);
std::size_t pair_size(std::size_t sumA, std::size_t sumB);
std::size_t ijk_index(std::size_t i, std::size_t j, std::size_t k);

// code gen utility

// 0 ~ 9, A -> 10, B -> 11, ...
std::string encode_number(std::size_t n);

// axbx
std::string array_name(std::size_t na, std::size_t nb);

// return the length of a number
std::size_t number_length(std::size_t num);

// print a number at fix len
std::string format_number(std::size_t num, std::size_t len);