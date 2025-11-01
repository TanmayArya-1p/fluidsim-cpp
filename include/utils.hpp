#pragma once
#include "defs.hpp"
#include <cstddef>

namespace utils {

template<size_t N, size_t M>
bool check_bounds(const defs::Point& point) {
	return point.x >= 0 && point.x < N && point.y >= 0 && point.y < M;
}

template<size_t N, size_t M>
void clamp_point(defs::Point& point) {
	point.x  = std::clamp<float>(point.x, 0, N - 1);
	point.y  = std::clamp<float>(point.y, 0, M - 1);
}

} //namespace utils
