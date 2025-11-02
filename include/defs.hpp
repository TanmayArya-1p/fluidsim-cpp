#pragma once
#include <bits/stdc++.h>
#include <cstddef>

#define x first
#define y second

namespace defs {

template <size_t N, size_t M> using grid = std::array<std::array<float, M>, N>;

using Point = std::pair<size_t, size_t>;

template <class T> using vec2D = std::pair<T, T>;

} // namespace defs
