#pragma once
#include "defs.hpp"
#include <cstddef>

namespace utils {

template<size_t N, size_t M>
inline bool check_bounds(const defs::Point& point) {
	return point.x >= 0 && point.x < N && point.y >= 0 && point.y < M;
}

template<size_t N, size_t M>
inline void clamp_point(defs::Point& point) {
	point.x  = std::clamp<float>(point.x, 0, N - 1);
	point.y  = std::clamp<float>(point.y, 0, M - 1);
}

template<size_t N , size_t M>
void gaussian_blur(defs::grid<N, M>& grid) {
	defs::grid<N, M> temp = grid;

	constexpr float kernel[3][3] = {
		{1/16.0f, 2/16.0f, 1/16.0f},
		{2/16.0f, 4/16.0f, 2/16.0f},
		{1/16.0f, 2/16.0f, 1/16.0f}
	};

	#pragma omp parallel for collapse(2)
	for(size_t x = 0; x < N; ++x) {
		for(size_t y = 0; y < M; ++y) {
			float acc = 0.0f;
			
			#pragma omp simd
			for(int kx = -1; kx <= 1; ++kx) {
				for(int ky = -1; ky <= 1; ++ky) {
					int nx = std::clamp<int>(x + kx, 0, N - 1);
					int ny = std::clamp<int>(y + ky, 0, M - 1);
					acc += temp[nx][ny] * kernel[kx + 1][ky + 1];
				}
			}

			grid[x][y] = acc;
		}
	}
}

} // namespace utils
