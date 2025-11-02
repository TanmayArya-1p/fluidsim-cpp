#pragma once
#include <bits/stdc++.h>
#include <cstddef>
#include "defs.hpp"
#include "utils.hpp"
#include "const.h"


namespace {
	std::mt19937 rng(std::random_device{}());
}


namespace solver {

template <size_t sizeX , size_t sizeY>
class FluidGrid {
	public:

	   	FluidGrid(float cell_size, float density=10.0f, float viscosity=0.01f) : viscosity(viscosity), cell_size(cell_size), density(density) {
			this->vel_x.fill({});
			this->vel_y.fill({});
			this->pressure.fill({});
	    }

		inline bool is_solid_cell(const defs::Point& point) const {
			return (point.x<0 || point.x>=sizeX || point.y<0 || point.y>=sizeY);
		}

		float divergence(const defs::Point& point) const {
			
			if(is_solid_cell(point)) {
				return 0.0f;
			}

			float vel_top = this->vel_y[point.x][point.y+1];
			float vel_bottom = this->vel_y[point.x][point.y];
			float vel_left = this->vel_x[point.x][point.y];
			float vel_right = this->vel_x[point.x+1][point.y];

			float div = (vel_right - vel_left + vel_top - vel_bottom)/this->cell_size;
			return div;
		}

		float get_pressure(const defs::Point& point) const {
			if constexpr(CLAMP_BOUNDS) {
				auto clamped_point = point;
				utils::clamp_point<sizeX,sizeY>(clamped_point);
				return this->pressure[clamped_point.x][clamped_point.y];
			} else {
				if(!utils::check_bounds<sizeX, sizeY>(point)) {
					throw std::out_of_range("Point out of bounds");
				}
				return this->pressure[point.x][point.y];
			}
		}

		float compute_cell_pressure(const defs::Point& point) {
			float p_sum = 0.0;
			float vel_out_flow = 0.0;
			uint ct = 0;

			if(!this->is_solid_cell({point.x, point.y+1})) {
				p_sum += this->get_pressure({point.x, point.y+1});
				vel_out_flow += this->vel_y[point.x][point.y+1];
				ct++;
			}

			if(!this->is_solid_cell({point.x, point.y-1})) {
				p_sum += this->get_pressure({point.x, point.y-1});
				vel_out_flow += -this->vel_y[point.x][point.y];
				ct++;
			}

			if(!this->is_solid_cell({point.x-1, point.y})) {
				p_sum += this->get_pressure({point.x-1, point.y});
				vel_out_flow += -this->vel_x[point.x][point.y];
				ct++;
			}

			if(!this->is_solid_cell({point.x+1, point.y})) {
				p_sum += this->get_pressure({point.x+1, point.y});
				vel_out_flow += this->vel_x[point.x+1][point.y];
				ct++;
			}

			return (p_sum - this->density * this->cell_size * vel_out_flow / DELTA_T) / ct;
		}


		void compute_grid_pressure() {

			for (int iter = 0; iter < GAUSS_SEIDEL_ITERATIONS; ++iter) {
				#pragma omp parallel for collapse(2)
				for (size_t x = 0; x < sizeX; ++x) {
					for (size_t y = 0; y < sizeY; ++y) {
						if (this->is_solid_cell({x, y})) {
							this->pressure[x][y] = 0.0f;
							continue;
						}

						const float new_p = this->compute_cell_pressure({x, y});
						this->pressure[x][y] += GAUSS_SEIDEL_SOR * (new_p - this->pressure[x][y]);
					}
				}
			}
		}

		float _interpolate_velocity_x(const defs::Point& pt) const {

			auto px = static_cast<float>(pt.x);
			auto py = static_cast<float>(pt.y);
			px = std::clamp(px, 0.0f, static_cast<float>(sizeX));
			py = std::clamp(py, 0.5f, static_cast<float>(sizeY) - 0.5f);
			
			int x0 = static_cast<int>(px);
			int y0 = static_cast<int>(py-0.5f);

			int x1 = std::min(x0+1, static_cast<int>(sizeX));
			int y1 = std::min(y0+1, static_cast<int>(sizeY-1));

			float fx = px - x0;
			float fy = (py - 0.5f) - y0;
			
			float vbl = vel_x[x0][y0];
			float vbr = vel_x[x1][y0];
			float vtl = vel_x[x0][y1];
			float vtr = vel_x[x1][y1];

			return (1-fx)*(1-fy)*vbl + fx*(1-fy)*vbr + (1-fx)*fy*vtl + fx*fy*vtr;
		}

		float _interpolate_velocity_y(const defs::Point& pt) const {

			auto px = static_cast<float>(pt.x);
			auto py = static_cast<float>(pt.y);
			px = std::clamp(px, 0.5f, static_cast<float>(sizeX)-0.5f);
			py = std::clamp(py, 0.0f, static_cast<float>(sizeY));
			
			int x0 = static_cast<int>(px-0.5f);
			int y0 = static_cast<int>(py);

			int x1 = std::min(x0 + 1, static_cast<int>(sizeX-1));
			int y1 = std::min(y0 + 1, static_cast<int>(sizeY));

			float fx = (px - 0.5f) - x0;
			float fy = py - y0;
			
			float vbl = vel_y[x0][y0];
			float vbr = vel_y[x1][y0];
			float vtl = vel_y[x0][y1];
			float vtr = vel_y[x1][y1];

			return (1-fx)*(1-fy)*vbl + fx*(1-fy)*vbr + (1-fx)*fy*vtl + fx*fy*vtr;
		}

		defs::vec2D<float> interpolate_velocity(const defs::Point& point) const {
			return {_interpolate_velocity_x(point), _interpolate_velocity_y(point)};
		}

		void advect() {
			defs::grid<sizeX+1, sizeY> new_vel_x;
			defs::grid<sizeX, sizeY+1> new_vel_y;


			#pragma omp parallel for collapse(2)
			for (size_t x = 0; x < sizeX+1; ++x) {
				for (size_t y = 0; y < sizeY; ++y) {
					defs::Point curr_pos = {x, y + 0.5f};
					defs::vec2D<float> vel = interpolate_velocity(curr_pos);
					defs::Point prev_pos = {curr_pos.x - vel.x * DELTA_T, curr_pos.y - vel.y * DELTA_T};

					new_vel_x[x][y] = _interpolate_velocity_x(prev_pos);
				}
			}


			#pragma omp parallel for collapse(2)
			for(size_t x = 0; x < sizeX ; ++x) {
				for (size_t y = 0; y < sizeY+1; ++y) {
					defs::Point curr_pos = {x + 0.5f, y};
					defs::vec2D<float> vel = interpolate_velocity(curr_pos);
					defs::Point prev_pos = {curr_pos.x - vel.x * DELTA_T, curr_pos.y - vel.y * DELTA_T};

					new_vel_y[x][y] = _interpolate_velocity_y(prev_pos);
				}
			}

			this->vel_x = std::move(new_vel_x);
			this->vel_y = std::move(new_vel_y);
		}

		void enforce_velocity_boundaries() {
			#pragma omp parallel for simd
			for(size_t y = 0; y < sizeY; ++y) {
				this->vel_x[0][y] = 0.0f;
				this->vel_x[sizeX][y] = 0.0f;
			}

			#pragma omp parallel for simd
			for(size_t x = 0; x < sizeX; ++x) {
				this->vel_y[x][0] = 0.0f;
				this->vel_y[x][sizeY] = 0.0f;
			}
		}

		void update_velocity() {
			float K = DELTA_T / (this->density * this->cell_size);
			
			#pragma omp parallel for collapse(2)
			for(size_t x = 1 ; x < sizeX; ++x) {
				for(size_t y = 0; y < sizeY; ++y) {
					defs::Point pt_right = {x, y};
					defs::Point pt_left = {x-1, y};

					this->vel_x[x][y] -= K*(this->get_pressure(pt_right) - this->get_pressure(pt_left));
				}
			}


			#pragma omp parallel for collapse(2)
			for(size_t x = 0 ; x < sizeX; ++x) {
				for(size_t y = 1; y < sizeY; ++y) {
					defs::Point pt_top = {x, y};
					defs::Point pt_bottom = {x, y-1};

					if(this->is_solid_cell(pt_top) || this->is_solid_cell(pt_bottom)) {
						this->vel_y[x][y] = 0.0f;
						continue;
					}
					this->vel_y[x][y] -= K*(this->get_pressure(pt_top) - this->get_pressure(pt_bottom));
				}
			}
		}

		void randomize() {
			std::uniform_real_distribution<float> dist(-1.0f, 1.0f);


			#pragma omp parallel for collapse(2)
			for(size_t x = 0; x < sizeX+1; ++x) {
				for(size_t y = 0; y < sizeY; ++y) {
					this->vel_x[x][y] = dist(rng)*this->cell_size;
				}
			}


			#pragma omp parallel for collapse(2)
			for(size_t x = 0; x < sizeX; ++x) {
				for(size_t y = 0; y < sizeY+1; ++y) {
					this->vel_y[x][y] = dist(rng)*this->cell_size;
				}
			}

		}

		void smoothen() {
			utils::gaussian_blur<sizeX+1, sizeY>(this->vel_x);
			utils::gaussian_blur<sizeX, sizeY+1>(this->vel_y);
			this->enforce_velocity_boundaries();
		}

		void add_force(defs::Point point, defs::vec2D<float> force) {
			if constexpr (CLAMP_BOUNDS) {
				utils::clamp_point<sizeX+1,sizeY>(point);
			} else {
				if(!utils::check_bounds<sizeX+1,sizeY>(point)) {
					throw std::out_of_range("Point out of bounds");
				}
			}
			this->vel_x[point.x][point.y] += force.x*DELTA_T/this->density;
			this->vel_y[point.x][point.y] += force.y*DELTA_T/this->density;
		}


		void step() {

			advect();
			enforce_velocity_boundaries();
			compute_grid_pressure();
			update_velocity();
			enforce_velocity_boundaries();

			#pragma omp parallel for collapse(2)
			for(size_t x = 0; x < sizeX+1; ++x) {
				for(size_t y = 0; y < sizeY; ++y) {
					vel_x[x][y] -= (vel_x[x][y] * this->viscosity);
				}
			}

			#pragma omp parallel for collapse(2)
			for(size_t x = 0; x < sizeX; ++x) {
				for(size_t y = 0; y < sizeY+1; ++y) {
					vel_y[x][y] -= (vel_y[x][y] * this->viscosity);
				}
			}
		}

		defs::grid<sizeX+1, sizeY> vel_x;
		defs::grid<sizeX, sizeY+1> vel_y;
	private:
		float density;
		float cell_size;
		float viscosity;

		defs::grid<sizeX, sizeY> pressure;

		std::map<std::pair<size_t, size_t>, bool> solid_cell_map;

};


} // namespace solver