#pragma once
#include <bits/stdc++.h>
#include <cstddef>
#include "defs.hpp"
#include "utils.hpp"
#include "const.h"


namespace {
	std::mt19937 rng(std::random_device{}());
}


template <size_t sizeX , size_t sizeY>
class FluidGrid {
	public:
	   	FluidGrid(float cell_size, float density=10.0f, float delta_t=0.01f) {
		   	this->cell_size = cell_size;
		   	this->density = density;
		   	this->delta_t = delta_t;
	    }



		inline bool is_solid_cell(const defs::Point& point) {
			return (point.x<0 || point.x>=sizeX || point.y<0 || point.y>=sizeY);
		}

		float divergence(const defs::Point& point) const {

			float vel_top = this->vel_y[point.x][point.y+1];
			float vel_bottom = this->vel_y[point.x][point.y];
			float vel_left = this->vel_x[point.x][point.y];
			float vel_right = this->vel_x[point.x+1][point.y];

			float div = (vel_right - vel_left + vel_top - vel_bottom)/this->cell_size;
			return div;
		}

		float get_pressure(const defs::Point& point) const {
			if constexpr(CLAMP_BOUNDS) {
				utils::clamp_point<sizeX,sizeY>(point);
			} else {
				if(!utils::check_bounds<sizeX, sizeY>(point)) {
					throw std::out_of_range("Point out of bounds");
				}
			}
			return this->pressure[point.x][point.y];
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

			return (p_sum - this->density * this->cell_size * vel_out_flow / this->delta_t) / ct;
		}


		void compute_grid_pressure(float sor_w = 1.5f, int iterations = 30) {
			for (int iter = 0; iter < iterations; ++iter) {
				for (size_t x = 0; x < sizeX; ++x) {
					for (size_t y = 0; y < sizeY; ++y) {
						const float new_p = this->compute_cell_pressure({x, y});
						this->pressure[x][y] += sor_w * (new_p - this->pressure[x][y]);
					}
				}
			}
		}

		float _interpolate_velocity_x(const defs::Point& point) const {
			utils::clamp_point<sizeX,sizeY>(point);

			size_t x0 = point.x;
			size_t y0 = point.y-0.5f;

			x0 = std::clamp<size_t>(x0, 0, sizeX);
			y0 = std::clamp<size_t>(y0, 0, sizeY-1);

			float fx = point.x - x0;
			float fy = point.y - 0.5 - y0;

			float u = this->vel_x[x0][y0]*(1-fx)*(1-fy)
			u += fx*(1-fy)*this->vel_x[std::min(x0+1, sizeX)][y0];
			u += fx*(1-fy)*this->vel_x[std::min(x0+1, sizeX)][y0];
			u += (1-fx)*fy*this->vel_x[x0][std::min(y0+1, sizeY-1)];
			u += fx*fy*this->vel_x[std::min(x0+1, sizeX)][std::min(y0+1, sizeY-1)];

			return u;
		}

		float _interpolate_velocity_y(const defs::Point& point) const {
			utils::clamp_point<sizeX,sizeY>(point);

			size_t x0 = point.x-0.5f;
			size_t y0 = point.y;

			x0 = std::clamp<size_t>(x0, 0, sizeX-1);
			y0 = std::clamp<size_t>(y0, 0, sizeY);

			float fx = point.x - 0.5 - x0;
			float fy = point.y - y0;

			float v = this->vel_y[x0][y0]*(1-fx)*(1-fy)
			v += fx*(1-fy)*this->vel_y[std::min(x0+1, sizeX-1)][y0];
			v += (1-fx)*fy*this->vel_y[x0][std::min(y0+1, sizeY)];
			v += fx*fy*this->vel_y[std::min(x0+1, sizeX-1)][std::min(y0+1, sizeY)];

			return v;
		}


		defs::vec2D<float> interpolate_velocity(const defs::Point& point) const {
			return {_interpolate_velocity_x(point), _interpolate_velocity_y(point)};
		}

		void advect() {
			defs::grid<sizeX+1, sizeY> new_vel_x;
			defs::grid<sizeX, sizeY+1> new_vel_y;

			for (size_t x = 0; x < sizeX+1; ++x) {
				for (size_t y = 0; y < sizeY; ++y) {
					defs::Point curr_pos = {x, y + 0.5f};
					defs::vec2D<float> vel = interpolate_velocity(curr_pos);
					defs::Point prev_pos = {curr_pos.x - vel.x * this->delta_t/this->cell_size, curr_pos.y - vel.y * this->delta_t/this->cell_size};

					new_vel_x[x][y] = _interpolate_velocity_x(prev_pos);
				}
			}

			for(size_t x = 0; x < sizeX ; ++x) {
				for (size_t y = 0; y < sizeY+1; ++y) {
					defs::Point curr_pos = {x + 0.5f, y};
					defs::vec2D<float> vel = interpolate_velocity(curr_pos);
					defs::Point prev_pos = {curr_pos.x - vel.x * this->delta_t/this->cell_size, curr_pos.y - vel.y * this->delta_t/this->cell_size};

					new_vel_y[x][y] = _interpolate_velocity_y(prev_pos);
				}
			}

			this->vel_x = std::move(new_vel_x);
			this->vel_y = std::move(new_vel_y);
		}

		void enforce_velocity_boundaries() {
			//TODO: PARALLELIZE THIS
			for(size_t y = 0; y < sizeY; ++y) {
				this->vel_x[0][y] = 0.0f;
				this->vel_x[sizeX][y] = 0.0f;
			}

			for(size_t x = 0; x < sizeX; ++x) {
				this->vel_y[x][0] = 0.0f;
				this->vel_y[x][sizeY] = 0.0f;
			}
		}

		void update_velocity() {
			K = this->delta_t / (this->density * this->cell_size);
			//TODO: parallelize this

			#pragma omp parallel for collapse(2)
			for(size_t x = 0 ; x < sizeX+1; ++x) {
				for(size_t y = 0; y < sizeY; ++y) {
					defs::Point pt_right = {x, y};
					defs::Point pt_left = {x-1, y};

					if(this->is_solid_cell(pt) || this->is_solid_cell(pt_left)) {
						this->vel_x[x][y] = 0.0f;
						continue;
					}

					this->vel_x[x][y] -= K*(this->get_pressure(pt_right) - this->get_pressure(pt_left));
				}
			}

			for(size_t x = 0 ; x < sizeX; ++x) {
				for(size_t y = 0; y < sizeY+1; ++y) {
					defs::Point pt_top = {x, y};
					defs::Point pt_bottom = {x, y-1};

					if(this->is_solid_cell(pt) || this->is_solid_cell(pt_bottom)) {
						this->vel_y[x][y] = 0.0f;
						continue;
					}

					this->vel_y[x][y] -= K*(this->get_pressure(pt_top) - this->get_pressure(pt_bottom));
				}
			}
		}


		void randomize() {
			std::uniform_real_distribution<float> dist(-1.0f, 1.0f);

			for(size_t x = 0; x < sizeX+1; ++x) {
				for(size_t y = 0; y < sizeY; ++y) {
					this->vel_x[x][y] = dist(rng)*this->cell_size;
				}
			}

			for(size_t x = 0; x < sizeX; ++x) {
				for(size_t y = 0; y < sizeY+1; ++y) {
					this->vel_y[x][y] = dist(rng)*this->cell_size;
				}
			}

		}

		//TODO: GAUSSIAN FILTERING C++




	private:
		float density;
		float cell_size;
		float delta_t;

		defs::grid<sizeX+1, sizeY> vel_x;
		defs::grid<sizeX, sizeY+1> vel_y;
		defs::grid<sizeX, sizeY> pressure;

		std::map<std::pair<size_t, size_t>, bool> solid_cell_map;



};
