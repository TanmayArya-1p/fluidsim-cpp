#pragma once
#include <bits/stdc++.h>
#include<cstddef>
#include <SDL2/SDL.h>
#include "solver.hpp"

namespace display {

template <size_t sizeX, size_t sizeY>
class FluidDisplay {
    public:
        FluidDisplay(std::shared_ptr<solver::FluidGrid<sizeX, sizeY>> grid, size_t window_width, size_t window_height): grid(grid), window_width(window_width), window_height(window_height)
        {
            SDL_Init(SDL_INIT_VIDEO);
            window = SDL_CreateWindow("", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,window_width, window_height, SDL_WINDOW_SHOWN);
            renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
            
            reset_mimx();
            cell_w = static_cast<float>(window_width)/sizeX;
            cell_h = static_cast<float>(window_height)/sizeY;

        }

        ~FluidDisplay() {
            SDL_DestroyRenderer(renderer);
            SDL_DestroyWindow(window);
            SDL_Quit();
        }

        void reset_mimx() {
            mi_speed = std::numeric_limits<float>::max();
            mx_speed = std::numeric_limits<float>::lowest();
        }

        void render() {
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderClear(renderer);

            defs::grid<sizeX, sizeY> speeds;

            for (size_t x = 0; x < sizeX; ++x) {
                for (size_t y = 0; y < sizeY; ++y) {
                    if (grid->is_solid_cell({x, y})) continue;
                    
                    float vx = 0.5f * (grid->vel_x[x][y] + grid->vel_x[x+1][y]);
                    float vy = 0.5f * (grid->vel_y[x][y] + grid->vel_y[x][y+1]);
                    float speed = std::sqrt(vx*vx + vy*vy);
                    speeds[x][y] = speed;

                    if(std::isfinite(speed)) {
                        mi_speed = std::min(mi_speed, speed);
                        mx_speed = std::max(mx_speed, speed);
                    }
                }
            }

            if (mx_speed - mi_speed < 0.001f) {
                mx_speed = mi_speed + 0.1f;
            }

            for (size_t x = 0; x < sizeX; ++x) {
                for (size_t y = 0; y < sizeY; ++y) {
                    if (grid->is_solid_cell({x, y})) {
                        SDL_SetRenderDrawColor(renderer, 50, 50, 50, 255);

                    } else {
                        float norm_speed = (speeds[x][y] - mi_speed) / (mx_speed - mi_speed);
                        norm_speed = std::clamp(norm_speed, 0.0f, 1.0f);
                        norm_speed = std::pow(norm_speed, 0.5f);

                        Uint8 intensity = static_cast<Uint8>(norm_speed * 255);
                        SDL_SetRenderDrawColor(renderer, intensity, intensity, intensity, 255);
                    }

                    SDL_Rect cell = {
                        static_cast<int>(x * cell_w),
                        static_cast<int>(y * cell_h),
                        static_cast<int>(cell_w)+1,
                        static_cast<int>(cell_h)+1
                    };
                    SDL_RenderFillRect(renderer, &cell);
                }
            }

            SDL_RenderPresent(renderer);
        }


    private:
        std::shared_ptr<solver::FluidGrid<sizeX, sizeY>> grid;
        SDL_Window* window;
        SDL_Renderer* renderer;
        float cell_w;
        float cell_h;

        size_t window_width;
        size_t window_height;

        float mi_speed = 0.0f;
        float mx_speed = 1.0f;
};

} // namespace display