#include<bits/stdc++.h>
#include<chrono>
#include<thread>
#include "solver.hpp"
#include "display.hpp"


int main() {

    auto f = std::make_shared<solver::FluidGrid<GRID_SIZE_X, GRID_SIZE_Y>>(CELL_SIZE);
    auto display = std::make_unique<display::FluidDisplay<GRID_SIZE_X, GRID_SIZE_Y>>(f, WINDOW_WIDTH, WINDOW_HEIGHT);

    bool running = true;
    bool paused = false;
    int frames = 0;

    // for (int i = 30; i < 50; ++i) {
    //     for (int j = 30; j < 50; ++j) {
    //         f->vel_x[i][j] = 5.0f;
    //     }
    // }

    auto start_time = std::chrono::high_resolution_clock::now();
    auto force_dir_distr = std::uniform_real_distribution<float>(-1, 1);

    while (running) {
        auto f_start_time = std::chrono::high_resolution_clock::now();

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_QUIT:
                    running = false;
                    break;

                case SDL_KEYDOWN:
                    switch (event.key.keysym.sym) {
                        case SDLK_ESCAPE:
                            running = false;
                            break;

                        case SDLK_SPACE:
                            paused = !paused;
                            std::cout << (paused ? "PAUSED" : "RUNNING") << std::endl;
                            break;
                    }
                    break;

                case SDL_MOUSEBUTTONDOWN:
                case SDL_MOUSEMOTION:
                    if (event.motion.state & SDL_BUTTON(SDL_BUTTON_LEFT)) {

                        int mx, my;
                        SDL_GetMouseState(&mx, &my);

                        float grid_x = (mx/static_cast<float>(WINDOW_WIDTH))*GRID_SIZE_X;
                        float grid_y = (my/static_cast<float>(WINDOW_HEIGHT))*GRID_SIZE_Y;

                        size_t cx = static_cast<size_t>(grid_x);
                        size_t cy = static_cast<size_t>(grid_y);
                        
                        if (cx < GRID_SIZE_X && cy < GRID_SIZE_Y) {
                            float force_x = MOUSE_FORCE * (force_dir_distr(rng));
                            float force_y = MOUSE_FORCE * (force_dir_distr(rng));

                            f->add_force({cx, cy}, {force_x, force_y});

                            for (int dx = -1; dx <= 1; ++dx) {
                                for (int dy = -1; dy <= 1; ++dy) {
                                    size_t nx = cx + dx;
                                    size_t ny = cy + dy;

                                    if (nx < GRID_SIZE_X && ny < GRID_SIZE_Y) {
                                        f->add_force({nx, ny}, {force_x, force_y});
                                    }
                                }
                            }
                        }
                    }
                    break;
            }
        }

        if (!paused) f->step();
        display->render();

        frames++;

        auto f_end_time = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration_cast<std::chrono::seconds>(f_end_time - start_time).count();
        
        if(frames>=100) {
            float fps = frames / static_cast<float>(dur);
            std::cout << "fps: " << fps << std::endl;
            frames = 0;
            start_time = f_end_time;
        }

        std::this_thread::sleep_for(std::chrono::microseconds(20));
    }
}