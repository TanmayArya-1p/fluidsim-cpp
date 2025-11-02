#pragma once
#include <cstddef>

#define CLAMP_BOUNDS 1

constexpr float DELTA_T = 0.001f;
constexpr int GAUSS_SEIDEL_ITERATIONS = 30;
constexpr float GAUSS_SEIDEL_SOR = 1.8f;

constexpr size_t GRID_SIZE_X = 80;
constexpr size_t GRID_SIZE_Y = 80;
constexpr float CELL_SIZE = 0.1f;
constexpr size_t WINDOW_WIDTH = 800;
constexpr size_t WINDOW_HEIGHT = 800;
constexpr float MOUSE_FORCE = 100.0f;

constexpr int TARGET_FPS = 125;