#include <bits/stdc++.h>
#include "solver.hpp"
#include "const.h"
#include "utils.hpp"

int main() {

    const float density = 10.0f;
    const float cell_size = 0.1f;

    solver::FluidGrid<80, 80> f(cell_size, density);

    f.randomize();

    while(true) {
        f.advect();
        f.compute_grid_pressure();
        f.update_velocity();
        f.smooth_velocity();



        //display
    }

    return 0;
}