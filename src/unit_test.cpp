#include <iostream>
#include "particle_filter.h"

using std::cout;
using std::endl;

int main()
{
    // Create particle filter
    ParticleFilter pf;

    // GPS measurement uncertainty [x [m], y [m], theta [rad]]
    double sigma_pos[3] = {0.3, 0.3, 0.01};

    //initialize particle
    pf.init(10.5, 5.0, 0.25, sigma_pos, 10);

    return 0;
}