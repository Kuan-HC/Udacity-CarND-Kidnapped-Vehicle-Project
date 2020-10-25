#include <iostream>
#include <cmath>
#include "particle_filter.h"

using std::cout;
using std::endl;

#define Motion_test false

int main()
{
#ifdef Motion_Test
    /* Motion Function Test*/
    Particle particle_i{0, 102, 65, 5 * M_PI / 8};
    double velocity = 110;
    double yaw_rate = M_PI / 8;
    double delta_t = 0.1;

    particle_i.x += velocity * (sin(particle_i.theta + yaw_rate * delta_t) - sin(particle_i.theta)) / yaw_rate;
    particle_i.y += velocity * (cos(particle_i.theta) - cos(particle_i.theta + yaw_rate * delta_t)) / yaw_rate;
    particle_i.theta += yaw_rate * delta_t;
/* result shall be :particle_i.x 97.59 particle_i.y 75.08 particle_i.theta 2.002700*/
#endif

    // Create particle filter
    ParticleFilter pf;

    // GPS measurement uncertainty [x [m], y [m], theta [rad]]
    double sigma_pos[3] = {0.3, 0.3, 0.01};

    //initialize particle
    pf.init(10.5, 5.0, 0.25, sigma_pos, 10);

    /* move yaw rate 0*/
    pf.prediction(1.0, sigma_pos, 1.0, 0.0);

    /* move yaw rate 1.0*/
    pf.prediction(1.0, sigma_pos, 1.0, 1.0);

    return 0;
}