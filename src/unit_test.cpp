#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

#include "particle_filter.h"
#include "map.h"

using std::cout;
using std::endl;

// GPS measurement uncertainty [x [m], y [m], theta [rad]]
static double sigma_pos[3] = {0.3, 0.3, 0.01};
// Landmark measurement uncertainty [x [m], y [m]]
static double sigma_landmark[2] = {0.3, 0.3};

static double sensor_range = 10; // Sensor range [m]

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

#ifdef Transformation_Association_Weight_Test
    double std_landmark[2] = {0.3, 0.3};

    Particle particle_i;
    particle_i.x = 4.0;
    particle_i.y = 5.0;
    particle_i.theta = -M_PI / 2;

    double tmp_cos_theta = cos(particle_i.theta);
    double tmp_sin_theta = sin(particle_i.theta);

    LandmarkObs obser_i;
    obser_i.x = 2.0;
    obser_i.y = 2.0;

    double transed_x = particle_i.x + tmp_cos_theta * obser_i.x - tmp_sin_theta * obser_i.y;
    double transed_y = particle_i.y + tmp_sin_theta * obser_i.x + tmp_cos_theta * obser_i.y;
    /* result shall be  transed_x = 6 transed_y=3 */

    // Association
    std::vector<LandmarkObs> observations = {{-1, 6, 3}, {-1, 2, 2}, {-1, 0, 5}};
    std::vector<LandmarkObs> predicted = {{0, 5, 3}, {1, 2, 1}, {2, 6, 1}, {3, 7, 4}, {4, 4, 7}};

    for (auto &observ_1 : observations)
    {
        double min_dist = std::numeric_limits<float>::max();
        for (const auto &predict_i : predicted)
        {
            double tmp_dist = dist(observ_1.x, observ_1.y, predict_i.x, predict_i.y);
            if (min_dist > tmp_dist)
            {
                min_dist = tmp_dist;
                observ_1.id = predict_i.id;
            }
        }
    }
    /* calculate weight */
    double total_weight = 1.0;
    double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
    for (const auto &t_observ_i : observations)
    {
        double exponent = pow(t_observ_i.x - predicted[t_observ_i.id].x, 2) / (2 * pow(std_landmark[0], 2))
                        + pow(t_observ_i.y - predicted[t_observ_i.id].y, 2) / (2 * pow(std_landmark[1], 2));
        double weight = gauss_norm * exp(-exponent);
        total_weight *= weight;
    }

#endif
    // Create particle filter
    ParticleFilter pf;

    //initialize particle
    pf.init(4.0, 5.0, -M_PI / 2, sigma_pos, 5);

    /* move yaw rate 0*/
    pf.prediction(0.1, sigma_pos, 0.2, 0.0);
    pf.prediction(0.1, sigma_pos, 0.2, 0.2);

    Map map;
    //Map::single_landmark_s mark_1{1,5,3};
    map.landmark_list.emplace_back(Map::single_landmark_s{0, 5.0, 3.0});
    map.landmark_list.emplace_back(Map::single_landmark_s{1, 2.0, 1.0});
    map.landmark_list.emplace_back(Map::single_landmark_s{2, 6.0, 1.0});
    map.landmark_list.emplace_back(Map::single_landmark_s{3, 7.0, 4.0});
    map.landmark_list.emplace_back(Map::single_landmark_s{4, 4.0, 7.0});

    std::vector<LandmarkObs> noisy_observations;

    noisy_observations.push_back(LandmarkObs{0, 2, 2});
    noisy_observations.push_back(LandmarkObs{1, 3, -2});
    noisy_observations.push_back(LandmarkObs{2, 0, -4});

    pf.updateWeights(sensor_range, sigma_landmark, noisy_observations, map);

    pf.resample();

    return 0;
}