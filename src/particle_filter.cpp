/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

static const float EPSINON = 0.00001;
static std::default_random_engine gen;

void ParticleFilter::init(const double &x, const double &y, const double &theta, const double std[], int value)
{
  /**
   * DONE: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * DONE: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  num_particles = value; // DONE: Set the number of particles
  
  /* create particles */
  for (int i = 0; i < num_particles; ++i)
  {
    Particle temp{i, gaussRandom(x, std[0]), gaussRandom(y, std[1]), gaussRandom(theta, std[2]), 1};
    particles.emplace_back(temp);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(const double &delta_t, const double std_pos[], const double &velocity, const double &yaw_rate)
{
  /**
   * DONE: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  for (auto &particle_i : particles)
  {
    /* yaw raw equals 0 */
    if ((yaw_rate >= -EPSINON) && (yaw_rate <= EPSINON))
    {
      particle_i.x += velocity * delta_t * cos(particle_i.theta);
      particle_i.y += velocity * delta_t * sin(particle_i.theta);
    }
    else
    {
      particle_i.x += velocity * (sin(particle_i.theta + yaw_rate * delta_t) - sin(particle_i.theta)) / yaw_rate;
      particle_i.y += velocity * (cos(particle_i.theta) - cos(particle_i.theta + yaw_rate * delta_t)) / yaw_rate;
      particle_i.theta += yaw_rate * delta_t;
    }
    /* adding Gaussian noise to each particleposition */
    particle_i.x += gaussRandom(0.0, std_pos[0]);
    particle_i.y += gaussRandom(0.0, std_pos[1]);
    particle_i.theta += gaussRandom(0.0, std_pos[2]);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs> &observations)
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks)
{
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
}

void ParticleFilter::resample()
{
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
}

void ParticleFilter::SetAssociations(Particle &particle,
                                     const vector<int> &associations,
                                     const vector<double> &sense_x,
                                     const vector<double> &sense_y)
{
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord)
{
  vector<double> v;

  if (coord == "X")
  {
    v = best.sense_x;
  }
  else
  {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}

double ParticleFilter::gaussRandom(const double &mean, const double mu)
{
  std::normal_distribution<double> gauss_dist(mean, mu);
  return gauss_dist(gen);
}