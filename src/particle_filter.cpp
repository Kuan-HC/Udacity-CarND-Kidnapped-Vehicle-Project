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
    particles.emplace_back(Particle{i, gaussRandom(x, std[0]), gaussRandom(y, std[1]), gaussRandom(theta, std[2]), 1});
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
   * DONE: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

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
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks)
{
  /**
   * DONE: Update the weights of each particle using a mult-variate Gaussian 
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

  /* 
      Steps:
      1. filter out the landmarks can not be seen from a specific particle, create a new vector stores detectable landmarks 
      2. fram transformation: from vehicle coordinates -> map coordinates
         Xmap = Xparticle + cos(theta) * Xcar - sin(theta) * Ycar
         Ymap = Yparticle + sin(theta) * Xcar + cos(theta) * Ycar
      3. Association part landmark to nearest observation
      4. Calculate the Particle's weight
   */
  weights.clear();

  for (auto &particle_i : particles)
  {     
    // Step 1 Select Landmarks
    vector<LandmarkObs> mark_avl;
    for (const auto &mark_i : map_landmarks.landmark_list)
    {
      double tmp_dist = dist(particle_i.x, particle_i.y, mark_i.x_f, mark_i.y_f);
      if (tmp_dist <= sensor_range)
        mark_avl.emplace_back(LandmarkObs{mark_i.id_i, mark_i.x_f, mark_i.y_f});
    }
    // End of Step 1

    // Step 2 Transormations
    vector<LandmarkObs> transed_observ;
    double tmp_cos_theta = cos(particle_i.theta);
    double tmp_sin_theta = sin(particle_i.theta);

    for (const auto &observ_i : observations)
    {
      double transed_x = particle_i.x + tmp_cos_theta * observ_i.x - tmp_sin_theta * observ_i.y;
      double transed_y = particle_i.y + tmp_sin_theta * observ_i.x + tmp_cos_theta * observ_i.y;
      
      transed_observ.emplace_back(LandmarkObs{-1,transed_x,transed_y});
    }
    // End of Step 2

    // step 3 Association
    dataAssociation(mark_avl, transed_observ);
    // End of Step 3

    // Step 4 Weight Calculate
    particle_i.weight = 1.0; 
    double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
    for (const auto &t_observ_i : transed_observ)
    {
      double exponent = pow(t_observ_i.x - map_landmarks.landmark_list[t_observ_i.id - 1].x_f, 2) / (2 * pow(std_landmark[0], 2))
                      + pow(t_observ_i.y - map_landmarks.landmark_list[t_observ_i.id - 1].y_f, 2) / (2 * pow(std_landmark[1], 2));
      double weight = gauss_norm * exp(-exponent);
      particle_i.weight *= weight;
    }
    
    // End of Step 4

    // step 5 Update Weight
    weights.push_back(particle_i.weight);
    // End of Step 5
  }
}

void ParticleFilter::resample()
{
  /**
   * DONE: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  /* resample base on http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution*/
  std::vector<Particle> resampled_particles;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> rand(weights.begin(), weights.end());

  for (int i = 0; i < num_particles; ++i)
  {
    int index = rand(gen);
    resampled_particles.emplace_back(particles[index]);
  }
  particles = std::move(resampled_particles);
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