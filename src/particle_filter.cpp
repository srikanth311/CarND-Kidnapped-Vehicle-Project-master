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
#include <random>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  // TODO: Set the number of particles
  // Here I set it to 10.
  num_particles = 10;
  // Creating default random engine - random number engine class to generate random numbers.
  default_random_engine gen;

  // Creating a Gaussian distribution for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  particles.resize(num_particles);
  weights.resize(num_particles);

  // Initializing the particles with x, y and theta using normal distribution values and weights to 1.

  for(int i=0; i<num_particles; i++)
  {
      particles[i].id = i;
      particles[i].x = dist_x(gen);
      particles[i].y = dist_y(gen);
      particles[i].theta = dist_theta(gen);
      particles[i].weight = 1.0;
      weights[i] = particles[i].weight;
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  // Gaussian distributions for noise

  normal_distribution<double> dist_x(0.0, std_pos[0]);
  normal_distribution<double> dist_y(0.0, std_pos[1]);
  normal_distribution<double> dist_theta(0.0, std_pos[2]);
  default_random_engine gen;

  // Now update particle x, y and theta

  for(int i = 0; i< num_particles; i++)
  {
      double x = particles[i].x;
      double y = particles[i].y;
      double theta = particles[i].theta;

      if (fabs(yaw_rate) > 0.001)
      {
          // If the yaw_rate value greater than 0 - for non zero yaw rate
          particles[i].x = x + ((velocity / yaw_rate) * (sin(theta + (yaw_rate * delta_t)) - sin(theta)));
          particles[i].y = y + ((velocity / yaw_rate) * (cos(theta) - cos(theta + (yaw_rate * delta_t))));
          particles[i].theta = theta + (yaw_rate * delta_t);
      } else {
          // For zero yaw rate
          particles[i].x = x + (velocity * delta_t * cos(theta));
          particles[i].y = y + (velocity * delta_t * sin(theta));
          particles[i].theta = theta;
      }

      // Adding random Gaussian noise

      particles[i].x += dist_x(gen);
      particles[i].y += dist_y(gen);
      particles[i].theta += dist_theta(gen);

  }


}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
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
                                   const Map &map_landmarks) {
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

  for(int i = 0; i< num_particles; i++)
  {
      vector<LandmarkObs> range_landmarks;
      Map::single_landmark_s coor_landmark;
      LandmarkObs landmarkObs;

      for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {

          coor_landmark = map_landmarks.landmark_list[j];

          // Check if landmarks in range of particles
          if (dist(coor_landmark.x_f , coor_landmark.y_f , particles[i].x,particles[i].y) <= sensor_range) {
              // Create a vector of landmarks in range
              landmarkObs.id = coor_landmark.id_i;
              landmarkObs.x = coor_landmark.x_f;
              landmarkObs.y = coor_landmark.y_f;
              range_landmarks.push_back(landmarkObs);
          }
      }

      vector<int> associations;
      vector<double> sense_x;
      vector<double> sense_y;

      // initialize weights
      particles[i].weight = 1.0;
      weights[i] = 1.0;

      // For every observation, transform coordinates, find the best landmark , calculate probability & update particle weights
      for (int jdx = 0; jdx < observations.size(); jdx++){

          LandmarkObs current_observation, best_landmark;
          landmarkObs = observations[jdx];

          // transform the coordinates to Map coordinates
          current_observation.x = particles[i].x + ((landmarkObs.x * cos(particles[i].theta)) - (landmarkObs.y * sin(particles[i].theta)));
          current_observation.y = particles[i].y + ((landmarkObs.x * sin(particles[i].theta)) + (landmarkObs.y * cos(particles[i].theta)));
          current_observation.id = observations[jdx].id;


          double distance_min = numeric_limits<double>::max();
          bool bestfound = false;

          //For this current observation, check every landmark in sensor range find the best landmark for this observation
          for (int kdx = 0; kdx < range_landmarks.size(); kdx++) {

              // Calculate distance between observation and landmark
              double distance = dist(range_landmarks[kdx].x, range_landmarks[kdx].y, current_observation.x, current_observation.y);

              // Note the landmark with the least distance delta as the best landmark
              if (distance < distance_min)
              {
                  distance_min = distance;
                  best_landmark = range_landmarks[kdx];
                  bestfound = true;
              }
          }

          if (bestfound) {
              // Calculate probability and update particle weights
              double normalizer = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
              double xterm = pow((current_observation.x - best_landmark.x) / std_landmark[0], 2) / 2.0;
              double yterm = pow((current_observation.y - best_landmark.y) / std_landmark[1], 2) / 2.0;
              double probw = normalizer * exp(-(xterm + yterm));
              particles[i].weight *= probw;
              weights[i] = particles[i].weight;
          }

          // Adding the details to associations
          associations.push_back(best_landmark.id);
          sense_x.push_back(current_observation.x);
          sense_y.push_back(current_observation.y);
      }
      // Assign associations, sense_x and sense_y to the particle
      SetAssociations(particles[i], associations, sense_x, sense_y);
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

    default_random_engine gen;

    // Create a discrete distribution for weights.
    discrete_distribution<int> dist_w(weights.begin(), weights.end());
    vector<Particle> r_particles;

    for (int i = 0; i < num_particles; i++)  {
        int widx = dist_w(gen);
        r_particles.push_back(particles[widx]);
    }

    particles = r_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;

  return particle;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}