/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  // std::cout << x << " " << y << " " << theta << " " << std[0] << std::endl;
  num_particles = 500;
  std::default_random_engine gen;
  std::normal_distribution<double> N_x_init(x, std[0]);
  std::normal_distribution<double> N_y_init(y, std[1]);
  std::normal_distribution<double> N_theta_init(theta, std[2]);
  for (int i=0; i<num_particles; i++) {
      Particle new_particle;
      new_particle.id = i;
      new_particle.weight = 1.0;
      new_particle.x = N_x_init(gen);
      new_particle.y = N_y_init(gen);
      new_particle.theta = N_theta_init(gen);
      particles.push_back(new_particle);
  }
  is_initialized = true;
  //std::cout << particles.size() << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  for (int i=0; i<num_particles; i++) {
    if (yaw_rate != 0.0) {
      double p_theta = particles[i].theta + delta_t * yaw_rate;
      particles[i].x += (velocity / yaw_rate) * (sin(p_theta)- sin(particles[i].theta));
      particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(p_theta));
      particles[i].theta = p_theta;    
    } else {
      particles[i].x += velocity * sin(particles[i].theta) * delta_t;
      particles[i].y += velocity * cos(particles[i].theta) * delta_t;
    }
    
    // Add noise
    std::default_random_engine gen;
    std::normal_distribution<double> N_x(particles[i].x, std_pos[0]);
    std::normal_distribution<double> N_y(particles[i].y, std_pos[1]);
    std::normal_distribution<double> N_theta(particles[i].theta, std_pos[2]);
    particles[i].x += N_x(gen);
    particles[i].y += N_y(gen);
    particles[i].theta += N_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation 
  //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
  //   for the fact that the map's y-axis actually points downwards.)
  //   http://planning.cs.uiuc.edu/node99.html
  for (int i=0; i<num_particles; i++) {
    // convert observed landmarks into world space from POV of particle
    std::vector<LandmarkObs> observations_particle_worldspace;
    for (int i_obs=0; i_obs<observations.size(); i_obs++) {
      LandmarkObs landmark;
      landmark.id = observations[i_obs].id;
      landmark.x = particles[i].x +
              (observations[i_obs].x*cos(particles[i].theta) - (observations[i_obs].y*sin(particles[i].theta)));
      landmark.y = particles[i].y +
            (observations[i_obs].x*sin(particles[i].theta) + (observations[i_obs].y*cos(particles[i].theta)));
      observations_particle_worldspace.push_back(landmark);
    }
    
    // associate landmark observation to map landmark in particle space
    // nearest-neighbor association
    std::vector<LandmarkObs> landmarks_to_obs;
    //dataAssociation(landmarks_obs_particle_worldspace, &landmarks_to_obs);
    for (int i=0; i<observations_particle_worldspace.size(); i++) {
      double max_distance = sensor_range+1;
      int max_i = 0;
      for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
        double map_x = map_landmarks.landmark_list[j].x_f;
        double map_y = map_landmarks.landmark_list[j].y_f;
        double obs_x = observations_particle_worldspace[i].x;
        double obs_y = observations_particle_worldspace[i].y;
        double dist = sqrt((map_x - obs_x) * (map_x - obs_x) +
                           (map_y - obs_y) * (map_y - obs_y));
        if (dist < max_distance)
          max_distance = dist;
          max_i = j;
      }
      LandmarkObs landmark;
      landmark.x = map_landmarks.landmark_list[max_i].x_f;
      landmark.y = map_landmarks.landmark_list[max_i].y_f;
      landmarks_to_obs.push_back(landmark);
    }
    // compute weight (Multivariate Gaussian)
    double weight = 1.0;
    for (int i_obs=0; i_obs<observations.size(); i_obs++) {
      double mvg_part1 = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);
      double x_dist = landmarks_to_obs[i].x - particles[i].x;
      double y_dist = landmarks_to_obs[i].y - particles[i].y;
      double mvg_part2 = ((x_dist*x_dist) / (2*std_landmark[0]*std_landmark[0]))
                        + ((y_dist*y_dist) / (2*std_landmark[1]*std_landmark[1]));
      double mvg = mvg_part1 * exp(-mvg_part2);
      weight *= mvg;
      
              
    }
    particles[i].weight = weight;
  }
  // normalize weights
  
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}

