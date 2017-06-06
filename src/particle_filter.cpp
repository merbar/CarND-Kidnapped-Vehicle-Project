#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"
#include "helper_functions.h"


void ParticleFilter::init(double x, double y, double theta, const std::vector<double> &std) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 20;
  
  // create normal distributions
  std::default_random_engine gen;
  std::normal_distribution<double> N_x_init(x, std[0]);
  std::normal_distribution<double> N_y_init(y, std[1]);
  std::normal_distribution<double> N_theta_init(theta, std[2]);
  
  particles.clear();
  weights.clear();
  particles.reserve(num_particles);
  weights.reserve(num_particles);
  
  for (int i=0; i<num_particles; i++) {
      Particle new_particle;
      new_particle.id     = i;
      new_particle.weight = 1.0;
      new_particle.x      = N_x_init(gen);
      new_particle.y      = N_y_init(gen);
      new_particle.theta  = N_theta_init(gen);
      particles.push_back(new_particle);
      weights.push_back(1.0);
  }
  is_initialized = true;
}


void ParticleFilter::prediction(double delta_t, const std::vector<double> &std_pos, double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  
  // create normal distributions for process noise
  std::default_random_engine gen;
  std::normal_distribution<double> N_x(0.0, std_pos[0]);
  std::normal_distribution<double> N_y(0.0, std_pos[1]);
  std::normal_distribution<double> N_theta(0.0, std_pos[2]);
  
  for (int i=0; i<num_particles; i++) {
    if (std::fabs(yaw_rate) > 0.00001) {
      const double p_theta = particles[i].theta + delta_t * yaw_rate;
      particles[i].x += (velocity / yaw_rate) * (std::sin(p_theta)- std::sin(particles[i].theta));
      particles[i].y += (velocity / yaw_rate) * (std::cos(particles[i].theta) - std::cos(p_theta));
      particles[i].theta = p_theta;    
    } else {
      particles[i].x += velocity * std::cos(particles[i].theta) * delta_t;
      particles[i].y += velocity * std::sin(particles[i].theta) * delta_t;
    }
    
    // Add noise
    particles[i].x     += N_x(gen);
    particles[i].y     += N_y(gen);
    particles[i].theta += N_theta(gen);
  }  
}


void ParticleFilter::updateWeights(double sensor_range, const std::vector<double> &std_landmark, 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
  
  double weights_sum = 0.0;
  
  for (int i=0; i<num_particles; i++) {      
    // convert observed landmarks into world space from POV of particle
    std::vector<LandmarkObs> observations_particle_worldspace;
    
    const double cos_particle_theta = std::cos(particles[i].theta);
    const double sin_particle_theta = std::sin(particles[i].theta);
    
    for (int i_obs=0; i_obs<observations.size(); i_obs++) {
      LandmarkObs landmark;
      landmark.x = particles[i].x + observations[i_obs].x*cos_particle_theta - observations[i_obs].y*sin_particle_theta;
      landmark.y = particles[i].y + observations[i_obs].x*sin_particle_theta + observations[i_obs].y*cos_particle_theta;
      observations_particle_worldspace.push_back(landmark);
    }
    
    // associate known map landmark to observed landmarks
    // nearest-neighbor association
    std::vector<LandmarkObs> landmarks_to_observed;
    double weight = 1.0;
    
    for (int j=0; j<observations_particle_worldspace.size(); j++) {
      double min_distance = sensor_range;
      const Map::single_landmark_s* closest_lm = nullptr;
      double obs_x = observations_particle_worldspace[j].x;
      double obs_y = observations_particle_worldspace[j].y;
      
      for (int k=0; k<map_landmarks.landmark_list.size(); k++) {
        const double lm_x = map_landmarks.landmark_list[k].x_f;
        const double lm_y = map_landmarks.landmark_list[k].y_f;
        const double distance = dist(obs_x, obs_y, lm_x,  lm_y);
        if (distance < min_distance) {
          min_distance = distance;
          closest_lm = &map_landmarks.landmark_list[k];
        }
      }

      if (closest_lm) {
        // Multi-Variate Gaussian
        const double mvg_part1 = 1/(2.0*M_PI*std_landmark[0]*std_landmark[1]);
        const double x_dist    = obs_x - closest_lm->x_f;
        const double y_dist    = obs_y - closest_lm->y_f;
        const double mvg_part2 = ((x_dist*x_dist) / (2*std_landmark[0]*std_landmark[0])) + ((y_dist*y_dist) / (2*std_landmark[1]*std_landmark[1]));
        const double mvg       = mvg_part1 * std::exp(-mvg_part2);
        weight *= mvg;
      }
    }   
    weights[i] = weight;
    particles[i].weight = weight;
    weights_sum += weight;
  }
    
  // normalize weights
  for (int i=0; i<num_particles; i++)
      particles[i].weight /= weights_sum;
}


void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight. 
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dis(0, 1);
    std::vector<Particle> resampled_particles;
    resampled_particles.reserve(particles.size());
    for (int i=0; i<particles.size(); i++) {
        double rand_num = dis(gen);
        double particle_weights = 0.0;
        int j = 0;
        while (particle_weights < rand_num) {
            // shouldn't happen, but make sure we don't go out of bounds
            if (j >= particles.size())
              j = 0;
            particle_weights += particles[j].weight;
            j++;
        }
        resampled_particles.push_back(particles[j-1]);
    }
    particles = resampled_particles;
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
