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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <cfloat>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    num_particles = 100;
    
    default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
    
    for (int i =0; i < num_particles; i++)
    {
        Particle p = {};
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        p.weight = 1;
        
        particles.push_back(p);
    }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    default_random_engine gen;
    
    for (int i =0; i < num_particles; i++)
    {
        if (yaw_rate != 0)
        {
            particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
            particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
        }
        else
        {
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);
        }
        particles[i].theta += yaw_rate * delta_t;
        
        normal_distribution<double> dist_x(particles[i].x , std_pos[0]);
        normal_distribution<double> dist_y(particles[i].y , std_pos[1]);
        normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
        
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    
    for (unsigned int o = 0; o < observations.size(); o++)
    {
        double min_dist = DBL_MAX;
        int min_id = 0;
        for (unsigned int p = 0; p < predicted.size(); p++)
        {
            double distance =  dist(predicted[p].x, predicted[p].y,
                                    observations[o].x, observations[o].y);
            if (distance < min_dist)
            {
                min_dist = distance;
                min_id = predicted[p].id;
            }
        }
        observations[o].id = min_id;
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    for (int i =0; i < num_particles; i++)
    {    
        // Get predictions 
        std::vector<LandmarkObs> predicted;
        for (unsigned int l = 0; l < map_landmarks.landmark_list.size(); l++) 
        {
            Map::single_landmark_s lm = map_landmarks.landmark_list[l];
            // Assign predictions inside sensor range
            if (dist(particles[i].x, particles[i].y, lm.x_f, lm.y_f) <= sensor_range) 
            {
                LandmarkObs pred = {};
                pred.x = lm.x_f;
                pred.y = lm.y_f;
                pred.id = lm.id_i;
                predicted.push_back(pred);
            }
        }
        
        // Transform observations
        std::vector<LandmarkObs> observations_transform;
        for (unsigned int o = 0; o < observations.size(); o++)
        { 
            LandmarkObs mark = {};
            mark.id = observations[o].id;
            mark.x = particles[i].x + cos(particles[i].theta) * observations[o].x - sin(particles[i].theta) * observations[o].y;
            mark.y = particles[i].y + sin(particles[i].theta) * observations[o].x + cos(particles[i].theta) * observations[o].y;
            observations_transform.push_back(mark);
        }        
        
        // Associate
        dataAssociation(predicted,observations_transform);
        
        std::vector<int> associations;
        std::vector<double> sense_x;
        std::vector<double> sense_y;
        
        // Reset weight
        particles[i].weight = 1;
        
        // Update weights
        for (unsigned int o = 0; o < observations_transform.size(); o++) 
        {
            //  Find prediction associated with each observation
            int associate_id = observations_transform[o].id;
            associations.push_back(associate_id);
            
            double p_x, p_y = 0;
            for (unsigned int p = 0; p < predicted.size(); p++) 
            {
                if (predicted[p].id == associate_id) 
                {
                  p_x = predicted[p].x;
                  sense_x.push_back(p_x);
                  p_y = predicted[p].y;
                  sense_y.push_back(p_y);
                  break;
                }
            }
            
            // MVG
            double gaus_norm = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);            
            double exponent = pow(observations_transform[o].x - p_x,2)/(2 * pow(std_landmark[0],2)) + 
                              pow(observations_transform[o].y - p_y,2)/(2 * pow(std_landmark[1],2));
            particles[i].weight *= gaus_norm * exp(-exponent);

        }
                
        SetAssociations(particles[i],associations, sense_x, sense_y);
        break;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::vector<Particle> particles_new;

    
    // setup weights vector
    vector<double> weights;
    for (int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
    }
    
    // Setup discrete distribution
    default_random_engine gen;
    discrete_distribution<int> dist_index(weights.begin(), weights.end());
    
    // resample particles based on weights from discrete distribution 
    for (int i =0; i < num_particles; i++)
    {
        particles_new.push_back(particles[dist_index(gen)]);
    }
    
    particles = particles_new;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
