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

#include "particle_filter.h"

using namespace std;

double sigma_landmark [2] = {0.3, 0.3};


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    num_particles = 250;
    default_random_engine gen;

    normal_distribution<double> N_x(x,std[0]);
    normal_distribution<double> N_y(y,std[1]);
    normal_distribution<double> N_theta(theta,std[2]);

    for (int i = 0; i < num_particles; i++)
    {
        Particle particle;
        particle.id = i;
        particle.x = N_x(gen);
        particle.y = N_y(gen);
        particle.theta = N_theta(gen);
        particle.weight = 1;
        particles.push_back(particle);
        weights.push_back(1);
    }
    is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	default_random_engine gen;

	for (int i = 0; i < num_particles; i++)
    {
        double new_x;
        double new_y;
        double new_theta;

        if (yaw_rate == 0)
        {
            new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
            new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
            new_theta =  particles[i].theta;
        }
        else
        {
            new_x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
            new_y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
            new_theta =  particles[i].theta+yaw_rate*delta_t;
        }

        normal_distribution<double> N_x(new_x,std_pos[0]);
        normal_distribution<double> N_y(new_y,std_pos[1]);
        normal_distribution<double> N_theta(new_theta,std_pos[2]);

        particles[i].x = N_x(gen);
        particles[i].y = N_y(gen);
        particles[i].theta = N_theta(gen);


    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {

    for (int tt = 0; tt < num_particles; tt++)
    {
        vector<LandmarkObs> trans_observations;
        LandmarkObs obs;

        for ( int i = 0; i < observations.size();i++)
        {
            LandmarkObs trans_obs;
            obs = observations[i];
            trans_obs.x = particles[tt].x +(obs.x*cos(particles[tt].theta)-obs.y*sin(particles[tt].theta));
            trans_obs.y = particles[tt].y +(obs.x*sin(particles[tt].theta)+obs.y*cos(particles[tt].theta));
            trans_observations.push_back(trans_obs);
        }

        particles[tt].weight = 1.0;

        for(int i = 0 ; i < trans_observations.size(); i++)
        {
            double range = sensor_range;
            int found_landmark=0;

            for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
            {

                double landmark_x = map_landmarks.landmark_list[j].x_f;
                double landmark_y = map_landmarks.landmark_list[j].y_f;
                double calc_dist = sqrt(pow(trans_observations[i].x-landmark_x,2.0)+pow(trans_observations[i].y-landmark_y,2.0));

                if(calc_dist < range)
                {
                    range = calc_dist;
                    found_landmark=j;
                }
            }

            if(found_landmark!=0)
                {
                double meas_x = trans_observations[i].x;
                double meas_y = trans_observations[i].y;
                double mu_x = map_landmarks.landmark_list[found_landmark].x_f;
                double mu_y = map_landmarks.landmark_list[found_landmark].y_f;
                long double multipler = 1/(2*M_PI*std_landmark[0]*std_landmark[1])*exp(-(pow(meas_x-mu_x,2)/(2*pow(std_landmark[0],2))+pow(meas_y-mu_y,2)/(2*pow(std_landmark[1],2))));
                if (multipler > 0)
                {
                    particles[tt].weight*=multipler;
                }
                }
            weights[tt] = particles[tt].weight;
        }
    }
}





void ParticleFilter::resample() {

	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(),weights.end());
	vector<Particle> resample_particles;

	for (int i = 0; i<num_particles; i++)
    {
        resample_particles.push_back(particles[distribution(gen)]);
    }
    particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
