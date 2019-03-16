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

#define NUMBER_OF_PARTICLES 500 // TODO: Set the number of particles

double ParticleFilter::gaussian_noise(double mean, double std)
{
    std::normal_distribution<double> norm_dist(mean, std);
    return norm_dist(gen);
}

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs, double mu_x,
                                   double mu_y) const
{
    // calculate normalization term
    auto gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

    auto exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
        + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

    return gauss_norm * exp(-exponent);
}

void ParticleFilter::init(double x, double y, double theta, double std[])
{
    /**
     * TODO: Set the number of particles. Initialize all particles to 
     *   first position (based on estimates of x, y, theta and their uncertainties
     *   from GPS) and all weights to 1. 
     * TODO: Add random Gaussian noise to each particle.
     * NOTE: Consult particle_filter.h for more information about this method 
     *   (and others in this file).
     */
    num_particles = NUMBER_OF_PARTICLES;

    for (auto i = 0; i < num_particles; ++i)
    {
        Particle p;

        p.id = i;
        p.x = gaussian_noise(x, std[0]);
        p.y = gaussian_noise(y, std[1]);
        p.theta = gaussian_noise(theta, std[2]);
        p.weight = 1.0;

        particles.push_back(p);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate)
{
    /**
     * TODO: Add measurements to each particle and add random Gaussian noise.
     * NOTE: When adding noise you may find std::normal_distribution 
     *   and std::default_random_engine useful.
     *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
     *  http://www.cplusplus.com/reference/random/default_random_engine/
     */
    auto theta_dt = delta_t * yaw_rate;
    auto velocity_dt = velocity * delta_t;

    for (auto particle : particles)
    {
        if (fabs(yaw_rate) < 0.00001)
        {
            particle.x = particle.x + (velocity_dt * cos(particle.theta));
            particle.y = particle.y + (velocity_dt * sin(particle.theta));
        }
        else
        {
            particle.x = particle.x + (velocity / particle.theta * (sin(particle.theta + theta_dt) - sin(particle.theta)));
            particle.y = particle.y + (velocity / particle.theta * (cos(particle.theta) - cos(particle.theta + theta_dt)));
            particle.theta = particle.theta + theta_dt;
        }

        particle.x += gaussian_noise(0, std_pos[0]);
        particle.y += gaussian_noise(0, std_pos[1]);
        particle.theta += +gaussian_noise(0, std_pos[2]);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> landmark_predictions,
                                     vector<LandmarkObs>& observations)
{
    /**
     * TODO: Find the predicted measurement that is closest to each 
     *   observed measurement and assign the observed measurement to this 
     *   particular landmark.
     * NOTE: this method will NOT be called by the grading code. But you will 
     *   probably find it useful to implement this method and use it as a helper 
     *   during the updateWeights phase.
     */

    // select the id from the nearest landmark for each observation
    for (auto observation : observations)
    {
        auto min_dist = std::numeric_limits<double>::max();

        for (auto landmark : landmark_predictions)
        {
            auto current_dist = dist(landmark.x, landmark.y, observation.x, observation.y);

            if (current_dist < min_dist)
            {
                min_dist = current_dist;
                observation.id = landmark.id;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs>& observations,
                                   const Map& map_landmarks)
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

    for (auto particle : particles)
    {
        // transform observations from vehicle coordinates to map coordinates
        auto cos_theta = cos(particle.theta);
        auto sin_theta = sin(particle.theta);

        vector<LandmarkObs> transformed_observations;
        for (auto observation : observations)
        {
            // homogenous transformation
            LandmarkObs transformed_observation{
                -1, // we do not know with which landmark to associate this observation yet
                particle.x + cos_theta * observation.x - sin_theta * observation.y,
                particle.y + sin_theta * observation.x + cos_theta * observation.y
            };

            transformed_observations.push_back(transformed_observation);
        }

        // select all landmarks that are within the sensor range, measured from the current particle
        vector<LandmarkObs> landmark_predictions;
        for (auto landmark : map_landmarks.landmark_list)
        {
            if (dist(particle.x, particle.y, landmark.x_f, landmark.y_f) <= sensor_range)
                landmark_predictions.push_back(LandmarkObs{landmark.id_i, landmark.x_f, landmark.y_f});
        }

        // select the id from the nearest landmark for each observation
        dataAssociation(landmark_predictions, transformed_observations);

        // calculate and update the particle's final weight
        particle.weight = 1.0; // set to 1 for multiplication in the end of the loop
        for (auto observation : transformed_observations)
        {
            const auto landmark = landmark_predictions[observation.id - 1];
            particle.weight *= multiv_prob(std_landmark[0], std_landmark[1], observation.x, observation.y, landmark.x,
                                           landmark.y);
        }
    }
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

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y)
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
