/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iterator>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
#include <cassert>

using std::string;
using std::vector;

#define NUMBER_OF_PARTICLES 100
#define EPSILON 0.000001

double ParticleFilter::gaussian_noise(double mean, double std)
{
    std::normal_distribution<double> norm_dist(mean, std);
    return norm_dist(gen);
}

double ParticleFilter::multiv_prob(double sig_x, double sig_y,
                                   double x_obs, double y_obs,
                                   double mu_x, double mu_y) const
{
    auto gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

    auto exponent = pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2))
        + pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2));

    return gauss_norm * exp(-exponent);
}

void ParticleFilter::init(double x, double y, double theta, double std[])
{
    num_particles = NUMBER_OF_PARTICLES;

    for (auto i = 0; i < num_particles; ++i)
    {
        particles.push_back(Particle{
            i, gaussian_noise(x, std[0]), gaussian_noise(y, std[1]), gaussian_noise(theta, std[2]), 1.0
        });
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate)
{
    auto theta_dt = yaw_rate * delta_t;
    auto velocity_dt = velocity * delta_t;

    for (auto& particle : particles)
    {
        // use simplified formula if yaw rate is ~0.0
        if (fabs(yaw_rate) < EPSILON)
        {
            particle.x += velocity_dt * cos(particle.theta);
            particle.y += velocity_dt * sin(particle.theta);
        }
        else
        {
            particle.x += velocity / yaw_rate * (sin(particle.theta + theta_dt) - sin(particle.theta));
            particle.y += velocity / yaw_rate * (cos(particle.theta) - cos(particle.theta + theta_dt));
            particle.theta += theta_dt;
        }

        // add noise
        particle.x += gaussian_noise(0, std_pos[0]);
        particle.y += gaussian_noise(0, std_pos[1]);
        particle.theta += gaussian_noise(0, std_pos[2]);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> landmark_predictions,
                                     vector<LandmarkObs>& observations)
{
    // associate the id from the nearest landmark for each observation
    for (auto& observation : observations)
    {
        assert(observation.id == -1);

        auto min_dist = std::numeric_limits<double>::max();

        for (auto landmark : landmark_predictions)
        {
            const auto current_dist = dist(landmark.x, landmark.y, observation.x, observation.y);

            if (current_dist < min_dist)
            {
                min_dist = current_dist;
                observation.id = landmark.id; // associate observation to landmark
            }
        }

        assert(observation.id != -1);
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs>& observations,
                                   const Map& map_landmarks)
{
    for (auto& particle : particles)
    {
        // transform observations from vehicle coordinates to map coordinates
        auto cos_theta = cos(particle.theta);
        auto sin_theta = sin(particle.theta);

        vector<LandmarkObs> transformed_observations;
        for (auto observation : observations)
        {
            // homogenous transformation
            LandmarkObs transformed_observation{
                -1, // no landmark associated
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

        // associate observations to landmarks
        dataAssociation(landmark_predictions, transformed_observations);

        // calculate and update the particle's final weight
        particle.weight = 1.0;
        for (auto& observation : transformed_observations)
        {
            assert(observation.id != -1);

            // find associated landmark
            LandmarkObs landmark{-1};
            for (auto landmark_prediction : landmark_predictions)
            {
                if (landmark_prediction.id == observation.id)
                {
                    landmark = landmark_prediction;
                    break;
                }
            }

            assert(landmark.id != -1);

            auto weight = multiv_prob(std_landmark[0], std_landmark[1], observation.x, observation.y, landmark.x, landmark.y);
            if (weight > 0)
                particle.weight *= weight;
        }
    }
}

void ParticleFilter::resample()
{
    auto index = rand() % static_cast<int>(num_particles);
    auto beta = 0.0;
    const auto max_weight = std::max_element(particles.begin(), particles.end(),
                                             [](const Particle& a, const Particle& b)
                                             {
                                                 return a.weight < b.weight;
                                             })->weight;

    std::uniform_real_distribution<double> dist_max_weight(0.0, max_weight);

    // resampling wheel
    vector<Particle> resampled_particles;
    for (auto i = 0; i < num_particles; i++)
    {
        beta += dist_max_weight(gen) * 2.0;
        while (beta > particles[index].weight)
        {
            beta -= particles[index].weight;
            index = (index + 1) % num_particles;
        }

        resampled_particles.push_back(particles[index]);
    }

    particles = resampled_particles;
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
