#include <vector>
#include <math.h>
#include <iostream>
#include "FSE.h"

using namespace std;

//
// FSE class definition implementation.
//

/**
 * Constructor
 * max_speed the speed limit of the road, in m/s
 * lane_change_dist expected distance over which lane change will occur, in m
 */
FSE::FSE(double _speed_limit, double _lane_change_distance) {
  speed_limit = _speed_limit;
  lane_change_distance = _lane_change_distance;
  d_target = 6; // setting the initial target at 6 m, i.e., center of middle lane.  It could be initialized to the center of the current lane, but that is unnecessary since we are 100% confident that the car will start near here.
  is_lane_change_ongoing = false;
}

// Destructor
FSE::~FSE() {}

/**
 * s_now current s frenet coordinate of the driven car, m
 * d_now current d frenet coordinate of the driven car, m
 * speed_now current absolute speed of the driven car, m/s
 * sf sensor fusion data
 * returns the new desired d frenet coordinate
 */
double FSE::update_d_target(double s_now, double d_now, double speed_now, const vector<vector<double>> &sf) {

  // Never change lanes if a lane change is already ongoing.  A lane change ends if the car gets close to its d_target.
  if (is_lane_change_ongoing && fabs(d_now - d_target) > max_d_error_lane_change_complete) {
    // Lane change ongoing
    return d_target;
  } else if (is_lane_change_ongoing) {
    // Lane change just ended
    is_lane_change_ongoing = false;
  }

  // Initialize cost, and calculate it for staying in the same lane
  double lowest_cost = calc_cost_from_slow_cars(d_target - 0.5*lane_width, d_target + 0.5*lane_width, s_now, sf);
  // cout << "stay_cost=" << lowest_cost << endl;

  // Set of new possible lanes
  double possible_new_lanes[] = {d_target - lane_width, d_target + lane_width};

  // Evaluate feasibility and costs of changing lanes

  // Define the safe zone
  double min_safe_s = s_now - safe_zone;
  double max_safe_s = s_now + safe_zone;

  // Evaluate each of the possible new lanes
  for (const double d_lane : possible_new_lanes) {
    // Check if lane exists
    if (d_lane <= 0 || d_lane >= right_edge_of_road) {
      continue; // If the lane doesn't exist, no further consideration necessary
    }

    // Initialize left and right bounds
    double left_bound = d_lane - 0.5*lane_width;
    double right_bound = d_lane + 0.5*lane_width;

    // Initialize cost
    double cost = 0;

    // Check if lane change would be safe
    // The concerns are:
    // 1. Is there another vehicle adjacent to the car
    // 2. Is there a vehicle coming up from behind in the new lane at high relative speed
    // Slow moving vehicles ahead in the new lane will stop us from making the lane change via the cost function, so needn't be considered here
    for (const vector<double> car : sf) {
      // Check if the car is even in the relevant lane
      if (left_bound < car[6] && car[6] < right_bound) {

        // If the car is in the safe zone, lane change impossible
        if (min_safe_s < car[5] && car[5] < max_safe_s) {
          cost = 10;
          break;
        }

        // If the car is expected to enter the safe zone over the course of the lane change, lane change impossible
        double car_speed = calc_hypotenuse(car[3], car[4]); // m/s
        double delta_speed = car_speed - speed_now; // m/s
        double lane_change_time = lane_change_distance/car_speed; // s
        double predicted_end_lane_change_position = car[5] + delta_speed*lane_change_time; // m
        if (car[5] < min_safe_s && predicted_end_lane_change_position > min_safe_s) {
          cost = 10;
          break;
        }
      }
    }

    // If we're still in this iteration of the loop, a lane change would be safe
    cost += lane_change_penalty_cost + calc_cost_from_slow_cars(left_bound, right_bound, s_now, sf);
    // cout << "Cost of d_target=" << d_lane << " is " << cost << endl;
    if (cost < lowest_cost) {
      lowest_cost = cost;
      d_target = d_lane;
      is_lane_change_ongoing = true;
    }
  }

  return d_target;
}

// Uses Pythagore to calculate the hypotenuse of a triangle; used here to calculate speeds
double FSE::calc_hypotenuse(double a, double b) {
  return sqrt(a*a + b*b);
}

// Cost is zero unless there's a close car in the same lane, in which case it's calculated based on the speed of the closest car
double FSE::calc_cost_from_slow_cars(double left_bound, double right_bound, double s_now, const vector<vector<double>> &sf) {
  double closest_car_s = s_now + max_look_distance; // Initialize high
  double cost = 0;
  for (const vector<double> car : sf) {
    if (s_now < car[5] && car[5] < closest_car_s && left_bound < car[6] && car[6] < right_bound) {
      closest_car_s = car[5];
      double car_speed = calc_hypotenuse(car[3], car[4]); // m/s
      cost = max(0.0, 1.0 - car_speed/speed_limit);
    }
  }
  return cost;
}
