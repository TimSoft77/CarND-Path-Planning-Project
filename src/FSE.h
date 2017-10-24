#ifndef FSE_H
#define FSE_H

#include <vector>

using namespace std;

class FSE {
private:
  // Variables
  double d_target;             // m; Essentially the state;
  bool is_lane_change_ongoing; // True once we start a lane change, false once complete.  Part of the state.
  double speed_limit;          // m/s - not mph or m in 0.02 s
  double lane_change_distance; // An estimate over how much distance the car will use to make a lane change.  Coupled with the speed of the car, this will provide an idea of how long lane changes take

  // Constants
  const double max_d_error_lane_change_complete = 0.1; // m
  const double lane_width = 4;                         // m
  const double max_look_distance = 100;                // m; never consider speeds of cars further away for lange changes
  const double lane_change_penalty_cost = 0.2;         // added to the cost of any lane change
  const double right_edge_of_road = 12;                // m; assuming left edge is d = 0
  const double safe_zone = 15;                         // m; distance on either side of s=0 (relative to car) that must remain clear during a lane change

  // Functions
  double calc_hypotenuse(double a, double b);
  double calc_cost_from_slow_cars(double left_bound, double right_bound, double s_now, const vector<vector<double>> &sf);

public:
  // Constructor
  FSE(double _speed_limit, double _lane_change_distance);

  // Destructor
  virtual ~FSE();

  // Functions
  double update_d_target(double s_now, double d_now, double speed_now, const vector<vector<double>> &sf);

};

#endif /* FSE_H */
