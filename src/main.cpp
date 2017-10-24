#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "FSE.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
	}
	
	// FSE is initialized here
	FSE fse = FSE(50.0*1.6/3.6, 30); // 50 mph speed limit in m/s, 30 m lane change distance

  h.onMessage([&fse,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

						// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

						// Configuration constants
						const double accel = 0.0018; // speed increment over 0.02 s; not entirely sure what it works out to...
						const int max_previous_pts = 10; // max number of points from previous path to use
						const double max_speed = 0.43; // slightly less than 50 mph
						const double s_increment = 30; // m; difference between points to fit the spline
						const int n_s_pts_future = 3; // number of points in the future to fit the spline
						const int n_trajectory_pts = 50; // number of points to include in each trajectory
						const double follow_time = 3; // s

						// Initialize spline points vectors
						vector<double> spline_ptsx;
						vector<double> spline_ptsy;

						// First, figure out where the car will be and will be going at the end of any previous points that will be followed
						// Determine how many of the previous points to use (lesser of max_previous_pts and the size of the vector)
						int prev_size = previous_path_x.size();
						int n_prev_pts = min(max_previous_pts, prev_size);

						// Reference conditions
						double speed = car_speed*1.6/3.6*0.02;
						double ref_x = car_x;
						double ref_y = car_y;
						double ref_yaw = deg2rad(car_yaw);

						// If previous size is empty, use the car as a starting reference // taken from Udacity walkthrough
						if (prev_size == 0) {
							double prev_car_x = car_x - cos(car_yaw);
							double prev_car_y = car_y - sin(car_yaw);

							spline_ptsx.push_back(prev_car_x);
							spline_ptsx.push_back(car_x);
							spline_ptsy.push_back(prev_car_y);
							spline_ptsy.push_back(car_y);

							// If previous size is 1, use the next point and the current position (this should never happen)
						} else if (prev_size == 1) {
							ref_x = previous_path_x[n_prev_pts - 1];
							ref_y = previous_path_y[n_prev_pts - 1];

							double ref_x_prev = car_x;
							double ref_y_prev = car_y;

							ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

							spline_ptsx.push_back(ref_x_prev);
							spline_ptsx.push_back(ref_x);
							spline_ptsy.push_back(ref_y_prev);
							spline_ptsy.push_back(ref_y);

							speed = distance(ref_x, ref_y, ref_x_prev, ref_y_prev);				

						// Otherwise, use the last two of the previous points that will be included in the new trajectory // taken from Udacity walkthrough
						} else {
							ref_x = previous_path_x[n_prev_pts - 1];
							ref_y = previous_path_y[n_prev_pts - 1];

							double ref_x_prev = previous_path_x[n_prev_pts - 2];
							double ref_y_prev = previous_path_y[n_prev_pts - 2];

							ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

							spline_ptsx.push_back(ref_x_prev);
							spline_ptsx.push_back(ref_x);
							spline_ptsy.push_back(ref_y_prev);
							spline_ptsy.push_back(ref_y);

							speed = distance(ref_x, ref_y, ref_x_prev, ref_y_prev);

						}
						// cout << "prev_size=" << prev_size << "; speed=" << speed << endl;

						/**
						 * Criteria for a lange change:
						 * 1. Car in center of lane (i.e., a lane change is not ongoing already)
						 * 2. speed_target < max_speed (i.e., a car ahead is slowing us down now)
						 * 3. Cost function based on:
						 *   a. Is lane clear
						 *   b. Lane change more expensive than doing nothing
						 *   c. Relative speeds of next cars ahead in each lane
						 */
						double car_speed_mps = car_speed*1.6/3.6;
						double d_target = fse.update_d_target(car_s, car_d, car_speed_mps, sensor_fusion);
						
						// Build a set of points in frenet coordinates, and convert to global xy coordinates
						for (int i = 0; i < n_s_pts_future; i++) {
							vector<double> next_xy_global = getXY(car_s + (1+i)*s_increment, d_target, map_waypoints_s, map_waypoints_x, map_waypoints_y);
							spline_ptsx.push_back(next_xy_global[0]);
							spline_ptsy.push_back(next_xy_global[1]);
						}

						// Determine the speed_target
						double speed_target = max_speed;
						double left_bound = min(car_d, d_target) - 2;
						double right_bound = max(car_d, d_target) + 2;
						double follow_distance = follow_time*car_speed_mps;
						for (vector<double> car : sensor_fusion) {
							if (left_bound < car[6] && car[6] < right_bound && 0 < car[5] - car_s && car[5] - car_s < follow_distance) {
								// If a car is in the current or future lane, ahead of, and within follow_distance of the driven car, control speed.  At follow_distance, target will be max_speed, at distance = 0 (collision), target will be 0.
								speed_target = min(speed_target, (car[5] - car_s)/follow_distance*max_speed);
							}
						}
						
						// Shift from global xy to car xy coordinates (car xy system has x-axis in direction of the car, and with origin at the end of the previous points)
						// Copied from Udacity project walkthrough
						for (int i = 0; i < spline_ptsx.size(); i++) {
							double shift_x = spline_ptsx[i] - ref_x;
							double shift_y = spline_ptsy[i] - ref_y;

							spline_ptsx[i] = (shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw));
							spline_ptsy[i] = (shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw));
						}

						// Create and fit a spline
						tk::spline s;
						s.set_points(spline_ptsx, spline_ptsy);

						// Initialize the vector of target points, and initialize with previous points up to the limit
						vector<double> next_x_vals;
						vector<double> next_y_vals;

						for (int i = 0; i < n_prev_pts; i++) {
							next_x_vals.push_back(previous_path_x[i]);
							next_y_vals.push_back(previous_path_y[i]);
							// These will be in the global xy coordinates, but they're not being used to generate new points anymore
						}

						// Compute the trajectory, accelerating along the spline
						double x_pos_car = 0; // Initialize, car coordinates
						for (int i = 0; i < n_trajectory_pts - n_prev_pts; i++) {
							// Accelerate to speed_target
							if (speed < speed_target) {
								speed = min(speed_target, speed + accel);
							} else {
								speed = max(speed_target, speed - accel);
							}
							x_pos_car += speed;
							double y_pos_car = s(x_pos_car);

							// Convert back to global coordinates
							double x_ref = x_pos_car;
							double y_ref = y_pos_car;

							double x_pos_global = (x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw));
							double y_pos_global = (x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw));

							x_pos_global += ref_x;
							y_pos_global += ref_y;

							// Add next point
							next_x_vals.push_back(x_pos_global);
							next_y_vals.push_back(y_pos_global); // Neglects horizontal turning speed - should be ok for non-sharp turns
						}

						// End TODO

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
