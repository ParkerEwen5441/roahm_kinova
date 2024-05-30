#include "kinova_dynamics/KinovaControllerBlock.hpp"
#include "kinova_dynamics/RobustControlBlock.hpp"
#include "kinova_dynamics/Messages.hpp"
#include "kinova_dynamics/Parser.hpp"
#include "kinova_dynamics/Model.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using AccTrajectory = ::Roahm::msgs::AccTrajectory;

// Function to load trajectory data from a CSV line.
void loadTrajectoryFromCsv(const std::string &line, AccTrajectory &trajectory) {
    std::stringstream ss(line);
    std::vector<double> row;
    std::string value;

    while (std::getline(ss, value, ',')) {
        try {
            row.push_back(std::stod(value));  // Convert the string to a double and add to the row.
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument: " << e.what() << std::endl;
            return;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range: " << e.what() << std::endl;
            return;
        }
    }

    // Assuming 22 data points: 1 for time, 7 for angles, 7 for velocities, and 7 for accelerations
    if (row.size() != 22) {
        std::cerr << "Unexpected number of values in line: " << row.size() << std::endl;
        return;
    }

    trajectory.start_time = row[0];
    Eigen::Map<Eigen::VectorXd> angles(row.data() + 1, 7);
    Eigen::Map<Eigen::VectorXd> velocities(row.data() + 8, 7);
    Eigen::Map<Eigen::VectorXd> accelerations(row.data() + 15, 7);

    trajectory.pos = angles;
    trajectory.vel = velocities;
    trajectory.acc = accelerations;
}

int main() {
    using TrajectoryDeque = std::deque<std::shared_ptr<const AccTrajectory>>;

    // Initialize trajectory info
    AccTrajectory trajectory(7);

    // Use deque to store shared pointers to const AccTrajectory
    TrajectoryDeque trajectories;

    // Open the file
    std::ifstream file("/home/pewen/git/roahm_kinova/traj.csv");
    if (!file.is_open()) {
        std::cerr << "Could not open file /home/pewen/git/roahm_kinova/traj.csv" << std::endl;
        return 1;
    }

    std::string line;
    // Read each line from CSV
    while (std::getline(file, line)) {
        loadTrajectoryFromCsv(line, trajectory);

        // Here, make a copy of the trajectory in a shared_ptr and add to the deque
        std::shared_ptr<AccTrajectory> trajectoryPtr = std::make_shared<AccTrajectory>(trajectory);
        trajectories.push_back(trajectoryPtr);
    }

    for(auto& traj : trajectories){
        std::cout << traj->acc[0] << std::endl;
    }

    // file will be closed here automatically due to RAII

    Roahm::Model::Parser parse();

    const double fric_var = -1;
    const double model_variance = 0.03;
    const std::string block_name("kinova_contoller");
    const std::string model_path("/home/pewen/git/roahm_kinova/kinova_urdfs/gen3.urdf");
    Roahm::Model::SpatialModel<double> model = Roahm::Model::Parser::parse(model_path);
    Roahm::Dynamics::KinovaControlBlock kinova_controller(block_name, model, model_variance, fric_var);
    kinova_controller.set_trajectories(trajectories);
    // kinova_controller.run();

    return 0;
}