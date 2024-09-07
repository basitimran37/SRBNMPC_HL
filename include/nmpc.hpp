#ifndef NMPC_C
#define NMPC_C

#include "loco_structs.hpp"
#include <chrono>

using namespace casadi;
namespace fs = std::filesystem;

class NMPC {

private:
    
    Args args;     // Arguments for the solver
    Function f;
    double curr_time;
    Function f_dyn_;
    Function solver_;
    std::string basePath;

    // NMPC variables
    DM qk_; // Current state (qk) in the other code
    DM qkAll_; // All states for N horizons (X0) in the other code

    DM q_des_; // Desired state (xd) in the other code
    DM uk_; // Current control input - n_inputs x 1 
    DM ukAll_; // All control inputs for N horizons (u0) in the other code
    
    // Logger variables
    DM uout;
    DM xout;
    DM feet;
    DM traj;
    DM tout;
    DM NMPC_solve_time;

    // Solver inputs
    DM g_ub_, g_lb_;

public:
    Params params;

    NMPC()
    {
        // Initialize the parameters
        params.mpc_iter = 0;
        params.current_time = 0.0;
        params.dt = 0.01;
        params.Fs = 1/params.dt;
        params.sim_length = SIM_TIME;
        // params.L = 1000; // not sure if this is needed
        params.N = HORIZON;   // Horizon length
        params.height = 0.28;   // stand height
        params.n_states = 12;   // states
        params.n_inputs = 12;   // control inputs
        params.n_controls = 12; // control inputs
        params.gait = 8;        // Trot gait
        params.G = 9.81;        // gravity
        params.mass = 12.4530;  // mass
        params.J = DM::zeros(3,3);  // inertia matrix
        params.pf = DM::zeros(3, 4);    // nominal foot position matrix
        params.p_hip = DM::zeros(3, 4); // nominal hip position matrix
        params.Pr_refined_ = Eigen::MatrixXd::Zero(4, 4);
        params.Prd_refined_ = Eigen::MatrixXd::Zero(4, 4);
        params.agentIndex = 1;  // agent index

        // Initialize the pf matrix within the constructor body
        params.pf(0,0) = 0.183;  params.pf(1,0) = -0.1321; params.pf(2,0) = 0.01;
        params.pf(0,1) = 0.183;  params.pf(1,1) = 0.1321;  params.pf(2,1) = 0.01;
        params.pf(0,2) = -0.183; params.pf(1,2) = -0.1321; params.pf(2,2) = 0.01;
        params.pf(0,3) = -0.183; params.pf(1,3) = 0.1321;  params.pf(2,3) = 0.01;

        params.p_hip(0,0) = 0.183;  params.p_hip(1,0) = -0.1321; params.p_hip(2,0) = 0.01;
        params.p_hip(0,1) = 0.183;  params.p_hip(1,1) = 0.1321;  params.p_hip(2,1) = 0.01;
        params.p_hip(0,2) = -0.183; params.p_hip(1,2) = -0.1321; params.p_hip(2,2) = 0.01;
        params.p_hip(0,3) = -0.183; params.p_hip(1,3) = 0.1321;  params.p_hip(2,3) = 0.01;

        
        params.J(0,0) = 0.01683993;  params.J(0,1) = 8.3902e-5;   params.J(0,2) = 0.000597679;
        params.J(1,0) = 8.3902e-5;   params.J(1,1) = 0.056579028; params.J(1,2) = 2.5134e-5;
        params.J(2,0) = 0.000597679; params.J(2,1) = 2.5134e-5;   params.J(2,2) = 0.064713601;

        basePath = "tmp/";
    
    };
    ~NMPC();

    casadi::DM bezier(const DM& afra, double s);
    casadi::DM bezierd(const DM& alpha, double s);

    casadi::SX Rotz(const casadi::SX& t);
    casadi::SX skewSym(const casadi::SX& x);
    void planner();
    void force_inequality(const casadi::SX& U, casadi::SX& g_force);
    void writeMatrixToFile(const casadi::SX& matrix, const std::string& filename);
    void writeMatrixToFile(const casadi::DM& matrix, const std::string& filename);
    void writeMatrixToFile(const std::string& filename);
    void equalityconstraints(const casadi::SX& X, const casadi::SX& U, const casadi::SX& P, casadi::SX& obj, casadi::SX& g_eq, const Function& f);
    void set_constraint_bounds();
    void set_statebounds();
    void force_inequality_bounds();
    void shift();

    void code_gen();
    void run_NMPC();
    void init();
    void solve_nlp();
    void logData();
    inline std::vector<double> get_com_pos() {
        casadi::DM first_three = qk_(casadi::Slice(0, 3));  // Extract first three elements
        return std::vector<double>(first_three->begin(), first_three->end());
    }

    // Modify this function to return std::vector<double> for uk_
    inline std::vector<double> get_GRF() {
        return std::vector<double>(uk_->begin(), uk_->end());  // Convert the whole uk_ to std::vector
    }

    inline std::vector<double> get_contacts() {
        return std::vector<double>(params.contact->begin(), params.contact->end());  // Convert the whole uk_ to std::vector
    }

    inline std::vector<double> get_Pf() {
        casadi::DM Pf = casadi::DM::vertcat({
            params.pf(Slice(), 0), // First column
            params.pf(Slice(), 1), // Second column
            params.pf(Slice(), 2), // Third column
            params.pf(Slice(), 3)  // Fourth column
        });
        return std::vector<double>(Pf->begin(), Pf->end());  // Convert the whole uk_ to std::vector
    }
};
#endif // NMPC_H