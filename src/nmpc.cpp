// Author: Basit Muhammad Imran 
// Date: 2024-04-17

#include "nmpc.hpp"

SX blkdiag(const std::vector<SX>& matrices) {
    if (matrices.empty()) return SX();

    // Start with the first matrix
    SX result = matrices[0];
    for (size_t i = 1; i < matrices.size(); ++i) {
        // Pad the existing result and the next matrix with zeros to make them block diagonal
        SX topRight = SX::zeros(result.size1(), matrices[i].size2());
        SX bottomLeft = SX::zeros(matrices[i].size1(), result.size2());

        // Concatenate to form the new result
        result = vertcat(horzcat(result, topRight), 
                         horzcat(bottomLeft, matrices[i]));
    }

    return result;
}

casadi::DM NMPC::bezier(const DM& afra, double s) {
    int n = afra.size1(); // Number of rows
    int m = afra.size2(); // Number of columns
    DM value = DM::zeros(n);
    int M = m - 1;

    std::vector<int> k; // Binomial coefficients
    switch (M) {
        case 3: k = {1, 3, 3, 1}; break;
        case 4: k = {1, 4, 6, 4, 1}; break;
        case 5: k = {1, 5, 10, 10, 5, 1}; break;
        case 6: k = {1, 6, 15, 20, 15, 6, 1}; break;
        case 7: k = {1, 7, 21, 35, 35, 21, 7, 1}; break;
        case 8: k = {1, 8, 28, 56, 70, 56, 28, 8, 1}; break;
        case 9: k = {1, 9, 36, 84, 126, 126, 84, 36, 9, 1}; break;
        case 10: k = {1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1}; break;
        case 20: k = {1, 20, 190, 1140, 4845, 15504, 38760, 77520, 125970, 167960, 184756, 167960, 125970, 77520, 38760, 15504, 4845, 1140, 190, 20, 1}; break;
        default: std::cerr << "Degree " << M << " not supported." << std::endl; return value; // Optionally handle other degrees dynamically
    }

    std::vector<double> x(M + 1, 1.0), y(M + 1, 1.0);
    for (int i = 1; i <= M; ++i) {
        x[i] = s * x[i - 1];
        y[i] = (1 - s) * y[i - 1];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= M; ++j) {
            value(i) += afra(i, j) * k[j] * x[j] * y[M - j];
        }
    }

    return value;
}

DM NMPC::bezierd(const DM& alpha, double s) {
    int n = alpha.size1(); // Number of rows
    int m = alpha.size2(); // Number of columns
    DM value = DM::zeros(n);
    int M = m - 1; // Degree of the Bezier curve

    std::vector<int> k; // Binomial coefficients for the derivative
    switch (M) {
        case 3: k = {3, 6, 3}; break;
        case 4: k = {4, 12, 12, 4}; break;
        case 5: k = {5, 20, 30, 20, 5}; break;
        case 6: k = {6, 30, 60, 60, 30, 6}; break;
        case 7: k = {7, 42, 105, 140, 105, 42, 7}; break;
        case 8: k = {8, 56, 168, 280, 280, 168, 56, 8}; break;
        case 9: k = {9, 72, 252, 504, 630, 504, 252, 72, 9}; break;
        case 10: k = {10, 45, 120, 210, 252, 210, 120, 45, 10, 1}; break;
        case 20: k = {20, 380, 3420, 19380, 77520, 232560, 542640, 1007760, 1511640, 1847560, 1847560, 1511640, 1007760, 542640, 232560, 77520, 19380, 3420, 380, 20}; break;
        default: std::cerr << "Degree " << M << " not supported for derivative calculation." << std::endl; return value;
    }

    std::vector<double> x(M, 1.0), y(M, 1.0); // Adjust for derivative size
    for (int i = 1; i < M; ++i) {
        x[i] = s * x[i - 1];
        y[i] = (1 - s) * y[i - 1];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < M; ++j) {
            value(i) += (alpha(i, j + 1) - alpha(i, j)) * k[j] * x[j] * y[M - 1 - j];
        }
    }

    return value;
}

SX NMPC::Rotz(const SX& t) {
    // Ensure t is a vector of size 3
    if (t.size1() != 3 || t.size2() != 1) {
        throw std::runtime_error("Input vector t must be of size 3.");
    }

    casadi::SX Rz = casadi::SX::zeros(3, 3);

    // Note: In CasADi, indexing is zero-based, similar to C++ and Eigen but different from MATLAB.
    Rz(0, 0) = cos(t(1)) * cos(t(2));
    Rz(0, 1) = -sin(t(2));
    Rz(0, 2) = 0;
    Rz(1, 0) = cos(t(1)) * sin(t(2));
    Rz(1, 1) = cos(t(2));
    Rz(1, 2) = 0;
    Rz(2, 0) = -sin(t(1));
    Rz(2, 1) = 0;
    Rz(2, 2) = 1;

    return Rz;
}

casadi::SX NMPC::skewSym(const casadi::SX& x) {
    using namespace casadi;

    // Determine the size based on the input x
    int cols = x.size2();
    SX mat = SX::zeros(3, 3 * cols);

    // Loop through each column of x to create skew-symmetric matrices
    for (int n = 0; n < cols; ++n) {
        SX slice = x(Slice(), n);  // Get the nth column of x
        SX skew = SX::zeros(3, 3);
        // Fill the skew-symmetric matrix
        skew(0, 1) = -slice(2);
        skew(0, 2) = slice(1);
        skew(1, 0) = slice(2);
        skew(1, 2) = -slice(0);
        skew(2, 0) = -slice(1);
        skew(2, 1) = slice(0);

        // Place the skew matrix into the correct position in mat
        mat(Slice(0, 3), Slice(3*n, 3*(n+1))) = skew;
    }

    return mat;
}

void NMPC::planner() {
    const int n_states = params.n_states; // Example: Number of states
    // const int N = params.N;
    
    DM xd = DM::zeros(n_states * params.N);

    if (params.current_time <= 1.0) { // Stand up duration
        params.gait = 1;
    }
    else {
        params.gait = 8;
    }

    // writeMatrixToFile(q_des_, basePath + "q_des_.txt");

    DM traj = q_des_(casadi::Slice(0, n_states), 0); // Second parameter 0 is to indicate column index for a vector


    switch (params.gait) {
        case 1: {
            double t0 = 0, tf = 1.0; // Start and end times for mode 1
            double s = std::min((params.current_time - t0) / (tf - t0), 1.0); // Ensure s is within [0, 1]
            
            // Extract relevant parameters for easy access
            double xf = 0.0; // Initial x position
            double yf = 0.0; // Initial y position
            double zf = params.height; // Assuming the height is constant or predefined
            double z0 = 0.08;

            // Using DM::horzcat to create each row of the alpha matrix
            casadi::DM row1 = casadi::DM::horzcat({DM(0), DM(0), DM(0), DM(xf / 4), DM(3 * xf / 4), DM(xf), DM(xf), DM(xf)});
            casadi::DM row2 = casadi::DM::horzcat({DM(0), DM(0), DM(0), DM(yf / 4), DM(3 * yf / 4), DM(yf), DM(yf), DM(yf)});
            casadi::DM row3 = casadi::DM::horzcat({DM(z0), DM(z0), DM(z0), DM(zf / 4), DM(3 * zf / 4), DM(zf), DM(zf), DM(zf)});

            casadi::DM alpha = casadi::DM::vertcat({row1, row2, row3});

            // Update trajectory position
            casadi::DM position_update = bezier(alpha, s);

            for (int i = 0; i < 3; ++i) {
                traj(i) = position_update(i);
            }

            // Update trajectory velocity
            casadi::DM velocity_update = bezierd(alpha, s);

            for (int i = 0; i < 3; ++i) {
                traj(i + 3) = velocity_update(i);
            } 
            break;
        }
        case 8: {
            double domLen = 200.0; // Domain length in milliseconds
            double domLenSec = domLen / 1000.0; // Domain length in seconds
            DM Rz = Rotz(vertcat(qk_(6), qk_(7), qk_(8))); // Assuming q contains Euler angles in the last 3 elements
            // size_t currentStep = std::floor(t / params.dt) + 1;

            if (params.current_time >= (params.stepcnt * domLenSec + 1)) {
                // Update contact sequence
                if (static_cast<double>(params.contact(0)) == 1.0) {
                    params.contact = {0, 1, 1, 0};
                } else {
                    params.contact = {1, 0, 0, 1};
                }

                // Increase velocity in body frame
                traj(3) += 0.1; // Increase command velocity in body X
                traj(4) += 0.1; // Increase command velocity in body Y

                // Saturate the command velocities
                traj(3) = fmin(traj(3), params.desVel(0));
                traj(4) = fmin(traj(4), params.desVel(1));

                // Calculate step length in body frame and update step foot position
                DM stepLen = vertcat(traj(3) * domLenSec / 2, traj(4) * domLenSec / 2, 0.0);
                DM hipTmp = params.p_hip + stepLen;
                hipTmp = mtimes(Rz, hipTmp); // Rotate to the world frame
                params.pf = hipTmp + vertcat(qk_(Slice(0, 2)), 0.0);

                // Set velocity in world frame
                DM cmdWrld = mtimes(Rz, traj(Slice(3, 6))); // Convert to world frame

                // std::string basePath = "tmp/pathplanner/";
                // writeMatrixToFile(stepLen, basePath + "stepLen.txt");
                // writeMatrixToFile(params.p_hip, basePath + "params_p_hip.txt");
                // writeMatrixToFile(hipTmp, basePath + "hipTmp.txt");
                // writeMatrixToFile(q(Slice(0, 3)), basePath + "q_COM.txt");
                // writeMatrixToFile(cmdWrld, basePath + "cmdWrld.txt");
                // writeMatrixToFile(params.pf, basePath + "params_pf.txt");
                traj(Slice(3, 6)) = cmdWrld;
                // Increment step counter
                params.stepcnt++;
            }

            float heading = 0;

            // Update position and height based on velocity
            traj(0) = qk_(0) + params.dt * traj(3);
            traj(1) = qk_(1) + params.dt * traj(4);
            traj(2) = params.height; // Assuming height is constant or predefined
            traj(5) = 0; // Assuming zero vertical velocity

            for (int n = 1; n <= (int)params.N; ++n) {
                // Compute the predicted state for each step in the horizon
                DM state = DM::zeros(n_states, 1);
                
                // Updating the state based on your MATLAB code logic
                state(0) = qk_(0) + n * params.dt * traj(3); // pos X
                state(1) = qk_(1) + n * params.dt * traj(4); // pos Y
                state(2) = params.height;                  // pos Z (constant height assumption)
                state(3) = traj(3);                        // vel X (constant velocity assumption)
                state(4) = traj(4);                        // vel Y
                state(5) = 0.0;                            // vel Z (assuming zero vertical velocity)
                // Assuming Euler angles (heading) and angular velocity are zero
                state(6) = 0.0; // Roll or phi
                state(7) = 0.0; // Pitch or theta
                state(8) = heading; // Yaw or psi
                state(9) = 0.0;  // Angular velocities (assuming zeros)
                state(10) = 0.0;
                state(11) = 0.0;

                // Place the computed state into the appropriate segment of xd
                //xd.segment((n-1)*n_states, n_states) = state;
                xd(Slice((n-1)*n_states, n*n_states)) = state;
            }
            break;
        }
        // Handle other cases as necessary
        default:
            std::cerr << "Mode '" << params.gait << "' is not defined in the motion planner.\n";
    }

    if (params.gait != 8) {
        xd = repmat(traj, params.N, 1);
    }

    q_des_ = xd; // Only xd needs to be returned
}

void NMPC::force_inequality(const casadi::SX& U, casadi::SX& g_force) {
    // Use SX for symbolic constraint expressions
    casadi::SX g_force_sx = casadi::SX::vertcat(std::vector<casadi::SX>{});

    // g_ub and g_lb as symbolic vectors with appropriate bounds
    casadi::DM g_ub_DM = casadi::DM::vertcat(std::vector<casadi::DM>{});
    casadi::DM g_lb_DM = casadi::DM::vertcat(std::vector<casadi::DM>{});


    for (size_t k = 0; k < params.N; ++k) {
        for (size_t leg = 0; leg < 4; ++leg) {
            size_t base_idx = leg * 3;
            // Directly construct and append symbolic expressions
            g_force_sx = casadi::SX::vertcat({g_force_sx,
                U(base_idx, k) - params.mu * U(base_idx + 2, k),
                U(base_idx + 1, k) - params.mu * U(base_idx + 2, k),
                U(base_idx, k) + params.mu * U(base_idx + 2, k),
                U(base_idx + 1, k) + params.mu * U(base_idx + 2, k)
            });
            g_ub_DM = casadi::SX::vertcat({g_ub_DM, 0.0, 0.0, casadi::inf, casadi::inf}); 
            g_lb_DM = casadi::SX::vertcat({g_lb_DM, -casadi::inf, -casadi::inf, 0, 0});
        }
    }


    // Convert symbolic SX expressions to DM for g_ub and g_lb if necessary
    // Note: This might not be needed if keeping them symbolic is preferable
    g_force = g_force_sx; // g_force is already SX and can stay symbolic
    g_ub_ = g_ub_DM;
    g_lb_ = g_lb_DM;
}

void NMPC::force_inequality_bounds() {

    // g_ub and g_lb as DM vectors with appropriate bounds
    casadi::DM g_ub_sx;
    casadi::DM g_lb_sx;

    for (size_t k = 0; k < params.N; ++k) {
        for (size_t leg = 0; leg < 4; ++leg) {
            // Directly construct and append symbolic expressions
            g_ub_sx = casadi::SX::vertcat({g_ub_sx, 0.0, 0.0, casadi::inf, casadi::inf}); 
            g_lb_sx = casadi::SX::vertcat({g_lb_sx, -casadi::inf, -casadi::inf, 0, 0});
        }
    }
    g_ub_ = g_ub_sx;
    g_lb_ = g_lb_sx;
}

void NMPC::writeMatrixToFile(const casadi::SX& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Print matrix dimensions
        // file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << std::endl;
        // Print matrix content preserving original shape
        for (int i = 0; i < matrix.size1(); ++i) {
            for (int j = 0; j < matrix.size2(); ++j) {
                file << matrix(i, j);
                if (j < matrix.size2() - 1) file << ", ";
            }
            file << std::endl;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

void NMPC::writeMatrixToFile(const casadi::DM& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Print matrix dimensions
        // file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << std::endl;
        // Print matrix content preserving original shape
        for (int i = 0; i < matrix.size1(); ++i) {
            for (int j = 0; j < matrix.size2(); ++j) {
                file << matrix(i, j);
                if (j < matrix.size2() - 1) file << ", ";
            }
            file << std::endl;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

void NMPC::writeMatrixToFile(const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::trunc);  // Open for output and truncate to clear existing content
    if (file.is_open()) {
        // Helper lambda to write a single matrix
        auto writeMatrix = [&](const casadi::DM& matrix) {
            file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << "\n" << std::endl;
            for (int i = 0; i < matrix.size1(); ++i) {
                for (int j = 0; j < matrix.size2(); ++j) {
                    file << matrix(i, j);
                    if (j < matrix.size2() - 1) file << ", ";
                }
                file << std::endl;
            }
        };

        // Write each matrix with a header
        file << "params.pf:" << std::endl;
        writeMatrix(params.pf);
        file << "\nparams.contact:" << std::endl;
        writeMatrix(params.contact);
        file << "\nparams.desVel:" << std::endl;
        writeMatrix(params.desVel);
        
        // Additionally write the size_t variable
        file << "\nparams.stepcnt: " << params.stepcnt << std::endl;
        file << "\n\nparams.current_time: " << params.current_time << std::endl;
        file << "\n\nparams.mpc_iter: " << params.mpc_iter << std::endl;

        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

void NMPC::equalityconstraints(const SX& X, const SX& U, const SX& P, SX& obj, SX& g_eq, const Function& f) {

    std::string basePath = "tmp/";
    int n_states = params.n_states;
    double dt = params.dt;

    // Initial condition constraints
    SX st = X(Slice(0, n_states), 0);
    g_eq = st - P(Slice(0, n_states));

    for (int k = 0; k < (int)params.N; ++k) {
        SX st = X(Slice(0, n_states), k);
        SX con = U(Slice(0, n_states), k);

        SX des_st = P(Slice(k * n_states, (k + 1) * n_states)); // Desired state

        // writeMatrixToFile(st, basePath + "st.txt");
        // writeMatrixToFile(des_st, basePath + "des_st.txt");
        // writeMatrixToFile(params.Q, basePath + "params_Q.txt");
        // writeMatrixToFile(con, basePath + "con.txt");
        // writeMatrixToFile(params.R, basePath + "params_R.txt");
        // Calculate objective
        SX st_err = st - des_st;
        // writeMatrixToFile(st_err, basePath + "st_err.txt");
        SX stage_cost = mtimes(mtimes(st_err.T(), params.Q), st_err);
        SX terminal_cost = mtimes(mtimes(con.T(), params.R), con);

        // writeMatrixToFile(stage_cost, basePath + "stage_cost.txt");
        // writeMatrixToFile(terminal_cost, basePath + "terminal_cost.txt");

        obj = obj + stage_cost + terminal_cost;

        // Calculate next state prediction
        SX st_next = X(Slice(0, n_states), k + 1);

        std::map<std::string, SX> arg;
        arg["st"] = st;  // Example random inputs
        arg["con"] = con;

        std::map<std::string, SX> result = f(arg);
        SX f_value = result.at("rhs");

        // writeMatrixToFile(st_next, basePath + "st_next.txt");
        // writeMatrixToFile(f_value, basePath + "f_value.txt");

        // Compute constraints
        SX st_next_euler = st + dt * f_value;
        g_eq = SX::vertcat({g_eq, st_next - st_next_euler});

        // writeMatrixToFile(st_next_euler, basePath + "st_next_euler.txt");
        // writeMatrixToFile(g_eq, basePath + "g_eq.txt");
    }
}

void NMPC::set_constraint_bounds() {
    // Calculate the start index for setting lower and upper bounds

    int startIndex = params.n_states * (params.N + 1); // Adjust for zero-based indexing

    // Resize the lbg and ubg if necessary
    if (args.lbg.is_empty() || args.lbg.size1() < startIndex + g_lb_.size1()) {
        args.lbg = DM::zeros(startIndex + g_lb_.size1(), 1); // Adjust total size if needed
    }
    if (args.ubg.is_empty() || args.ubg.size1() < startIndex + g_ub_.size1()) {
        args.ubg = DM::zeros(startIndex + g_ub_.size1(), 1); // Adjust total size if needed
    }

    // Set values for lower bounds
    for (int i = 0; i < g_lb_.size1(); ++i) {
        args.lbg(startIndex + i) = g_lb_(i);
    }

    // Set values for upper bounds
    for (int i = 0; i < g_ub_.size1(); ++i) {
        args.ubg(startIndex + i) = g_ub_(i);
    }
}

void NMPC::set_statebounds() {
    int n_states = params.n_states;  // Number of states
    int n_inputs = params.n_inputs;  // Number of inputs
    int N = params.N;
    double zmax = 150;

    // Initialize LB and UB
    args.lbx = DM::zeros(n_states*(N+1) + n_inputs*N, 1);
    args.ubx = DM::zeros(n_states*(N+1) + n_inputs*N, 1);
    args.lbx(Slice()) = -DM::inf();
    args.ubx(Slice()) = DM::inf();

    for (int k = 0; k < N; ++k) {
        for (int leg = 0; leg < 4; ++leg) {
            int index = n_states*(N+1) + n_inputs*k + 3*leg + 2; // Adjusted for 0-based indexing
            args.lbx(index) = 0;
            args.ubx(index) = params.contact(leg) * zmax; // Assuming contact is std::vector<double>
        }
    }
}

void NMPC::shift() {
    // Prepare input for the function f
    std::map<std::string, DM> arg;
    arg["st"] = qk_;  // Current state
    arg["con"] = ukAll_(0, Slice()).T();  // First control input, transposed to match MATLAB

    // Call the function f
    std::map<std::string, DM> result = f_dyn_(arg);
    DM f_value = result["rhs"];  // Assuming the output is named 'rhs'

    // Update state
    qk_ = qk_ + DM(params.dt) * f_value;

    // Increment time
    params.current_time += params.dt;

    DM prev_ukAll = ukAll_;

    // Shift control inputs
    DM u_new = vertcat(prev_ukAll(Slice(1, prev_ukAll.size1()), Slice()),  // Skip the first control set
                       prev_ukAll(Slice(prev_ukAll.size1()-1, prev_ukAll.size1()), Slice()));  // Duplicate the last control set
    ukAll_ = u_new;
}

void NMPC::code_gen() {

    const double dt = params.dt; // Sample time
    const double sim_time = params.sim_length; // Simulation time
    const float sim_length = sim_time/dt; // Simulation length
    const int n_states = params.n_states; // Number of states
    const int n_inputs = params.n_inputs; // Number of inputs
    const int n_controls = params.n_controls; // Number of inputs
    const int gait = params.gait; // Gait type
    const double G = params.G; // Gravity
    const double mass = params.mass; // Mass

    init(); // Initialize the variables

    // Symbolic variables 
    casadi::SXVector states_vec = {
        SX::sym("Xx"), SX::sym("y"), SX::sym("z"),
        SX::sym("dx"), SX::sym("dy"), SX::sym("dz"),
        SX::sym("phi"), SX::sym("theta"), SX::sym("psi"),
        SX::sym("dphi"), SX::sym("dtheta"), SX::sym("dpsi")
    };

    // Define controls similarly
    casadi::SXVector controls_vec = {
        SX::sym("u1x"), SX::sym("u1y"), SX::sym("u1z"),
        SX::sym("u2x"), SX::sym("u2y"), SX::sym("u2z"),
        SX::sym("u3x"), SX::sym("u3y"), SX::sym("u3z"),
        SX::sym("u4x"), SX::sym("u4y"), SX::sym("u4z")
    };

    casadi::SX states = casadi::SX::vertcat(states_vec);
    casadi::SX controls = casadi::SX::vertcat(controls_vec);

    SX U = SX::sym("U", n_controls, params.N);
    SX P = SX::sym("P", n_states + n_states * params.N);
    SX X = SX::sym("X", n_states, params.N + 1);
    SX OPT_variables = vertcat(reshape(X, n_states * (params.N + 1), 1), reshape(U, n_inputs * params.N, 1));


    DM qk = DM::zeros(n_states);
    qk(2) = 0.08; // Assuming z position is constant or predefined

    double t0 = 0.0;
    DM xx(n_states, sim_length);
    xx(Slice(0, n_states), 0) = qk;  

    DM u0 = DM::zeros(params.N, params.n_inputs); // Control inputs for N horizon steps, not sim_steps

    // States decision variables initialization
    DM X0 = repmat(qk, 1, params.N+1);
    DM xd = repmat(qk, params.N, 1);

    double current_time = 0.0;
    size_t mpciter = 0;

    Args args;
    args.lbg = DM::zeros(n_states * (params.N + 1) + 4*4*params.N, 1);
    args.ubg = DM::zeros(n_states * (params.N + 1) + 4*4*params.N, 1);

    DM uout = DM::zeros(n_controls, sim_length); // Vector to store control outputs
    DM xout = DM::zeros(n_states, sim_length);
    DM feet = DM::zeros(12, sim_length); // Vector to store foot positions
    DM traj = DM::zeros(12, sim_length); // Vector to store desired trajectory
    DM tout = DM::zeros(1, sim_length);


    params.current_time = current_time;
    params.mpc_iter = mpciter;
    planner(); // Assuming gait and other details are handled inside
    args.p = vertcat(qk, xd);

    SX t_rot = casadi::SX::vertcat({states(6), states(7), states(8)});  // Extract phi, theta, and psi into a column vector for Rotz
    SX Rz = Rotz(t_rot);    // Call Rotz with the column vector
    SX z3 = casadi::SX::zeros(3, 3);    // Create a 3x3 zero matrix
    SX i3 = casadi::SX::eye(3);         // Create a 3x3 identity matrix

    DM J = DM::zeros(3,3); // Assuming a paramse-5; J(2,2) = 0.064713601;
    SX omega = casadi::SX::vertcat({states(9), states(10), states(11)}); // Extracting angular velocities (dphi, dtheta, dpsi)
    SX skew_omega = skewSym(omega);
    SX J_sx = SX(J); // Convert J from DM to SX
    SX IwI = mtimes(mtimes(mtimes(Rz, J_sx), Rz.T()), skew_omega) * mtimes(mtimes(Rz, J_sx), Rz.T());

    SX A = SX::zeros(12, 12); // A is a 12x12 matrix given your block structure
    // Fill in the blocks of A
    A(Slice(0, 3), Slice(0, 3)) = z3;
    A(Slice(0, 3), Slice(3, 6)) = i3;
    A(Slice(0, 3), Slice(6, 9)) = z3;
    A(Slice(0, 3), Slice(9, 12)) = z3;

    A(Slice(3, 6), Slice(0, 3)) = z3;
    A(Slice(3, 6), Slice(3, 6)) = z3;
    A(Slice(3, 6), Slice(6, 9)) = z3;
    A(Slice(3, 6), Slice(9, 12)) = z3;

    A(Slice(6, 9), Slice(0, 3)) = z3;
    A(Slice(6, 9), Slice(3, 6)) = z3;
    A(Slice(6, 9), Slice(6, 9)) = z3;
    A(Slice(6, 9), Slice(9, 12)) = Rz.T();

    A(Slice(9, 12), Slice(0, 3)) = z3;
    A(Slice(9, 12), Slice(3, 6)) = z3;
    A(Slice(9, 12), Slice(6, 9)) = z3;
    A(Slice(9, 12), Slice(9, 12)) = z3; //-IwI
    
    SX D = SX::zeros(12,1); // D is a 12x1 vector
    D(Slice(3,6)) = SX::vertcat({SX::zeros(2,1), -G}); // Fill the 4th to 6th rows

    SX p_foot_sx = SX(params.pf);
    SX com_vector = SX::vertcat({states(0), states(1), states(2)});
    SX com_vector_expanded = repmat(com_vector, 1, 4); // Replicate across columns to match p_foot_sx dimensions

    SX rd = p_foot_sx - com_vector_expanded;
    
    SX skew_rd = skewSym(rd);
    SX B_tmp = mtimes(mtimes(Rz, J_sx), Rz.T()) * skew_rd;

    SX B = SX::vertcat({
        SX::zeros(3, 12),
        SX::repmat(1/mass * i3, 1, 4),
        SX::zeros(3, 12),
        mtimes(mtimes(Rz, J_sx), Rz.T()) * skew_rd}); //mtimes(mtimes(Rz, J_sx), Rz.T()) * skew_rd
    SX rhs = mtimes(A, states) + mtimes(B, controls) + D;

    // writeMatrixToFile(rhs, basePath + "rhs.txt");

    Function f_dyn = Function("f_dyn", {states, controls}, {rhs}, {"st", "con"}, {"rhs"});

    // Declare variables to hold the outputs
    SX g_force;

    // Call the function
    force_inequality(U, g_force);

    SX obj = 0, g_eq; 
    equalityconstraints(X, U, P, obj, g_eq, f_dyn);
    SX g = vertcat(g_eq, g_force);
    set_statebounds();
    set_constraint_bounds();
    SXDict nlp_prob = {{"f", obj}, {"x", OPT_variables}, {"g", g}, {"p", P}};

    // Define options for the solver
    Dict opts;

    // For ipopt
    opts["ipopt.max_iter"] = 10;  // Replace Max_mpciter with its actual value
    opts["ipopt.print_level"] = 0;  // Can be changed to 0 or 3 based on verbosity required
    opts["print_time"] = 0;  // Disable printing solver time
    opts["ipopt.acceptable_tol"] = 1e-2;  // Tolerance for stopping criterion
    opts["ipopt.acceptable_obj_change_tol"] = 1e-2;  // Objective change tolerance for stopping

    /// for snopt
    // opts["snopt.Iterations limit"] = 10;
    // opts["snopt.Major Print level"] = 3;
    // opts["snopt.Minor Print level"] = 3;
    // opts["snopt.Derivative option"] = 0;
    // opts["snopt.Major feasibility tolerance"] = 1.0e-3;
    // opts["snopt.Minor feasibility tolerance"] = 1.0e-3;
    // opts["snopt.Major optimality tolerance"] = 1.0e-3;
    // opts["snopt."] = 10;
    // opts["print_level"] = 0;
    // // opts["print_file"] = "snopt.out";
    // opts["major_feasibility_tolerance"] = 1e-2;
    // opts["major_optimality_tolerance"] = 1e-2;

    // Create an NLP solver instance
    Function solver = nlpsol("solver", "ipopt", nlp_prob, opts);
    // Function solver = nlpsol("solver", "snopt", {{"f", obj}, {"x", OPT_variables}, {"g", g}, {"p", P}});

    // Reshape X0 and u0 transposed into column vectors
    DM X0_col = reshape(X0.T(), n_states * (params.N + 1), 1);
    DM u0_col = reshape(u0.T(), n_inputs * params.N, 1);

    // Concatenate reshaped vectors vertically
    DM x0u0= vertcat(X0_col, u0_col);
    args.x0 = x0u0;
    args.p = vertcat(qk, xd);
    // Solve the problem
    DMDict arg = {{"x0", args.x0}, {"lbx", args.lbx}, {"ubx", args.ubx}, {"lbg", args.lbg}, {"ubg", args.ubg}, {
        "p", args.p}};
    DMDict res = solver(arg);


    DM sol = DM::zeros(n_states * (params.N + 1) + n_controls * params.N, 1); // sol.x should be initialized with actual solution data
    sol = res.at("x");


    // =========================== CODE GEN =========================== //

    std::string file_name = "nlp_code";
    // code predix
    std::string prefix_code = fs::current_path().string() + "/";
    casadi::Dict code_gen_opts = casadi::Dict();
    code_gen_opts["cpp"] = false;
    code_gen_opts["with_header"] = false;
    solver.generate_dependencies(file_name + ".c");

    CodeGenerator codegen = CodeGenerator(file_name + "_f" + ".c", code_gen_opts);
    codegen.add(f_dyn);
    codegen.generate();
    std::string prefix_lib = fs::current_path().string() + "/";
    std::string compile_command_solver = "gcc -fPIC -shared -O3 " + 
        prefix_code + file_name + ".c -o " +
        prefix_lib + file_name + ".so -L/usr/local/lib -lipopt";
    
    std::string compile_command_f = "gcc -fPIC -shared -O3 " + 
        prefix_code + file_name + "_f" + ".c -o " +
        prefix_lib + file_name + "_f" + ".so";
    std::cout << compile_command_solver << std::endl;
    std::cout << compile_command_f << std::endl;

    int compile_flag_solver = std::system(compile_command_solver.c_str());
    casadi_assert(compile_flag_solver==0, "Compilation failed");
    std::cout << "Compilation successed!" << std::endl;

    int compile_flag_f = std::system(compile_command_f.c_str());
    casadi_assert(compile_flag_f==0, "Compilation failed");
    std::cout << "Compilation successed!" << std::endl;
}

void NMPC::init()
{
    qk_ = DM::zeros(params.n_states);
    qk_(2) = 0.08; // sitting com height
    q_des_ = DM::zeros(params.n_states, 1);
    ukAll_ = DM::zeros(params.N, params.n_inputs);

    qkAll_ = repmat(qk_, 1, params.N+1); //contains initial state, and solution for all N horizon steps
    q_des_ = repmat(qk_, params.N, 1);   //contains desired state for all N horizon steps

    args.lbg = DM::zeros(params.n_states * (params.N + 1) + 4*4*params.N, 1);
    args.ubg = DM::zeros(params.n_states * (params.N + 1) + 4*4*params.N, 1);

    uout = DM::zeros(params.n_controls, params.sim_length/params.dt);
    xout =  DM::zeros(params.n_states, params.sim_length/params.dt);
    feet =  DM::zeros(12, params.sim_length/params.dt);   // Vector to store foot positions   
    traj =  DM::zeros(12, params.sim_length/params.dt); 
    tout =  DM::zeros(1, params.sim_length/params.dt);
    NMPC_solve_time = DM::zeros(1, params.sim_length/params.dt);  

    std::string file_name = "nlp_code";
    std::string prefix_lib = fs::current_path().string() + "/";
    
    std::string lib_name = prefix_lib + file_name + ".so";
    // Load the nlp solver
    // Define options for the solver
    Dict opts;
    opts["ipopt.max_iter"] = 10;  // Replace Max_mpciter with its actual value
    opts["ipopt.print_level"] = 0;  // Can be changed to 0 or 3 based on verbosity required
    opts["print_time"] = 0;  // Disable printing solver time
    opts["ipopt.acceptable_tol"] = 1e-2;  // Tolerance for stopping criterion
    opts["ipopt.acceptable_obj_change_tol"] = 1e-2;  // Objective change tolerance for stopping
    solver_ = casadi::nlpsol("solver","ipopt", lib_name, opts);

    // Load nonlinear dynamics function 
    lib_name = prefix_lib + file_name + "_f" ".so";
    f_dyn_ = external("f_dyn", lib_name);
}

void NMPC::run_NMPC() {

    // call the planner function to update q_des_ for N horizon steps
    planner();

    // writeMatrixToFile(q_des_, basePath + "q_des_.txt");
    // writeMatrixToFile(basePath + "params.txt");

    force_inequality_bounds();     // Updates g_ub_ and g_lb_ with appropriate bounds
    set_statebounds();                     // Updates args.lbx and args.ubx with appropriate bounds
    set_constraint_bounds(); // Updates args.lbg and args.ubg with g_lb_ and g_ub_
    // writeMatrixToFile(g_ub, basePath + "g_ub.txt");
    // writeMatrixToFile(g_lb, basePath + "g_lb.txt");
    // writeMatrixToFile(params, basePath + "params.txt");

    solve_nlp(); // Solve the NMPC problem

}

void NMPC::solve_nlp() {

    DM q0_col = reshape(qkAll_.T(), params.n_states * (params.N + 1), 1);
    DM u0_col = reshape(ukAll_.T(), params.n_inputs * params.N, 1);

    // Concatenate reshaped vectors vertically
    DM x0u0= vertcat(q0_col, u0_col);
    args.x0 = x0u0;
    args.p = vertcat(qk_, q_des_);

    std::map<std::string, casadi::DM> arg, res;
    arg["lbx"] = args.lbx;
    arg["ubx"] =  args.ubx;
    arg["lbg"] =  args.lbg;
    arg["ubg"] =  args.ubg;
    arg["x0"] = args.x0;
    arg["p"] = args.p;

    // writeMatrixToFile(args.lbx, basePath + "args_lbx.txt");
    // writeMatrixToFile(args.ubx, basePath + "args_ubx.txt");
    // writeMatrixToFile(args.lbg, basePath + "args_lbg.txt");
    // writeMatrixToFile(args.ubg, basePath + "args_ubg.txt");
    // writeMatrixToFile(args.x0, basePath + "args_x0.txt");
    // writeMatrixToFile(args.p, basePath + "args_p.txt");
    
    auto start = std::chrono::high_resolution_clock::now();
    res = solver_(arg);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);


    DM sol = DM::zeros(params.n_states * (params.N + 1) + params.n_controls * params.N, 1); // sol.x should be initialized with actual solution data
    sol = res.at("x");


    casadi_int start_index_u = static_cast<casadi_int>(params.n_states * (params.N + 1));
    casadi_int end_index = static_cast<casadi_int>(sol.numel());

    DM ukAll_sol = reshape(sol(Slice(start_index_u, end_index)).T(), params.n_controls, params.N).T();
    DM sol_x_N = reshape(sol(Slice(0, static_cast<casadi_int>(params.n_states * (params.N + 1)))).T(), params.n_states, params.N + 1).T();

    ukAll_ = ukAll_sol;

    // writeMatrixToFile(ukAll_sol, basePath + "ukAll_.txt");
    // writeMatrixToFile(sol_x_N, basePath + "sol_x_N.txt");

    DM uk = ukAll_sol(0, Slice());
    uk_ = uk;
    uout(Slice(), params.mpc_iter) = uk;     // Store current control outputs
    shift();
    qkAll_ = vertcat(sol_x_N(Slice(1, sol_x_N.size1()), Slice()), sol_x_N(Slice(sol_x_N.size1() - 1, sol_x_N.size1()), Slice()));
    
    // Data Logging
    xout(Slice(), params.mpc_iter) = qk_.T(); // Store current states, transposed
    traj(Slice(), params.mpc_iter) = q_des_(casadi::Slice(0, 12), 0).T(); // Store reshaped desired trajectory
    feet(Slice(), params.mpc_iter) = reshape(params.pf, 1, 12); // Store reshaped foot positions
    tout(0, params.mpc_iter) = params.current_time - params.dt;          // Store current time
    NMPC_solve_time(0, params.mpc_iter) = DM(duration.count()); // Store NMPC solve time

    // writeMatrixToFile(qk_, basePath + "qk_.txt");
    // writeMatrixToFile(uk.T(), basePath + "uk_.txt");

    params.mpc_iter++;
    std::cout << "Current time: " << params.current_time << std::endl;

    // params.current_time += params.dt; is already done in shift()
}

void NMPC::logData()
{
    std::string basePath = "tmp/";
    writeMatrixToFile(xout.T(), basePath + "xout.txt");
    writeMatrixToFile(traj.T(), basePath + "traj.txt");
    writeMatrixToFile(uout.T(), basePath + "uout.txt");
    writeMatrixToFile(feet.T(), basePath + "feet.txt");
    writeMatrixToFile(tout.T(), basePath + "tout.txt");
    writeMatrixToFile(NMPC_solve_time, basePath + "NMPC_solve_time.txt");
}