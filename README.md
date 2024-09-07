
# SRB_NMPC Project

## Overview

This project consists of a main project for the SRB_NMPC and a nested project located in the `SRB_NMPC_Casadi_Formulation` folder. The nested project is used for a one-time symbolic MPC setup and code generation. Once the solver is exported to shared library objects, it is utilized by the main project to run the NMPC problem in real-time.

## Environment Variables

Before proceeding, ensure the following environment variables are set:

- **LOCAL_INSTALL**: Points to where `raisimLib`, `raisimOgre`, and `ogre` build files are located.
- **WORKSPACE**: Points to the common workspace folder where `raisimLib` is cloned.

Example:
```bash
export LOCAL_INSTALL=/path/to/local/install
export WORKSPACE=/path/to/workspace
```

## Code Generation Steps

The code generation is automated by a bash script named `runCodegen.sh`. Follow these steps to generate the code:

Ensure `runCodegen.sh` has execute permissions:
```bash
chmod +x runCodegen.sh
```

Run the script:
```bash
./runCodegen.sh
```

### What `runCodegen.sh` Does

This script automates the following tasks:

1. **Change Directory**:
   - Navigates to the nested project `SRB_NMPC_Casadi_Formulation`.

2. **Create and Navigate to Build Folder**:
   - Creates the build folder if it doesn’t exist and navigates into it.

3. **CMake Configuration**:
   - Writes CMake configuration files by invoking `cmake ...`

4. **Build Project**:
   - Compiles the project using `make -j16`.

5. **Run Code Generation**:
   - Executes the exported target `./nmpccodegen`, which sets up the symbolic NMPC solver and exports it to shared library objects (.so files). The following files are generated:
     - `nlp_code.so`: Contains the solver.
     - `nlp_code_f.so`: Contains the dynamic function to calculate the next state.

6. **Copy Shared Libraries**:
   - Copies the generated .so files to the main folder’s build directory.

## Running the Main Raisim Project

After completing the code generation steps, follow these steps to run the main Raisim project:

1. **Ensure Scripts Have Execute Permissions**:
   - Set execute permissions for the following scripts:
     ```bash
     chmod +x cmakemake.sh
     chmod +x run_Raisim.sh
     ```

2. **Build the Main Project**:
   - Run the build script:
     ```bash
     ./cmakemake.sh
     ```

3. **Run the Simulation**:
   - Execute the Raisim run script:
     ```bash
     ./run_Raisim.sh
     ```

## Tuning NMPC Gains

The gains of the NMPC can be tuned by modifying the P, Q, and other gain matrices in the `/SRB_NMPC_Casadi_Formulation/include/nmpc.h` file. After making changes to the solver configuration, including gains or any symbolic elements, re-run `./runCodegen.sh` to regenerate the solver.