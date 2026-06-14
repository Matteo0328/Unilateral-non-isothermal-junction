# 1D Semiconductor Device Simulator (Scharfetter-Gummel)

## Description
This MATLAB code performs a one-dimensional simulation of charge transport within a semiconductor device. It solves the coupled system of the Poisson equation, the Continuity equation, and (optionally) the Energy Balance equation using the Gummel iterative method.

The simulator supports two physical regimes:
* **Isothermal Model (Drift-Diffusion):** Assumes a constant lattice and electron temperature.
* **Non-Isothermal Model (Hydrodynamic):** Solves for energy transport, accounting for electron heating caused by high electric fields.

---

## Prerequisites
To run the script correctly, ensure your working directory contains the following auxiliary files called by the main code:
* `data.m`: Script containing physical parameters, spatial mesh definition, voltage sweeps, and doping parameters.
* **Scharfetter-Gummel solver functions:**
    * `continuity_n.m`
    * `electron_current_sg.m`
    * `thermal_volt.m`
    * `energy_flux_sg.m`
* `exportPDF.m` (required only if you intend to export plots as PDFs).

---

## Configuration in `data.m`
Before running the simulation, you can define the physical properties, device geometry, and bias sweep inside the `data.m` script. 

The device features a smooth **n+-n-n+** doping profile. The transitions between regions are modeled using hyperbolic tangent functions, defined formally as:
$S_{L,R} = 0.5 \left( 1 \pm \tanh\left(\frac{x - L_C}{\sigma}\right) \right)$

Key parameters you can modify include:

| Category | Parameter | Description | Default Values |
| :--- | :--- | :--- | :--- |
| **Geometry** | `small_device` | Toggles between device sizes. `1` for small, `0` for large. | `0` (Large device) |
| | `L` | Total length of the device. | 3 $\mu$m (Large) / 1 $\mu$m (Small) |
| | `LC` | Length of the highly doped contact regions. | 0.45 $\mu$m (Large) / 0.15 $\mu$m (Small) |
| | `sigma` | Smoothing parameter for the doping profile transition. | 30 nm (Large) / 10 nm (Small) |
| **Doping** | `ND_low` | Doping concentration in the intrinsic/low-doped channel. | 1e22 m^-3 (1e16 cm^-3) |
| | `ND_high` | Doping concentration in the n+ contact regions. | 1e24 m^-3 (1e18 cm^-3) |
| **Physics** | `tau_e` | Energy relaxation time (crucial for hydrodynamic modeling). | 0.2e-12 s |
| | `vsat` | Electron saturation velocity. | 1e5 m/s |
| | `mun_0` | Low-field electron mobility. | 0.10 m^2/(V·s) |
| **Simulation**| `N` | Total number of spatial mesh nodes. | 1000 |
| | `Va_list` | Array of applied voltage biases to step through during simulation. | `[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]` V |

---

## Flag Configuration (Main Script)
At the beginning of the main file (in Section `1. Initialization and Configuration`), you will find a series of boolean variables (`0` = False, `1` = True). Configure them to customize the simulation logic and output:

* `ISOTHERMAL`: Set to `1` to run the constant temperature model, or `0` for the full hydrodynamic model.
* `MOBILITY_T`: Set to `1` to enable electron-temperature-dependent mobility.
* `MOBILITY_E`: Set to `1` to enable local-electric-field-dependent mobility. 
  > **Note:** Activate only one of the two mobility models at a time.
* `PLOTBIAS`: Set to `1` to display full plots (density, potential, currents, energy terms) at each bias step.
* `PLOTTEMP`: Set to `1` to display the temperature profile plot.
* `PRINT`: Set to `1` to export plots to PDF (requires an existing folder named `Es_1`).

---

## How to Run the Code (Simulation Workflow)
Non-isothermal systems are highly non-linear and require an accurate initial guess to converge. Therefore, the recommended workflow is divided into two phases:

### Phase 1: Isothermal Baseline Calculation
1. Open the main script and set `ISOTHERMAL = 1;` (ensure advanced mobility flags are set to `0`).
2. Run the script.
3. The simulator will compute the solutions in the drift-diffusion regime and save a fundamental file named `soluzione_iso.mat` in your current folder.

### Phase 2: Hydrodynamic Simulation (Non-Isothermal)
1. Modify the script by setting `ISOTHERMAL = 0;`. The code will automatically set `initial_guess = 2`, instructing the solver to use the previously calculated isothermal solutions as a starting point.
2. Enable mobility corrections if desired (`MOBILITY_T = 1` or `MOBILITY_E = 1`).
3. Run the script.
4. At each bias step defined in `Va_list`, the console will print the progress of the Gummel iterations and the final diagnostics (residual errors, energy balance violations, drift velocities).

---

## Output Structure and Post-Processing
The program will automatically create three folders (if they do not already exist) to organize and save the results as `.mat` data structures:
* `Variazione tau/`
* `Variazione Nint/`
* `Variazione mu/`

**Post-Processing Scripts:** Within these subfolders, you will find additional MATLAB scripts dedicated to post-processing the exported data. Run these to generate custom tables and comparative plots required for further analysis.

At the end of the main simulation, the code will automatically generate a baseline plot displaying the **I-V Characteristic** (Current-Voltage). This allows you to visually compare the difference between the isothermal and non-isothermal model predictions.