import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from collections import deque
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq # Added for finding intersections

# --- Helper Functions (Logic remains the same) ---

import numpy as np

def UQ(exact_p, p1_lower, p1_upper, p2_lower, p2_upper):
    """
    Compares Analytical Error Propagation vs Monte Carlo for transition pressure.
    Points are tuples: (x, y, sigma)
    """
    # Unpack for Phase A (current_lowest)
    xA1, yA1, sA1 = p1_lower
    xA2, yA2, sA2 = p1_upper
    
    # Unpack for Phase B (best_phase)
    xB1, yB1, sB1 = p2_lower
    xB2, yB2, sB2 = p2_upper

    dA = xA2 - xA1
    dB = xB2 - xB1

    # --- Analytical Propagation ---
    mA = (yA2 - yA1) / dA
    bA = (xA2 * yA1 - xA1 * yA2) / dA
    mB = (yB2 - yB1) / dB
    bB = (xB2 * yB1 - xB1 * yB2) / dB

    N = bB - bA
    D = mA - mB
    xstar = N / D

    # Partials
    dN_dyA1, dN_dyA2 = -xA2 / dA,  xA1 / dA
    dN_dyB1, dN_dyB2 =  xB2 / dB, -xB1 / dB
    dD_dyA1, dD_dyA2 = -1.0 / dA,  1.0 / dA
    dD_dyB1, dD_dyB2 =  1.0 / dB, -1.0 / dB

    def dx_dy(dN, dD):
        return (dN * D - N * dD) / (D**2)

    dxA1 = dx_dy(dN_dyA1, dD_dyA1)
    dxA2 = dx_dy(dN_dyA2, dD_dyA2)
    dxB1 = dx_dy(dN_dyB1, dD_dyB1)
    dxB2 = dx_dy(dN_dyB2, dD_dyB2)

    sigma_x = np.sqrt(dxA1**2 * sA1**2 + dxA2**2 * sA2**2 +
                     dxB1**2 * sB1**2 + dxB2**2 * sB2**2)

    # --- Monte Carlo Check ---
    Nmc = 200000
    ys = np.random.normal([yA1, yA2, yB1, yB2],
                          [sA1, sA2, sB1, sB2],
                          size=(Nmc, 4))
    
    # Vectorized MC for speed
    mA_mc = (ys[:, 1] - ys[:, 0]) / dA
    bA_mc = (xA2 * ys[:, 0] - xA1 * ys[:, 1]) / dA
    mB_mc = (ys[:, 3] - ys[:, 2]) / dB
    bB_mc = (xB2 * ys[:, 2] - xB1 * ys[:, 3]) / dB
    
    D_mc = mA_mc - mB_mc
    x_mc = (bB_mc - bA_mc) / D_mc
    
    # Filter for stability
    xmc_valid = x_mc[np.isfinite(x_mc)]

    # print(f"\n--- UQ Comparison at {exact_p:.4f} GPa ---")
    print(f"Analytical: x* = {xstar:.4f}, σ_x = {sigma_x:.4f}")
    # print(f"Monte Carlo: Mean = {np.mean(xmc_valid):.4f}, Std = {np.std(xmc_valid):.4f}")
    # print(f"Difference in σ: {abs(sigma_x - np.std(xmc_valid)):.2e}")


def get_last_thermo_values(logfile):
    press_deque, vol_deque = deque(maxlen=500), deque(maxlen=500)
    with open(logfile) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "Step" in line and "Temp" in line and "Press" in line and "Volume" in line:
            headers = line.strip().split()
            p_idx, v_idx = headers.index("Press"), headers.index("Volume")
            j = i + 1
            while j < len(lines) and lines[j].strip() and (lines[j].strip()[0].isdigit() or lines[j].startswith('-')):
                vals = lines[j].split()
                press_deque.append(float(vals[p_idx]))
                vol_deque.append(float(vals[v_idx]))
                j += 1
            break
    if not press_deque: raise ValueError(f"No thermo data in {logfile}")
    return sum(press_deque)/len(press_deque), sum(vol_deque)/len(vol_deque)

def get_num_atoms(logfile):
    with open(logfile) as f:
        for line in f:
            if "reading atoms" in line:
                return int(next(f).strip().split()[0])
    raise ValueError("Atoms count not found")

# --- Data Aggregation ---

n = range(5)
replicate_dirs = []
for i in n:
    replicate_dirs.append("rep_"+str(i))
# Group data: data_store[phase][press_folder] = {'P': [vals], 'G': [vals]}
data_store = {}

for rep in replicate_dirs:
    print(f"Scanning {rep}...")
    for f_path in glob.glob(f"{rep}/CON*/Press_*/forward_fl.dat"):
        press_dir = os.path.dirname(f_path)
        press_name = os.path.basename(press_dir)
        phase = os.path.basename(os.path.dirname(press_dir))
        b_path = os.path.join(press_dir, "backward_fl.dat")
        log_f = os.path.join(press_dir, "log.lammps")

        if not os.path.exists(b_path) or not os.path.exists(log_f): continue

        try:
            # TI Integrand Calculation
            def get_df(path):
                d = np.loadtxt(path, skiprows=2)
                TI = (-d[:,1] + d[:,2]) * d[:,3] * 2 / 23.062730
                return np.trapz(TI, x=d[:,3])
            f_data = np.loadtxt(f_path, skiprows=2)
            b_data = np.loadtxt(b_path, skiprows=2)
            if len(f_data) != 5001 or len(b_data) != 5001:
                print("Skipping: Data files do not have 5001 points.")
                continue

            delta_F = (get_df(f_path) - get_df(b_path)) / 2
            p_val, v_val = get_last_thermo_values(log_f)
            n_atoms = get_num_atoms(log_f)
            
            G = (delta_F + (6.32e-7 * p_val * v_val)) / n_atoms

            if phase not in data_store: data_store[phase] = {}
            if press_name not in data_store[phase]: data_store[phase][press_name] = {'P':[], 'G':[]}
            
            data_store[phase][press_name]['P'].append(p_val)
            data_store[phase][press_name]['G'].append(G)
        except Exception as e:
            print(f"Error in {press_dir}: {e}")

# --- Statistics & Baseline Subtraction ---

processed = {}
for phase, p_dict in data_store.items():
    p_means, g_means, g_errs = [], [], []
    for p_name in p_dict:
        p_means.append(np.mean(p_dict[p_name]['P']))
        g_means.append(np.mean(p_dict[p_name]['G']))
        # Standard Error of the Mean (SEM)
        g_errs.append(np.std(p_dict[p_name]['G']) / np.sqrt(len(replicate_dirs)))
    
    idx = np.argsort(p_means)
    processed[phase] = {'p': np.array(p_means)[idx], 'g': np.array(g_means)[idx], 'err': np.array(g_errs)[idx]}

target_phase = 'CON_opt_HCP'
if target_phase in processed:
    a, b = np.polyfit(processed[target_phase]['p'], processed[target_phase]['g'], 1)
    for phase in processed:
        processed[phase]['g_adj'] = processed[phase]['g'] - (a * processed[phase]['p'] + b)
else:
    for phase in processed: processed[phase]['g_adj'] = processed[phase]['g']

# --- Plotting with Shaded Uncertainty ---

plt.figure(figsize=(10, 7))
ATM_TO_GPA = 9870
EV_TO_KCAL = 23.06
splines = {}

# 1. Create Splines for all phases
for phase, data in processed.items():
    p_gpa = data['p'] / ATM_TO_GPA
    g_kcal = data['g_adj'] * EV_TO_KCAL
    splines[phase] = PchipInterpolator(p_gpa, g_kcal)

# 2. Determine the overarching pressure range
print(processed.values())
all_p = np.concatenate([data['p'] / ATM_TO_GPA for data in processed.values()])
p_min, p_max = all_p.min(), all_p.max()
p_fine = np.linspace(p_min, p_max, 1000)
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

# ... [Ensure 'processed' dictionary is already populated above this point] ...

# --- Phase Transition Detection (Integrated Fix) ---

ATM_TO_GPA = 9870
EV_TO_KCAL = 23.06
interpolators = {}  # Initialize the dictionary here

# 1. Create Pchip Interpolators for all phases
for phase, data in processed.items():
    p_gpa = data['p'] / ATM_TO_GPA
    g_kcal = data['g_adj'] * EV_TO_KCAL
    # Using Pchip to prevent artificial oscillations
    interpolators[phase] = PchipInterpolator(p_gpa, g_kcal)

# 2. Define the pressure range for scanning
all_p = np.concatenate([data['p'] / ATM_TO_GPA for data in processed.values()])
p_min, p_max = all_p.min(), all_p.max()
p_fine = np.linspace(5, 95, 2000)
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

# ... [Previous setup code] ...

print("\n--- Ground State Transition & Bracketing Data Points + Uncertainty ---")

current_lowest = None

for p in p_fine:
    min_g = float('inf')
    best_phase = None
    
    for phase, interp in interpolators.items():
        if p < interp.x.min() or p > interp.x.max(): continue
        g = interp(p)
        if g < min_g:
            min_g = g
            best_phase = phase
            
    if best_phase is not None and best_phase != current_lowest:
        if current_lowest is not None:
            # 1. Find the exact transition pressure
            def crossover_fn(px):
                return interpolators[current_lowest](px) - interpolators[best_phase](px)
            
            step_size = p_fine[1] - p_fine[0]
            
            try:
                # print(p)
                exact_p = brentq(crossover_fn, p - step_size, p)
                
            except:
                current_lowest = best_phase
                continue
            
            # 2. Extract Bracketing "True Points" (P, G, Err)
            def get_bracketing_points_with_err(phase_name, target_p):
                raw_p = processed[phase_name]['p'] / ATM_TO_GPA
                raw_g = processed[phase_name]['g_adj'] * EV_TO_KCAL
                # Get errors and convert units
                raw_err = processed[phase_name]['err'] * EV_TO_KCAL
                
                below_mask = raw_p < target_p
                above_mask = raw_p > target_p
                
                if not np.any(below_mask) or not np.any(above_mask):
                    return None, None 

                # Extract the single closest points (assuming raw_p is sorted)
                p_below = raw_p[below_mask][-1]
                g_below = raw_g[below_mask][-1]
                err_below = raw_err[below_mask][-1]
                
                p_above = raw_p[above_mask][0]
                g_above = raw_g[above_mask][0]
                err_above = raw_err[above_mask][0]
                
                # Returns tuples of (P, G, Uncertainty)
                return (p_below, g_below, err_below), (p_above, g_above, err_above)

            left_old, right_old = get_bracketing_points_with_err(current_lowest, exact_p)
            left_new, right_new = get_bracketing_points_with_err(best_phase, exact_p)

            if left_old and left_new:
                # print(f"Transition at {exact_p:.4f} GPa")
                print(f"  - [{current_lowest}] to {best_phase} at {exact_p:.4f}")
                # print(f"  - [{best_phase}] Lower: {left_new}, Upper: {right_new}")
                
                UQ(exact_p, left_old, right_old, left_new, right_new)
            
        current_lowest = best_phase


for phase, data in processed.items():
    p_gpa = data['p'] / ATM_TO_GPA
    g_kcal = data['g_adj'] * EV_TO_KCAL
    err_kcal = data['err'] * EV_TO_KCAL
    
    # 1. Calculate the discrete upper and lower bounds
    upper_raw = g_kcal + err_kcal
    lower_raw = g_kcal - err_kcal
    
    # 2. Fit individual splines to the Mean, Max, and Min
    cs_mean = PchipInterpolator(p_gpa, g_kcal)
    cs_upper = PchipInterpolator(p_gpa, upper_raw)
    cs_lower = PchipInterpolator(p_gpa, lower_raw)
    
    # 3. Generate smooth curves
    p_fine = np.linspace(p_gpa.min(), p_gpa.max(), 300)
    g_fine = cs_mean(p_fine)
    upper_fine = cs_upper(p_fine)
    lower_fine = cs_lower(p_fine)

    color = plt.gca()._get_lines.get_next_color()
    
    # Shaded Region (Between the Max and Min splines)
    plt.fill_between(p_fine, lower_fine, upper_fine, color=color, alpha=0.2)
    
    # Central Spline Line
    plt.plot(p_fine, g_fine, '-', color=color, linewidth=2, label=f"{phase}")
    
    # Original Data Points (Means)
    plt.plot(p_gpa, g_kcal, 'o', color=color, markersize=4, alpha=0.5)

plt.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.7)
plt.xlabel("Pressure (GPa)", fontsize=12)
plt.ylabel("$\Delta G$ relative to " + target_phase + " (kcal/mol/atom)", fontsize=12)
plt.title(f"Phase Stability at 1000K (Mean of {len(replicate_dirs)} Replicates)", fontsize=14)
plt.legend(frameon=False)
plt.grid(True, linestyle=':', alpha=0.6)
plt.xlim(70, 95) 
plt.ylim(-1.5, 1.5) 
plt.tight_layout()
plt.savefig("gibbs_stability_shaded.png", dpi=300)
plt.show()

print("Plot with shaded uncertainty saved.")