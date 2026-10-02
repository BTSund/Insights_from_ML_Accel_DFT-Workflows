import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.spatial import cKDTree 
from matplotlib.ticker import FormatStrFormatter

plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({
    'font.size': 16,          # base size (tick labels, etc.)
    'axes.labelsize': 20,     # x and y axis labels
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'axes.titlesize': 20,
})
def read_lammps_data(filename):
    """Parses a standard LAMMPS data file (write_data format)."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    xlo = xhi = ylo = yhi = zlo = zhi = 0.0
    coords = []
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if "xlo xhi" in line:
            xlo, xhi = float(line.split()[0]), float(line.split()[1])
        elif "ylo yhi" in line:
            ylo, yhi = float(line.split()[0]), float(line.split()[1])
        elif "zlo zhi" in line:
            zlo, zhi = float(line.split()[0]), float(line.split()[1])
        elif line.startswith("Atoms"):
            i += 1
            # Skip blank lines and section comments
            while i < len(lines) and (not lines[i].strip() or lines[i].strip().startswith("#")):
                i += 1
            # Parse atom entries (id type x y z ...)
            while i < len(lines) and lines[i].strip():
                parts = lines[i].split()
                if len(parts) >= 5:
                    coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
                i += 1
        i += 1
        
    box_min = np.array([xlo, ylo, zlo])
    box_max = np.array([xhi, yhi, zhi])
    return box_min, box_max, np.array(coords)

def read_last_frame_trajectory(filename):
    """Parses the final frame of a LAMMPS trajectory file."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    last_frame_start = None
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].startswith("ITEM: TIMESTEP"):
            last_frame_start = i
            break
            
    if last_frame_start is None:
        raise ValueError(f"ITEM: TIMESTEP header not found in {filename}.")
        
    frame_lines = lines[last_frame_start:]
    box_bounds, atoms_data, atoms_header = [], [], None
    
    idx = 0
    while idx < len(frame_lines):
        line = frame_lines[idx].strip()
        if line.startswith("ITEM: BOX BOUNDS"):
            for j in range(1, 4):
                box_bounds.append([float(val) for val in frame_lines[idx + j].split()])
            idx += 4
        elif line.startswith("ITEM: ATOMS"):
            atoms_header = line.split()[2:]
            idx += 1
            while idx < len(frame_lines) and not frame_lines[idx].startswith("ITEM:"):
                atoms_data.append(frame_lines[idx].split())
                idx += 1
        else:
            idx += 1
            
    box_min = np.array([box_bounds[0][0], box_bounds[1][0], box_bounds[2][0]])
    box_max = np.array([box_bounds[0][1], box_bounds[1][1], box_bounds[2][1]])
    
    x_col = "xu" if "xu" in atoms_header else "x"
    y_col = "yu" if "yu" in atoms_header else "y"
    z_col = "zu" if "zu" in atoms_header else "z"
    
    x_idx = atoms_header.index(x_col)
    y_idx = atoms_header.index(y_col)
    z_idx = atoms_header.index(z_col)
    
    coords = []
    for row in atoms_data:
        coords.append([float(row[x_idx]), float(row[y_idx]), float(row[z_idx])])
        
    return box_min, box_max, np.array(coords)

def load_structure():
    """Hierarchical file detection logic: New_final.data -> last.lammpstrj -> traj.lammpstrj"""
    if os.path.exists("New_final.data"):
        print("Found 'New_final.data'. Reading structure...")
        return read_lammps_data("New_final.data")
    elif os.path.exists("last.lammpstrj"):
        print("Found 'last.lammpstrj'. Reading final frame...")
        return read_last_frame_trajectory("last.lammpstrj")
    elif os.path.exists("traj.lammpstrj"):
        print("Found 'traj.lammpstrj'. Reading final frame...")
        return read_last_frame_trajectory("traj.lammpstrj")
    else:
        raise FileNotFoundError(
            "None of the target files ('New_final.data', 'last.lammpstrj', 'traj.lammpstrj') were found."
        )

def plot_q6_profile(r_cutoff=2.9, num_bins=50):
    box_min, box_max, coords = load_structure()
    
    box_lengths = box_max - box_min
    axes = ['X', 'Y', 'Z']
    
    longest_idx = np.argmax(box_lengths)
    longest_axis = axes[longest_idx]
    
    # Wrap coordinates safely into [0, length)
    coords_shifted = (coords - box_min) % box_lengths
    coords_wrapped = np.clip(coords_shifted, 0.0, box_lengths - 1e-12)
    
    # Query neighbor bonds
    tree = cKDTree(coords_wrapped, boxsize=box_lengths)
    pairs = tree.query_pairs(r=r_cutoff, output_type='ndarray')
    
    N = len(coords)
    q6m = np.zeros((N, 13), dtype=complex)
    neighbor_counts = np.zeros(N, dtype=int)
    
    if len(pairs) > 0:
        i_idx, j_idx = pairs[:, 0], pairs[:, 1]
        
        # Minimum image displacement
        dr = coords_wrapped[j_idx] - coords_wrapped[i_idx]
        dr = dr - box_lengths * np.round(dr / box_lengths)
        
        dr_all = np.vstack([dr, -dr])
        i_all = np.concatenate([i_idx, j_idx])
        
        r = np.linalg.norm(dr_all, axis=1)
        valid = (r > 1e-5)
        
        dr_all, i_all, r = dr_all[valid], i_all[valid], r[valid]
        x, y, z = dr_all[:, 0], dr_all[:, 1], dr_all[:, 2]
        
        theta = np.arccos(np.clip(z / r, -1.0, 1.0))
        phi = np.arctan2(y, x) % (2 * np.pi)
        
        # Spherical harmonics q6m calculation
        for m_idx, m in enumerate(range(-6, 7)):
            Ylm = sp.sph_harm(m, 6, phi, theta)
            np.add.at(q6m[:, m_idx], i_all, Ylm)
            
        np.add.at(neighbor_counts, i_all, 1)
        
    # Local atom q6m average
    valid_atoms = neighbor_counts > 0
    q6m_avg = np.zeros_like(q6m)
    q6m_avg[valid_atoms] = q6m[valid_atoms] / neighbor_counts[valid_atoms, None]
    
    # Lechner-Dellago neighbor-averaged q6m
    bar_q6m = np.zeros_like(q6m)
    bar_counts = np.zeros(N, dtype=int)
    
    bar_q6m[valid_atoms] += q6m_avg[valid_atoms]
    bar_counts[valid_atoms] += 1
    
    if len(pairs) > 0:
        j_all = np.concatenate([j_idx, i_idx])[valid]
        np.add.at(bar_q6m, i_all, q6m_avg[j_all])
        np.add.at(bar_counts, i_all, 1)
        
    valid_bar = bar_counts > 0
    bar_q6m[valid_bar] /= bar_counts[valid_bar, None]
    
    # Calculate scalar bar_q6
    bar_q6 = np.sqrt((4 * np.pi / 13.0) * np.sum(np.abs(bar_q6m)**2, axis=1))
    
    # Binning along longest axis
    L_longest = box_lengths[longest_idx]
    lo_longest = box_min[longest_idx]
    atom_positions = lo_longest + coords_wrapped[:, longest_idx]
    
    bin_edges = np.linspace(lo_longest, lo_longest + L_longest, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    q6_binned = np.zeros(num_bins)
    
    for k in range(num_bins):
        mask = (atom_positions >= bin_edges[k]) & (atom_positions < bin_edges[k+1])
        q6_binned[k] = np.mean(bar_q6[mask]) if np.any(mask) else np.nan

    # Plot profile
    plt.figure(figsize=(8, 5))
    plt.plot(bin_centers, q6_binned, color='#2ca02c', linewidth=2, marker='o', markersize=3)
    
    plt.xlabel(f'{longest_axis} Position (Å)')
    plt.ylabel('Steinhardt Order Parameter $\\bar{q}_6$')
    plt.xlim(lo_longest, lo_longest + L_longest)
    # plt.ylim(0.0, 0.6)
    plt.tight_layout()
    
    plt.savefig('q6_order_parameter.png', dpi=300)
    print(f"Calculated q6 along {longest_axis}-axis. Saved plot to q6_order_parameter.png")

if __name__ == "__main__":
    plot_q6_profile()
