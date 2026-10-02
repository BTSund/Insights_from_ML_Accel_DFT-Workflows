# Deep_Blue_TI — Gibbs free energies by thermodynamic integration

Computes the Gibbs free energy of each Si phase as a function of pressure at
fixed temperature, using Frenkel–Ladd thermodynamic integration between the
ChIMES crystal and an Einstein crystal (LAMMPS `fix ti/spring`). The script
then finds where the lowest-G phase changes (solid–solid transition pressures)
and propagates replica-to-replica scatter into an uncertainty on each
transition pressure.

## Contents

| File | Purpose |
|---|---|
| `configs/CON_<phase>` | Starting LAMMPS data files, one per phase. |
| `starterpack/lmp.in` | **Stage 1:** NPT at `set_temp` / `set_press` for 5000 steps (1 fs). Writes the box dimensions, averaged over the last 2500 steps, to `averages.txt`. |
| `starterpack/TI.in` | **Stage 2:** NVT (Langevin) in the averaged box with `fix ti/spring` (k = 900 kcal/mol/Å², 5000-step switch, 10000-step equilibration). Writes `forward_fl.dat` (crystal → Einstein) and `backward_fl.dat` (Einstein → crystal). |
| `run_TI.sh` | Builds the run tree for every temperature × replica × phase × pressure and fills in the placeholders. |
| `sbatch.cmd` | SLURM job template. Runs 5 replicas in parallel (48 MPI ranks each) and loops over the pressures for one phase. |
| `run_sbatch.sh` | Submits one copy of `sbatch.cmd` per phase per temperature. |
| `Spline_uncertainty.py` | Analysis: ΔG(P) per phase, transition pressures, uncertainties, and plot. |
| `stdoutmsg` | Leftover SLURM output from an earlier run. Safe to delete. |

### Placeholders filled by the scripts

| Token | File | Replaced with | By |
|---|---|---|---|
| `set_temp` | `lmp.in`, `TI.in` | Temperature in K | `run_TI.sh` |
| `set_press` | `lmp.in` | Pressure in atm (GPa × 9869) | `run_TI.sh` |
| `rand_1`, `rand_2` | `lmp.in`, `TI.in` | Velocity seeds that differ per replica | `run_TI.sh` |
| `set_lx` … `set_yz` | `TI.in` | Averaged box from `averages.txt` | `sbatch.cmd` |

## Workflow

### 1. Build the directory tree

```bash
cd Deep_Blue_TI
bash run_TI.sh
```

This creates:

```
TI_<T>/                      # T = 300, 650, 1000 K
├── starterpack/             # temperature already substituted
└── rep_<0..4>/              # 5 independent replicas
    └── CON_<phase>/
        └── Press_<PPP>/     # e.g. Press_010 = 10 GPa
            ├── data.lammps
            ├── lmp.in
            └── TI.in
```

Edit `temps=(...)` at the top of `run_TI.sh` to change temperatures. Edit the
`case` block to change the pressures sampled for each phase (in GPa).

### 2. Submit

Before the first submission, edit `sbatch.cmd`:

- `#SBATCH -p [queue]`: set your partition.
- `module load ...`: set the modules your LAMMPS build needs.
- The `ibrun ... lmp_mpi_chimes` lines: replace the hard-coded binary path, and
  on non-TACC machines replace `ibrun -n 48 -o $OFFSET` with your launcher
  (e.g. `srun -n 48 --exact`).
- If you change the number of replicas or ranks per replica, update `-N`,
  `--ntasks-per-node`, `{0..4}`, and `OFFSET` together.

Then:

```bash
bash run_sbatch.sh
```

Each job checks for finished pressures and skips them. A pressure counts as
finished when the last line of `backward_fl.dat` contains step 30000, so a
timed-out job can be resubmitted as-is.

### 3. Analyse

Run the analysis from inside one temperature directory:

```bash
cd TI_1000
python ../Spline_uncertainty.py
```

What the script does:

1. For each `rep_*/CON_*/Press_*`, integrates the forward and backward
   `*_fl.dat` files and averages them to get the free-energy difference ΔF to
   the Einstein crystal (eV).
2. Adds PV, using the average P and V from the last 500 thermo lines of the
   first block in `log.lammps`, and divides by the atom count to get G per
   atom.
3. Averages over replicas. Error bars are the standard error of the mean.
4. Subtracts a linear fit of the reference phase (`target_phase`, default
   `CON_opt_HCP`), so curves are ΔG relative to that phase.
5. Interpolates each phase with PCHIP, scans 5–95 GPa for changes in the
   lowest-G phase, and finds each crossing with `brentq`. For each transition
   it prints the analytical uncertainty on the crossing pressure, from linear
   interpolation between the bracketing data points.
6. Saves `gibbs_stability_shaded.png` (ΔG vs. P with ±1 SEM bands).

**Settings to edit for each run** (all hard-coded in the script):

| Setting | Where | Default |
|---|---|---|
| Reference phase | `target_phase` | `CON_opt_HCP` |
| Number of replicas | `n = range(5)` | 5 |
| Expected rows per `*_fl.dat` | `len(f_data) != 5001` | 5001 (must match the `run 5001` in `TI.in`) |
| Pressure scan range | `p_fine = np.linspace(5, 95, 2000)` | 5–95 GPa |
| Plot title temperature | `plt.title(...)` | "1000K" |
| Plot window | `plt.xlim`, `plt.ylim` | 70–95 GPa, ±1.5 kcal/mol |



