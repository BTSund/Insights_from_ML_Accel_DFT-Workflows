# Deep_Blue_NPH — melting point by solid–liquid coexistence

Estimates the melting temperature of diamond-phase Si at a fixed pressure with
the two-phase coexistence method. Half of the supercell is melted, and the
solid–liquid system then evolves at constant enthalpy and pressure (NPH). The
temperature drifts towards T_m as the interface moves, and the plateau
temperature during NPH is the melting-point estimate.

## Contents

| File | Purpose |
|---|---|
| `data.lammps` | 5184-atom cubic-diamond Si supercell (18 × 6 × 6 conventional cells, about 98 × 33 × 33 Å), elongated along x. |
| `lmp.in` | Four-stage coexistence run (below). |
| `q.py` | Computes the Lechner–Dellago averaged Steinhardt order parameter q̄₆ for each atom, then plots its profile along the longest box axis. The plot separates the crystal (high q̄₆) from the liquid (low q̄₆). |

## Simulation stages (`lmp.in`)

Timestep is 0.5 fs. Pressure is 49346.16 atm (5 GPa), applied anisotropically.

| # | Ensemble | Conditions | Steps (time) | Output |
|---|---|---|---|---|
| 1 | NPT | 1500 K | 10 000 (5 ps) | `Thermalization_1.data` |
| 2 | NVT on x < 48.88 Å only, rest frozen | 3000 K | 10 000 (5 ps) | `Half_melted.data` |
| 3 | NPT | 1500 K | 10 000 (5 ps) | `Thermilization_2.data` |
| 4 | NPH | 5 GPa | 200 000 (100 ps) | `final_structure.data` |

A trajectory (`traj.lammpstrj`, every 100 steps) and restart files (every 50
steps) are written throughout.

## Running

```bash
cd Deep_Blue_NPH
mpirun -np <N> /path/to/lmp_mpi_chimes -i lmp.in > out.lammps
```

Before running, point `pair_coeff` at your `params.txt.reduced` (see the
top-level README).

### Changing conditions

- **Pressure:** change `49346.16` (atm) on both `fix npt` lines and the
  `fix nph` line. 1 GPa = 9869.23 atm.
- **Initial temperature guess:** change `1500.0` in the `velocity` and
  `fix npt` lines. Start near the expected T_m. If one phase takes over the
  whole box during NPH, restart with a better guess.
- **Different supercell:** set the `region 1 block 0 <x_half> ...` bound to
  about half the new box length along x, so that only half is melted.

## Analysis

### Melting temperature

Average the `Temp` column of the Stage 4 (NPH) thermo output after it
plateaus, while both phases are still present.

### Interface profile (`q.py`)

```bash
python q.py
```

The script looks for its input in this order: `New_final.data`, then
`last.lammpstrj`, then `traj.lammpstrj` (using the last frame). Because
`lmp.in` writes `final_structure.data`, by default it reads the last frame of
`traj.lammpstrj`. To analyse a specific data file, copy or link it to
`New_final.data`.

Output: `q6_order_parameter.png`. Defaults are a 2.9 Å neighbour cutoff and 50
bins, set in `plot_q6_profile(r_cutoff=2.9, num_bins=50)`.

The script assumes an orthogonal box (it ignores tilt factors), which holds for
this setup.
