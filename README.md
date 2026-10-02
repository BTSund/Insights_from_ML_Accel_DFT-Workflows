# ChIMES Silicon Phase-Diagram Workflows

LAMMPS workflows for mapping the high-pressure phase behaviour of silicon with
[ChIMES](https://github.com/rk-lindsey/chimes_calculator) machine-learned
interatomic potentials. The three folders cover the three calculations needed
to build a P–T phase diagram:

| Folder | What it computes | Main outputs |
|---|---|---|
| [`Example_TI/`](Example_TI/README.md) | Finite-temperature Gibbs free energies by Frenkel–Ladd thermodynamic integration (TI) to an Einstein crystal | ΔG(P) curves per phase, solid–solid transition pressures with uncertainties |
| [`Example_NPH/`](Example_NPH/README.md) | Melting points by solid–liquid two-phase coexistence (NPH) | Coexistence temperature at a given pressure, plus a q̄₆ order-parameter profile across the interface |

All simulations use `units real` (energies in kcal/mol, pressures in atm,
lengths in Å, time in fs).

## Requirements

**Simulation**

- LAMMPS built with the ChIMES pair style (`pair_style chimesFF`). See the
  [ChIMES calculator](https://github.com/rk-lindsey/chimes_calculator) for
  build instructions.
- A ChIMES parameter file for Si. The inputs expect files named `params.txt`
  (E–V) and `params.txt.reduced` (TI, NPH). **These files are not included in
  this repository**, so supply your own (see [Parameter files](#parameter-files)).
- For TI: a SLURM cluster. `sbatch.cmd` was written for TACC Stampede3
  (`ibrun`, Intel MPI) and needs editing for other machines.

**Analysis**

Python 3.9+ with the packages in [`requirements.txt`](requirements.txt):

```bash
pip install -r requirements.txt
```

## Parameter files

Every LAMMPS input loads the potential with a relative path, e.g.

```
pair_coeff * * ../../../params.txt.reduced
```

Performed using lammps 'stable_29Aug2024_update1' version

These paths match the directory layout the workflows were first run in and
differ between folders. The most reliable fix is to point them at an absolute
path before running.


## Repository layout

```
.
├── Example_TI/
│   ├── configs/            # LAMMPS data files, one per phase, diamond example provided
│   ├── starterpack/        # lmp.in (NPT box averaging), TI.in (Frenkel–Ladd TI)
│   ├── run_TI.sh           # builds the T / replica / phase / pressure directory tree
│   ├── run_sbatch.sh       # submits one SLURM job per phase per temperature
│   ├── sbatch.cmd          # SLURM job template
│   └── Spline_uncertainty.py  # free-energy analysis, transitions, plots
├── Example_NPH/
│   ├── data.lammps         # 5184-atom diamond Si supercell (18×6×6 conventional cells)
│   ├── lmp.in              # equilibrate → half-melt → NPH coexistence
│   └── q.py                # q̄₆ order-parameter profile along the long axis
└── requirements.txt
```

## Citation

If you use these workflows, please cite the associated manuscript (reference to
be added) and the ChIMES papers listed in the
[ChIMES calculator](https://github.com/rk-lindsey/chimes_calculator) repository.

## License

No license has been chosen yet. Until one is added, all rights are reserved by
the authors.
