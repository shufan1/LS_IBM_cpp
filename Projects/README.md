# IMPORTANT

The `Projects/Example_Project` directory only exists to show how to
structure an external (i.e., local) project directory for a simulation.
Do not create real project directories inside this repo if you plan to
push branches to GitHub -- simulation config/input/output belongs
outside the solver's own git history (it accumulates data files that
have no business being version-controlled alongside the source code).

Create real project directories as siblings of this repo, e.g.:

```
/home/groups/ibattiat/sxia/LS_IBM/LS_IBM_sxia/
    LS_IBM_cpp/              <- this repo: solver source + one shared build
    Runs/
        fracture_rough_N120/ <- copy of Example_Project's layout
        fracture_rough_N480/
```

Each project directory holds that simulation's config/input/output and
calls the one shared binary built in `LS_IBM_cpp/bin/` -- see
`Example_Project/scripts/config.sh` for how the path is wired up.

# Notes

- `config/`  -- run parameters for this simulation (grid size, geometry
  case, physical parameters).
- `input/`   -- fixed geometry / IBM coefficients this run needs. For
  Milestone 1/2 (frozen geometry), this is the reference data exported
  from `beta_AD_correctu/RunADRE_no_LS.m`.
- `output/`  -- checkpoints and results this run produces.
- `scripts/` -- `config.sh` defines the paths; `run.sbatch` submits the
  job, pointing at the shared `LS_IBM_cpp/bin/` executable with this
  project's config/input/output paths.
