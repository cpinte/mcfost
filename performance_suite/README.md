# Wall-clock benchmark

Tracks the run time of mcfost per phase, to detect performance regressions and
measure optimisations. The physics benchmarks (Pascucci et al. 2004, Pinte et
al. 2009) are in `test_suite/benchmarks`. This directory is not part of the
pytest test suite; it uses the models of `test_suite/test_data`.

## Files

| File | Content |
|---|---|
| `run_benchmark.sh` | Runs the workloads with one or more binaries and writes a results file |
| `benchmark.sbatch` | Slurm job running `run_benchmark.sh` on a full node |
| `compare_benchmark.py` | Compares a baseline and a candidate binary from one or several results files |
| `submit_ozstar.sh`, `ozstar_remote.sh` | Build commits on the OzSTAR login node and submit the job (run by you, over ssh) |
| `fetch_ozstar.sh` | Shows the queued jobs and copies the results and logs back |
| `ozstar.conf.example` | Settings of the OzSTAR scripts, to copy to `ozstar.conf` (not under version control) |
| `results/` | Results files, one per job |

The times come from `mcfost -timing`, which prints the wall-clock time of each
phase (`TIMING <phase> <seconds> <calls>`). Phases can be nested: `init` and
`step2_setup` include `dust_prop` and `opacity`, `total` includes everything.
`compare_benchmark.py` adds `other`: the part of `total` that no phase covers.

| Phase | Content |
|---|---|
| `setup` | Reading the parameters, data directories (before `init`) |
| `init` | Grid, density, dust properties and opacities of the first step |
| `dust_prop` | Grain properties (Mie, DHS, ...) |
| `opacity` | Opacities per cell |
| `mc_temperature` | Monte Carlo packets of the temperature step |
| `temperature` | Equilibrium temperatures (LTE and nLTE) |
| `temperature_nRE` | Stochastic heating: temperature probabilities and re-emission (`Temp_nRE`, `update_proba_abs_nRE`, `emission_nRE`) |
| `step2_setup` | SED step: wavelength grid, dust properties, opacities and, per wavelength, scattering matrices, dark zone and emission probabilities before the packets |
| `mc_step2` | Monte Carlo packets of the image or SED step |
| `source_fct` | Source functions of the ray tracing (including `calc_Jth`) |
| `ray_tracing` | Ray tracing of the images and SEDs (`dust_map`) |
| `output` | Writing the temperature, SED and image files of the dust transfer |
| `mol`, `atom` | Molecular and atomic line transfer |

## Workloads

Built from the models of the test suite (`test_suite/get_test_data.sh`).

| ID | Model | Run | Main phases |
|---|---|---|---|
| W1 | `ref4.1_hydrostatic`, no hydrostatic equilibrium, 1.28e7 packets | temperature, `-no_SED` | `mc_temperature` |
| W2 | same | 1 um image | `mc_step2`, `source_fct`, `ray_tracing` |
| W3 | `ref3.0`, 100 velocity channels | `-mol -no_SED` | `mol` |
| W4 | `discF_00500` (Phantom, Voronoi mesh) | temperature, `-no_SED` | `init` (tessellation), `mc_temperature` |
| W5 | `ref4.1_PAH`, 10 PAH sizes, 1.28e6 packets | temperature and SED | `temperature_nRE`, `source_fct` |

On an Apple M3 Max with 8 threads: W1 8.5 s, W2 6 s, W3 2.5 s, W4 29 s, W5 6 s.
A pass with 3 repeats of two binaries takes about 10 minutes there.

## Usage

Locally, e.g. a quick check on 4 threads with one repeat of the PAH workload:

    export MCFOST_UTILS=...
    ./run_benchmark.sh -t 4 -b main=../../mcfost_main/src/mcfost -b new=../src/mcfost -n 1 -w W5

On OzSTAR (copy `ozstar.conf.example` to `ozstar.conf` and edit it first):

    ./submit_ozstar.sh main pah_nre      # commits must be pushed to GitHub
    ./fetch_ozstar.sh                    # once the job has finished
    ./compare_benchmark.py results/<file>.tsv

`submit_ozstar.sh` resolves the references locally, copies the benchmark scripts
of this checkout, then on the login node builds each commit once in
`<REMOTE_BASE>/build/<sha>` (`make static=yes`, compute nodes have no compiler),
prints the `ldd` of the candidate (libraries must be reachable from the compute
nodes, see `JOB_SETUP`), installs the test data of the candidate (cached) and
submits `benchmark.sbatch`. Results and logs are in `<REMOTE_BASE>/results` and
`<REMOTE_BASE>/logs`.

The remote commands run in a login shell (`bash -l`), which reads
`~/.bash_profile` (or `~/.profile`), and `~/.bashrc` only if the profile sources
it. The shell is not interactive, so the part of `~/.bashrc` after a guard such
as `[[ $- != *i* ]] && return` is skipped. To check what the scripts will see:

    ssh ozstar 'bash -l -c "module list; which gfortran; echo \$MCFOST_INSTALL"'

If this shows the compiler and libraries used to build mcfost, `BUILD_SETUP` can
be empty. The job is submitted with `--export=ALL`, so it inherits this
environment (e.g. `LD_LIBRARY_PATH`) and `JOB_SETUP` is usually not needed.

On another slurm cluster, by hand from this directory:

    sbatch --account=<account> benchmark.sbatch -b main=<path>/mcfost -b new=<path>/mcfost -n 3

Binaries older than the `-timing` option (October 2026) only get the total
wall time of each workload, measured by the script.

Comparison:

    ./compare_benchmark.py [results/<file>.tsv] [--baseline main --candidate new]

Without a file, the latest one in `results/` is used. The baseline and candidate
default to the first two binaries of the file.

A phase is flagged `SLOWER` when the candidate median is slower by more than
10 % (`--threshold`) and by more than the spread of the repeats; phases shorter
than 0.5 s (`--min-seconds`) are never flagged. The exit status is 1 if any
phase is flagged.

## Rules for comparable timings

- Compare binaries in the same job: the repeats of the binaries are interleaved,
  so they see the same node.
- Same node type, compiler, flags and `OMP_NUM_THREADS` for every run that is
  compared; `compare_benchmark.py` warns when the host, CPU, modules or number of
  threads differ between files.
- No other job on the node: `--exclusive` where allowed (not on OzSTAR), otherwise
  `--cpus-per-task` equal to the number of cores of a node. Threads pinned
  (`OMP_PROC_BIND=close`, `OMP_PLACES=cores`).
- Timings of a binary built from a modified tree are marked `-dirty` in the results.
- Each run has its own working directory, so that cached dust properties or
  Voronoi tessellations are never reused by another binary. The mcfost output of
  the first repeat is kept in `results/<file>_logs/` (not under version control).
